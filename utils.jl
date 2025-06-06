using Gurobi, JuMP, Printf, Oil

const ∞ = Inf

times = Vector{Float64}()
uppers = Vector{Float64}()
lowers = Vector{Float64}()

push!(times, 0.0)
push!(uppers, ∞)
push!(lowers, -∞)
cb_calls = Cint[]  # just for debugging
# Setup callback for storing the upper and lower bounds over time

function clean_bounds!(uppers::Vector{Float64}, lowers::Vector{Float64}; upper_placeholder=1e5, lower_placeholder=-1e5)
    for i in eachindex(uppers)
        # Cap and replace invalid values in upper bounds
        if isnan(uppers[i]) || uppers[i] == Inf || uppers[i] > upper_placeholder
            uppers[i] = upper_placeholder
        elseif uppers[i] == -Inf || uppers[i] < lower_placeholder
            uppers[i] = lower_placeholder
        end
        
        # Cap and replace invalid values in lower bounds
        if isnan(lowers[i]) || lowers[i] == -Inf || lowers[i] < lower_placeholder
            lowers[i] = lower_placeholder
        elseif lowers[i] == Inf || lowers[i] > upper_placeholder
            lowers[i] = upper_placeholder
        end
    end
end


function sanitize_bound(value, placeholder)
    if isnan(value) || isinf(value)
        println("Sanitizing bound: Original = $value, Replaced with $placeholder")
        return placeholder
    end
    return value
end


function store_upperbound(cb_data, cb_where::Cint)
    push!(cb_calls, cb_where)
    if cb_where == GRB_CB_MIP
        objbnd = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBND, objbnd)
        
        current_time = time() - start_time

        # Replace Inf/NaN with placeholders for the upper bound
        upper_bound = sanitize_bound(objbnd[], 1e5)

        push!(times, current_time)
        push!(lowers, lowers[end])  # Keep lower bound unchanged
        push!(uppers, upper_bound)  # Update upper bound with cleaned value

        # println("store_upperbound: Time = $current_time, Upper = $upper_bound, Lower = $(lowers[end])")
    end
end


function store_lowerbound(cb_data, cb_where::Cint)
    push!(cb_calls, cb_where)
    if cb_where == GRB_CB_MIP
        objbst = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, objbst)

        current_time = time() - start_time

        # Replace -Inf/NaN with placeholders for the lower bound
        lower_bound = sanitize_bound(objbst[], -1e5)

        push!(times, current_time)
        push!(lowers, lower_bound)  # Update lower bound with cleaned value
        push!(uppers, uppers[end])  # Keep upper bound unchanged

        # println("store_lowerbound: Time = $current_time, Upper = $(uppers[end]), Lower = $lower_bound")
    end
end

function store_both_bounds_gurobi(cb_data, cb_where::Cint; tol=1e-3)
    if cb_where == GRB_CB_MIP
        objbnd = Ref{Cdouble}()
        objbst = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBND, objbnd)
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, objbst)

        current_time = time() - start_time

        # Replace Inf/-Inf with placeholders
        upper_bound = sanitize_bound(objbnd[], 1e5)
        lower_bound = sanitize_bound(objbst[], -1e5)

        # Only store if bounds differ significantly from the last stored values
        if isempty(times) || 
           abs(upper_bound - uppers[end]) > tol || 
           abs(lower_bound - lowers[end]) > tol

            push!(times, current_time)
            push!(uppers, upper_bound)
            push!(lowers, lower_bound)

            # Debug print (optional)
            println("Time: $current_time, Upper: $upper_bound, Lower: $lower_bound")
        end
    end
end



function solve_minlp_gurobi(model::Model; time_limit::Float64 = ∞)
    # Setup Gurobi optimizer
    optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => time_limit)
    set_optimizer(model, optimizer)

    # Set up callback
    MOI.set(model, Gurobi.CallbackFunction(), store_both_bounds_gurobi)

    # Solve the model
    optimize!(model)

    # Extract results
    term_status = termination_status(model)
    obj_value = if has_values(model) objective_value(model) else ∞ end

    # Record final bounds if not already captured
    best_bound = MOI.get(model, MOI.ObjectiveBound())
    best_obj = obj_value
    if isempty(uppers) || isempty(lowers) || isempty(times)
        push!(times, time() - start_time)
        push!(uppers, best_bound)
        push!(lowers, best_solution)
    end

    gap = if best_bound != ∞ && obj_value != ∞
            abs(obj_value - best_bound) / abs(obj_value)
    else
         ∞
    end
    # Print summary
    @printf("Termination Status: %s\n", term_status)
    @printf("Objective Value: %.4f\n", best_obj)
    @printf("Best Bound: %.4f\n", best_bound)
    @printf("Relative Gap: %.4f%%\n", gap * 100)

    return best_obj, best_bound, gap, term_status, times, uppers, lowers
end

function solve(P; time_limit = ∞, save_bounds=nothing)
    set_silent(P)
    set_time_limit_sec(P, time_limit)
    if save_bounds == "upper"
        MOI.set(P, Gurobi.CallbackFunction(), store_upperbound)
    elseif save_bounds == "lower"
        MOI.set(P, Gurobi.CallbackFunction(), store_lowerbound)
    end

    optimize!(P)

    term_status = termination_status(P)

    if solution_summary(P).result_count > 0
        C = -objective_value(P)
    else
        C = ∞
    end

    vars = all_variables(P)

    return C, vars, term_status
end

function milp_solver(P::GenericModel; time_limit = Inf, save_bounds = false)
    if save_bounds
        return solve(P, time_limit = time_limit, save_bounds = "upper")
    else
        return solve(P, time_limit = time_limit)
    end
end

function nlp_solver(P::GenericModel; time_limit = Inf, save_bounds = false)
    if save_bounds
        return solve(P, time_limit = time_limit, save_bounds = "lower")
    else
        return solve(P, time_limit = time_limit)
    end
end

function get_vlp_fixing_values(vars::Vector{VariableRef}, name::String, dimensions::Vector{String}; int_tol=1e-6)::Dict{String, Int64}
    fixing_values = Dict{String, Int64}()

    for x in dimensions
        ξ = [v for v in vars if startswith(string(v), "ξ_$(x)_$(name)")]
        bps = [parse(Float64, match(r"\[([0-9.]+)\]$", string(ξ_i))[1]) for ξ_i in ξ]

        # sort according to breakpoints
        sort_ix = sortperm(bps)
        ξ = ξ[sort_ix]
        bps = bps[sort_ix]

        # filter breakpoints NOT to fix (considering solver tolerance)
        ξ_not_to_fix = ξ[abs.(value.(ξ)) .> int_tol]
        if length(ξ_not_to_fix) == 1
            bp_not_to_fix = only(bps[value.(ξ) .> int_tol])

            if bp_not_to_fix == bps[1]
                # add to not-to-fix the breakpoint right after the one in `model1` solution
                push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
            else
                # add to not-to-fix the breakpoint just before the one in `model1` solution
                push!(ξ_not_to_fix, last(ξ[bps .< bp_not_to_fix]))
            end

            @assert length(ξ_not_to_fix) == 2
        end

        for ξ_i in setdiff(ξ, ξ_not_to_fix)
            fixing_values[string(ξ_i)] = 0
        end
    end

    return fixing_values
end

function get_fixing_values(vars::Vector{VariableRef}, platform::Platform; int_tol = 1e-6)::Dict{String, Int64}
    fixing_values = Dict{String, Int64}()

    # Get well variables
    for well in all_wells(platform)
        # strategic decisions
        y_value = first(value(v) for v in vars if name(v) == "y_$(well.name)")
        fixing_values["y_$(well.name)"] = round(y_value)

        t_gl_value = first(value(v) for v in vars if name(v) == "t_gl_$(well.name)")
        fixing_values["t_gl_$(well.name)"] = round(t_gl_value)

        # VLP
        merge!(fixing_values, get_vlp_fixing_values(vars, well.name, ["iglr", "whp", "qliq_vlp"], int_tol=int_tol))
    end

    # Get manifold variables
    for manifold in platform.manifolds
        m = manifold.name

        # strategic decisions
        y_value = first(value(v) for v in vars if name(v) == "y_$m")
        fixing_values["y_$m"] = round(y_value)

        # VLP
        merge!(fixing_values, get_vlp_fixing_values(vars, m, ["qliq", "gor", "wct", "iglr"], int_tol=int_tol))
    end

    return fixing_values
end

function fix_model!(P::GenericModel, fixing_values::Dict{String, Int64})
    # unfix ξ variables
    for v in all_variables(P)
        if is_fixed(v)
            unfix(v)

            # TODO: this is not necessary anymore
            # if startswith(name(v), "ξ")
            #     set_lower_bound(v, 0.0)
            # end
        end
    end

    # fix variables
    for (var_name, value) in fixing_values
        var = variable_by_name(P, var_name)
        fix(var, value, force=true)
    end
end

function exclude!(P::GenericModel, fixing_values::Dict{String, Int64}; int_tol=1e-6)
    constraint_lhs = 0
    # println("Fixing_values; ", fixing_values)
    for (var_name, value) in fixing_values
        var = variable_by_name(P, var_name)

        if startswith(var_name, "y")
            constraint_lhs += value * (1 - var) + (1 - value) * var
        end

        if startswith(var_name, "t_gl")
            constraint_lhs += value * (1 - var) + (1 - value) * var
        end

        if startswith(var_name, "z")
            constraint_lhs += value * (1 - var) + (1 - value) * var
        end
    end

    @constraint(P, constraint_lhs >= 1)
end

function check_points_is_feasible(P, vars)
    P_, _ = copy_model(P)

    set_silent(P_)
    set_optimizer(P_, () -> Gurobi.Optimizer(GRB_ENV))

    for v in vars
        if ~(startswith(name(v), "ξ") || startswith(name(v), "λ"))
            fix(variable_by_name(P_, name(v)), value(v), force=true)
        end
    end

    @objective(P_, Max, 0)

    optimize!(P_)

    if termination_status(P_) == MOI.OPTIMAL
        return true
    else
        return false
    end
end

function load_instance(scenario_id::Int, instance_id::Int; base_path="scenarios")
    path = "$base_path/scenario_$scenario_id/instance_$instance_id.json"
    try
        data = JSON3.read(read(path, String), Dict)
        
        # Load units
        kgf, g, m3, d = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day)
        
        # Rebuild satellite wells
        satellite_wells = [
            Well(
                w["name"],
                w["gor"],
                w["wct"],
                w["q_min"],
                w["q_max"],
                VLP("data/Well_$(w["vlp_type"])_VLP.Ecl"),
                IPR(w["p_res"] * kgf + g, w["IP"] * (m3/d)/kgf)
            ) for w in data["satellite_wells"]]
        
        # println(satellite_wells[1].name)
        
        # Rebuild manifolds
        manifolds = [
            Manifold(
                VLP("data/MSP_UEP_VFP.Ecl"),
                [
                    Well(
                        w["name"],
                        w["gor"],
                        w["wct"],
                        w["q_min"],
                        w["q_max"],
                        VLP("data/Well_SubseaManifold_VLP.Ecl"),
                        IPR(w["p_res"] * kgf + g, w["IP"] * (m3/d)/kgf)
                    ) for w in m["wells"]
                ],
                choke_enabled = false
            ) for m in data["manifolds"]]
        
        platform = Platform(
            data["metadata"]["p_sep"] * kgf + g,
            satellite_wells,
            manifolds,
            data["metadata"]["q_inj_max_plat"],   
            nothing,                              
            nothing,                              
            data["metadata"]["q_liq_max_plat"]   
        )
        # println("Created platform")

        return platform
    
    catch e
        @error "Error loading instance $instance_id for scenario $scenario_id" exception=(e, catch_backtrace())
        rethrow()
    end
end

function replace_nan!(data)
    if data isa AbstractDict
        for (k, v) in data
            if v isa AbstractDict || v isa AbstractArray
                replace_nan!(v)
            elseif v isa Number
                if isnan(v)
                    data[k] = nothing
                end
            end
        end
    elseif data isa AbstractArray
        for i in eachindex(data)
            item = data[i]
            if item isa AbstractDict || item isa AbstractArray
                replace_nan!(item)
            elseif item isa Number
                if isnan(item)
                    data[i] = nothing
                end
            end
        end
    end
end

function safe_variable_value(model, var_name)
    try
        var = variable_by_name(model, var_name)
        isnothing(var) ? NaN : value(var)
    catch e
        @warn "Error retrieving variable $var_name: $e"
        NaN
    end
end

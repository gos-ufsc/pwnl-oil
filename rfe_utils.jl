using Gurobi, JuMP, Printf, Oil


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

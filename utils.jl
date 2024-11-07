using Gurobi, JuMP, Printf

const ∞ = Inf

times = Vector{Float64}()
uppers = Vector{Float64}()
lowers = Vector{Float64}()

push!(times, 0.0)
push!(uppers, ∞)
push!(lowers, -∞)
cb_calls = Cint[]  # just for debugging
# Setup callback for storing the upper and lower bounds over time
function store_upperbound(cb_data, cb_where::Cint)
    # You can reference variables outside the function as normal
    push!(cb_calls, cb_where)
    # You can select where the callback is run
    if cb_where == GRB_CB_MIP
        objbnd = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBND, objbnd)
        push!(times, time() - start_time)
        push!(lowers, lowers[end])
        push!(uppers, objbnd[])
    end
    return
end
function store_lowerbound(cb_data, cb_where::Cint)
    # You can reference variables outside the function as normal
    push!(cb_calls, cb_where)
    # You can select where the callback is run
    if cb_where == GRB_CB_MIP
        objbst = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, objbst)
        push!(times, time() - start_time)
        push!(lowers, objbst[])
        push!(uppers, uppers[end])
    end
    return
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

    if solution_summary(P).result_count > 0
        C = -objective_value(P)
    else
        C = ∞
    end

    vars = all_variables(P)

    return C, vars
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

function fix_model!(P::GenericModel, vars)
    fix_tol = 1e-6
    # get all wells in platform
    all_wells = platform.satellite_wells
    for manifold in platform.manifolds
        all_wells = all_wells ∪ manifold.wells
    end

    ## FIXING SATELLITE WELL VARIABLES
    for well in all_wells
        # strategic decisions
        y_value = first(value(v) for v in vars if name(v) == "y_$(well.name)")
        fix(variable_by_name(P, "y_$(well.name)"), round(y_value), force = true)

        t_gl_value = first(value(v) for v in vars if name(v) == "t_gl_$(well.name)")
        fix(variable_by_name(P, "t_gl_$(well.name)"), round(t_gl_value), force = true)

        # VLP
        for x in ["iglr", "whp", "qliq_vlp"]
            ξ = [v for v in vars if startswith(string(v), "ξ_$(x)_$(well.name)")]
            bps = [parse(Float64, match(r"\[([0-9.]+)\]$", string(ξ_i))[1]) for ξ_i in ξ]

            # sort according to breakpoints
            sort_ix = sortperm(bps)
            ξ = ξ[sort_ix]
            bps = bps[sort_ix]
            # filter breakpoints NOT to fix (considering solver tolerance)
            ξ_not_to_fix = ξ[abs.(value.(ξ)) .> fix_tol]
            println("Before TREATMENT:ξ_not_to_fix para $(x)_$(well.name): ", ξ_not_to_fix)
            if length(ξ_not_to_fix) == 1
                bp_not_to_fix = only(bps[value.(ξ) .> fix_tol])

                if bp_not_to_fix > bps[1]
                    # add to not-to-fix the breakpoint just before the one in `model1` solution
                    push!(ξ_not_to_fix, last(ξ[bps .< bp_not_to_fix]))
                else
                    push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                end

                # if bp_not_to_fix < bps[end]
                #     # add to not-to-fix the breakpoint right after the one in `model1` solution
                #     push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                # end

                @assert length(ξ_not_to_fix) == 2
            end
            println("AFTER TREATMENT:ξ_not_to_fix para $(x)_$(well.name): ", ξ_not_to_fix)
            # fix breakpoints to 0 => restrict to the segment with the nonzero bps
            for ξ_i in ξ
                # reset previous fixing
                if is_fixed(ξ_i)
                    unfix(ξ_i)
                    set_lower_bound(ξ_i, 0.0)
                end

                if ~(ξ_i in ξ_not_to_fix)
                    fix(variable_by_name(P, string(ξ_i)), 0, force=true)
                end
            end
        end
    end

    ## FIXING RISER VARIABLES
    for manifold in platform.manifolds
        m = manifold.name

        # strategic decisions
        y_value = first(value(v) for v in vars if name(v) == "y_$m")
        fix(variable_by_name(P, "y_$m"), round(y_value))

        # VLP
        for x in ["qliq", "gor", "wct", "iglr"]
            ξ = [v for v in vars if startswith(string(v), "ξ_$(x)_$m")]
            bps = [parse(Float64, match(r"\[([0-9.]+)\]$", string(ξ_i))[1]) for ξ_i in ξ]

            # sort according to breakpoints
            sort_ix = sortperm(bps)
            ξ = ξ[sort_ix]
            bps = bps[sort_ix]
            
            # filter breakpoints NOT to fix
            ξ_not_to_fix = ξ[value.(ξ) .> fix_tol]
            println("BEFORE TREATMENT FIXING: ξ_not_to_fix para $(x)_$(m): ", ξ_not_to_fix)
            if length(ξ_not_to_fix) == 1
                bp_not_to_fix = only(bps[value.(ξ) .> fix_tol])

                if bp_not_to_fix > bps[1]
                    # add to not-to-fix the breakpoint just before the one in `model1` solution
                    push!(ξ_not_to_fix, last(ξ[bps .< bp_not_to_fix]))
                else
                    push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                end

                # if bp_not_to_fix < bps[end]
                #     # add to not-to-fix the breakpoint right after the one in `model1` solution
                #     push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                # end

                @assert length(ξ_not_to_fix) == 2
                # @assert length(ξ_not_to_fix) <= 3
            end

            println("AFTER TREATMENT FIXING: ξ_not_to_fix para $(x)_$(m): ", ξ_not_to_fix)
            # fix breakpoints to 0 => restrict to the segment with the nonzero bps
            for ξ_i in ξ
                if ~(ξ_i in ξ_not_to_fix)
                    fix(variable_by_name(P, string(ξ_i)), 0, force=true)
                end
            end
        end
    end
end

function exclude!(P::GenericModel, vars)
    constraint_lhs = 0
    fix_tol = 1e-6

    # get all wells in platform
    all_wells = platform.satellite_wells
    for manifold in platform.manifolds
        all_wells = all_wells ∪ manifold.wells
    end

    ## ACTUAL FIXING
    for well in all_wells
        # strategic decisions
        y_value = first(value(v) for v in vars if name(v) == "y_$(well.name)")
        y_var = variable_by_name(P, "y_$(well.name)")

        constraint_lhs += y_value * (1 - y_var) + (1-y_value) * y_var

        t_gl_value = first(value(v) for v in vars if name(v) == "t_gl_$(well.name)")
        t_gl_var = variable_by_name(P, "t_gl_$(well.name)")

        constraint_lhs += t_gl_value * (1 - t_gl_var) + (1 - t_gl_value) * t_gl_var

        # VLP
        for x in ["iglr", "whp", "qliq_vlp"]
        # for x in ["iglr",]
            ξ = [v for v in vars if startswith(string(v), "ξ_$(x)_$(well.name)")]
            bps = [parse(Float64, match(r"\[([0-9.]+)\]$", string(ξ_i))[1]) for ξ_i in ξ]

            # sort according to breakpoints
            sort_ix = sortperm(bps)
            ξ = ξ[sort_ix]
            bps = bps[sort_ix]

            # filter breakpoints NOT to fix
            ξ_not_to_fix = ξ[value.(ξ) .> fix_tol]
            if length(ξ_not_to_fix) == 1
                bp_not_to_fix = only(bps[value.(ξ) .> fix_tol])

                if bp_not_to_fix > bps[1]
                    # add to not-to-fix the breakpoint just before the one in `model1` solution
                    push!(ξ_not_to_fix, last(ξ[bps .< bp_not_to_fix]))
                else
                    push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                end

                # if bp_not_to_fix < bps[end]
                #     # add to not-to-fix the breakpoint right after the one in `model1` solution
                #     push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                # end

                # @assert length(ξ_not_to_fix) >= 2
                @assert length(ξ_not_to_fix) == 2
            end

            for ξ_i in ξ
                if ~(ξ_i in ξ_not_to_fix)
                    constraint_lhs += variable_by_name(P, string(ξ_i))
                end
            end
        end
    end

    ## FIXING RISER VARIABLES
    for manifold in platform.manifolds
        m = manifold.name

        # strategic decisions
        y_value = first(value(v) for v in vars if name(v) == "y_$m")
        y_var = variable_by_name(P, "y_$m")

        constraint_lhs += y_value * (1 - y_var) + (1-y_value) * y_var

        # VLP
        for x in ["qliq", "gor", "wct", "iglr"]
            ξ = [v for v in vars if startswith(string(v), "ξ_$(x)_$m")]
            bps = [parse(Float64, match(r"\[([0-9.]+)\]$", string(ξ_i))[1]) for ξ_i in ξ]

            # sort according to breakpoints
            sort_ix = sortperm(bps)
            ξ = ξ[sort_ix]
            bps = bps[sort_ix]

            # filter breakpoints NOT to fix
            ξ_not_to_fix = ξ[value.(ξ) .> fix_tol]
            if length(ξ_not_to_fix) == 1
                bp_not_to_fix = only(bps[value.(ξ) .> fix_tol])

                if bp_not_to_fix > bps[1]
                    # add to not-to-fix the breakpoint just before the one in `model1` solution
                    push!(ξ_not_to_fix, last(ξ[bps .< bp_not_to_fix]))
                else
                    push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                end

                # if bp_not_to_fix < bps[end]
                #     # add to not-to-fix the breakpoint right after the one in `model1` solution
                #     push!(ξ_not_to_fix, first(ξ[bps .> bp_not_to_fix]))
                # end
# s
                # @assert length(ξ_not_to_fix) >= 2
                # @assert length(ξ_not_to_fix) <= 3
                @assert length(ξ_not_to_fix) == 2
            end

            # fix breakpoints to 0 => restrict to the segment with the nonzero bps
            for ξ_i in ξ
                if ~(ξ_i in ξ_not_to_fix)
                    constraint_lhs += variable_by_name(P, string(ξ_i))
                end
            end
        end
    end

    @constraint(P, constraint_lhs >= fix_tol)
end

function check_points_is_feasible(P, vars)
    P_, _ = copy_model(P)

    set_silent(P_)
    set_optimizer(P_, () -> Gurobi.Optimizer(GRB_ENV))

    for v in vars
        if ~(startswith(name(v), "ξ") || startswith(name(v), "λ"))
            println(name(v), " = ", value(v))
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



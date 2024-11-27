using Oil, CSV

include("models.jl")
include("utils.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)
time_budget = 60.0 * 3  # 3 minutes budget


## Profiling
time_milp_solver = 0
time_nlp_solver = 0
time_model_manipulation = 0

# SCENARIO 1
platform = Platform(
    10.001 * kgf + g,
    # satellite_wells = Vector{Well}(),
    satellite_wells = [
        Well("P05", 45.0, 0.25, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(175.0 * kgf + g, 54.123932 * (m3 / d) / kgf)),
        Well("P13", 70.0, 0.26, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(220.0 * kgf + g, 107.64259 * (m3 / d) / kgf)),
    ],
    manifolds = [
        Manifold(
            VLP("data/MSP_UEP_VFP.Ecl"),
            [
                Well("P01", 70.0, 0.20, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(180*kgf + g, 98.11 * (m3 / d) / kgf)),
                Well("P02", 100.0, 0.40, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(200.0 * kgf + g, 48.33 * (m3 / d) / kgf)),
            ],
            choke_enabled = false
        ),
    ],
)

# ALGORITHM

## 1. Build MINLP problem P
P_minlp = get_minlp_problem(platform, sos2_with_binary = true)
C_minlp = ∞  # following the minimization standard

### Start counting time after original problem is built
start_time = time()

## 2. Build the MILP relaxation \tilde{P}
P_relax = get_milp_relaxation(platform, sos2_with_binary = true)
C_relax = ∞

## 3. Solve the MILP relaxation \tilde{P}, get solution and cost
t = time()
C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time))
time_milp_solver += time() - t

@printf("| %-5s | %-15s | %-15s | %-15s | %-10s |\n", "It.", "Lower Bound", "Upper Bound", "Gap", "Time")

global i = 0
while (C_relax < C_minlp) & (time() - start_time < time_budget)
    t = time()
    # fixing values dictate the region which we are going to explore in this iteration
    fixing_values = get_fixing_values(vars_relax, platform)

    @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)
    # 4. Build P_fixed by constraining P to x_relax
    fix_model!(P_minlp, fixing_values)
    C_fixed = ∞

    # add lower-bound valid inequality
    valid_ineq = constraint_by_name(P_minlp, "valid_ineq")
    if isnothing(valid_ineq)
        @constraint(P_minlp, -variable_by_name(P_minlp, "q_oil_total") >= C_relax)
    else
        set_normalized_rhs(valid_ineq, C_relax)
    end
    global time_model_manipulation += time() - t

    # 5. Solve P_fixed, update x and C
    t = time()
    C_fixed, vars_fixed = nlp_solver(P_minlp, time_limit = time_budget - (time() - start_time))
    global time_nlp_solver += time() - t

    if C_fixed < C_minlp
        global C_minlp = C_fixed;
        global solution = Dict(name(v) => value(v) for v in vars_fixed)
    end

    if time() - start_time >= time_budget
        break
    end

    # Get best solution from previous P_relax that respects the fixing values
    warm_start_names = []
    warm_start_values = []
    for j = 1:result_count(P_relax)
        diff = 0
        for (v_name, v_value) in fixing_values
            diff += abs(v_value - value.(variable_by_name(P_relax, v_name), result = j))
        end
        if diff > 0
            warm_start = [(name(v), value(v, result = j)) for v in all_variables(P_relax)]
            break
        end
    end

    # 6. Exclude x_relax from P_relax
    exclude!(P_relax, fixing_values)

    # Just for debugging purposes!
    # TODO: this is breaking because it tries to access the values of the variables after
    # the model has been modified.
    # @assert ~check_points_is_feasible(P_relax, vars_relax)

    # add lower-bound valid inequality
    valid_ineq = constraint_by_name(P_relax, "valid_ineq")
    if isnothing(valid_ineq)
        @constraint(P_relax, -variable_by_name(P_relax, "q_oil_total") >= C_relax)
    else
        set_normalized_rhs(valid_ineq, C_relax)
    end
    global time_model_manipulation += time() - t

    # add warm-start feasible solution
    for (v_name, v_value) in warm_start
        set_start_value(variable_by_name(P_relax, v_name), v_value)
    end

    # 7. Solve P_relax, get solution x_relax and C_relax
    t = time()
    global C_relax_new, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time))
    global time_milp_solver += time() - t

    if is_solved_and_feasible(P_relax)  # this is not the case if the time limit is reached
        C_relax = C_relax_new
    end

    global i += 1
end
final_time = time() - start_time

C_relax = min(C_relax, C_minlp)  # in case the last relax solution was worse than the candidate
@printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)

# Store final results
push!(times, final_time)
push!(lowers, -C_minlp)
push!(uppers, -C_relax)

CSV.write("FnX_latest.csv", Dict(
    "times"=>times,
    "lowers"=>lowers,
    "uppers"=>uppers
))

println("Objective: ", C_minlp)
println("Lower bound: ", C_relax)
if abs(C_minlp) > 0
    gap = max(0, abs(C_minlp - C_relax) / abs(C_minlp))
else
    gap = abs(C_minlp - C_relax)
end
println("Gap: ", gap)
println("Runtime: ", final_time)

println("\nTime spent in MILP solver: ", time_milp_solver)
println("Time spent in NLP solver: ", time_nlp_solver)
println("Time spent in model manipulation: ", time_model_manipulation)

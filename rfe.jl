using Oil, CSV

include("models.jl")
include("utils.jl")
include("rfe_utils.jl")

const GRB_ENV = Gurobi.Env()

# (re)set global variables of the callbacks
times = Vector{Float64}()
uppers = Vector{Float64}()
lowers = Vector{Float64}()

push!(times, 0.0)
push!(uppers, ∞)
push!(lowers, -∞)
cb_calls = Cint[]  # just for debugging

start_time = time()


kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)
time_budget = 60.0 * 15  # 15 minutes budget


## Profiling
time_milp_solver = 0
time_nlp_solver = 0
time_model_manipulation = 0

# PROBLEM INSTANCE
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
set_optimizer(P_minlp, () -> Gurobi.Optimizer(GRB_ENV))
C_minlp = ∞  # following the minimization standard

### Start counting time after original problem is built
start_time = time()

## 2. Build the MILP relaxation \tilde{P}
P_relax = get_milp_relaxation(platform, sos2_with_binary = true)
set_optimizer(P_relax, () -> Gurobi.Optimizer(GRB_ENV))
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

    # 6. Exclude x_relax from P_relax
    t = time()
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

    # 7. Solve P_relax, get solution x_relax and C_relax
    t = time()
    global C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time))
    global time_milp_solver += time() - t

    global i += 1
end
final_time = time() - start_time

C_relax = min(C_relax, C_minlp)  # in case the last relax solution was worse than the candidate
@printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)

# Store final results
push!(times, final_time)
push!(lowers, -C_minlp)
push!(uppers, -C_relax)

CSV.write("RFE_latest.csv", Dict(
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

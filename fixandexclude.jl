using Oil, CSV

include("models.jl")
include("utils.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)
time_budget = 60.0 * 5  # 5 minutes budget

platform = Platform(
    10.001 * kgf + g,
    satellite_wells = [
        Well("P05", 65.0, 0.55, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(175.0 * kgf + g, 44.123932 * (m3 / d) / kgf)),
        Well("P13", 70.0, 0.25, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(220.0 * kgf + g, 117.64259 * (m3 / d) / kgf)),
    ],
    riser = Riser(
        VLP("data/MSP_UEP_VFP.Ecl"),
        [
            Well("P01", 70.0, 0.10, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(200*kgf + g, 98.11 * (m3 / d) / kgf)),
            # Well("P02", 120.0, 0.50, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(180.0 * kgf + g, 48.33 * (m3 / d) / kgf)),
            # Well("P03", 100.0, 0.05, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(170.0 * kgf + g, 28.08 * (m3 / d) / kgf)),
            # Well("P04", 150.0, 0.75, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(220.0 * kgf + g, 73.01 * (m3 / d) / kgf)),
        ],
        choke_enabled = false
    ),
    q_inj_max=10*1e3
)

# ALGORITHM

## 1. Build MINLP problem P
P_minlp = get_minlp_problem(platform)
C_minlp = ∞  # following the minimization standard

### Start counting time after original problem is built
start_time = time()

## 2. Build the MILP relaxation \tilde{P}
P_relax = get_milp_relaxation(platform)
C_relax = ∞

## 3. Solve the MILP relaxation \tilde{P}, get solution and cost
C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time))

@printf("| %-5s | %-15s | %-15s | %-15s | %-10s |\n", "It.", "Lower Bound", "Upper Bound", "Gap", "Time")

global i = 0
while (C_relax < C_minlp) & (time() - start_time < time_budget)
    @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)
    # 4. Build P_fixed by constraining P to x_relax
    fix_model!(P_minlp, vars_relax)
    C_fixed = ∞

    # 5. Solve P_fixed, update x and C
    C_fixed, vars_fixed = nlp_solver(P_minlp, time_limit = time_budget - (time() - start_time))

    if C_fixed < C_minlp
        global C_minlp = C_fixed;
        global solution = Dict(name(v) => value(v) for v in vars_fixed)
    end

    if time() - start_time >= time_budget
        break
    end

    # 6. Exclude x_relax from P_relax
    exclude!(P_relax, vars_relax)

    # Just for debugging purposes!
    # @assert ~check_points_is_feasible(P_relax, vars_relax)

    # 7. Solve P_relax, get solution x_relax and C_relax
    global C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time))

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

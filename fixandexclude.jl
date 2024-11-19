using Oil, CSV

include("models.jl")
include("utils.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)
time_budget = 60.0 * 60.0 * 24  # 5 minutes budget

platform = Platform(
    10.001 * kgf + g,
    # satellite_wells = Vector{Well}(),
    satellite_wells = [
        Well("P21", 45.0, 0.30, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_UEP_VLP.Ecl"), IPR(190.0 * kgf + g, 59.55 * (m3 / d) / kgf)),
        Well("P22", 85.0, 0.10, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_UEP_VLP.Ecl"), IPR(210.0 * kgf + g, 102.33 * (m3 / d) / kgf)),
        Well("P23", 60.0, 0.44, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_UEP_VLP.Ecl"), IPR(185.0 * kgf + g, 76.22 * (m3 / d) / kgf)),
    ],
    manifolds = [
        Manifold(
            VLP("data/MSP_UEP_VFP.Ecl"),
            [
                Well("P31", 90.0, 0.15, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(195 * kgf + g, 88.21 * (m3 / d) / kgf)),
                Well("P32", 130.0, 0.47, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(215.0 * kgf + g, 55.34 * (m3 / d) / kgf)),
                Well("P33", 105.0, 0.21, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(205.0 * kgf + g, 34.12 * (m3 / d) / kgf)),
            ],
            choke_enabled = false
        ),
        Manifold(
            VLP("data/MSP_UEP_VFP.Ecl"),
            [
                Well("P41", 140.0, 0.42, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(225.0 * kgf + g, 72.45 * (m3 / d) / kgf)),
                Well("P42", 115.0, 0.27, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(175.0 * kgf + g, 43.29 * (m3 / d) / kgf)),
                Well("P43", 155.0, 0.70, 100 * 1e3 * m3 / d, 200 * 1e3 * m3 / d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(230.0 * kgf + g, 67.88 * (m3 / d) / kgf)),
            ],
            choke_enabled = false
        ),
    ],
)

# ALGORITHM
sos2_with_binary = true

## 1. Build MINLP problem P
P_minlp = get_minlp_problem(platform, sos2_with_binary)
C_minlp = ∞  # following the minimization standard

### Start counting time after original problem is built
start_time = time()

## 2. Build the MILP relaxation \tilde{P}
P_relax = get_milp_relaxation(platform, sos2_with_binary)
C_relax = ∞

## 3. Solve the MILP relaxation \tilde{P}, get solution and cost
C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time))

@printf("| %-5s | %-15s | %-15s | %-15s | %-10s |\n", "It.", "Lower Bound", "Upper Bound", "Gap", "Time")

global i = 0
while (C_relax < C_minlp) & (time() - start_time < time_budget)
    @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)
    # 4. Build P_fixed by constraining P to x_relax
    if sos2_with_binary
        fix_model_with_binary!(P_minlp, vars_relax)
    else
        fix_model!(P_minlp, vars_relax)
    end
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
    if sos2_with_binary
        exclude_with_binary!(P_relax, vars_relax)
    else
        exclude!(P_relax, vars_relax)
    end

    # Just for debugging purposes!
    # TODO: this is breaking because it tries to access the values of the variables after
    # the model has been modified.
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

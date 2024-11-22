using Oil, CSV
using DataFrames
using FileIO
include("models.jl")
include("utils.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)

function solve_minlp_with_constraints(platform::Platform, time_budget::Float64, filename::String)
    ## Profiling
    time_milp_solver = 0
    time_nlp_solver = 0
    time_model_manipulation = 0

    ## 1. Build MINLP problem P
    P_minlp = get_minlp_problem(platform, sos2_with_binary=true)
    C_minlp = ∞  # Following the minimization standard

    ## 2. Build the MILP relaxation \tilde{P}
    P_relax = get_milp_relaxation(platform, sos2_with_binary=true)
    C_relax = ∞

    ## Start the optimization process
    start_time = time()
    t = time()
    C_relax, vars_relax = milp_solver(P_relax, time_limit=time_budget - (time() - start_time))
    time_milp_solver += time() - t

    @printf("| %-5s | %-15s | %-15s | %-15s | %-10s |\n", "It.", "Lower Bound", "Upper Bound", "Gap", "Time")

    q_liq_value = 0
    q_wat_value = 0
    q_oil_value = 0
    q_inj_value = 0

    global i = 0
    while (C_relax < C_minlp) & (time() - start_time < time_budget)
        t = time()
        fixing_values = get_fixing_values(vars_relax, platform)
        @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100 * abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)

        # Fix and solve P_minlp
        fix_model!(P_minlp, fixing_values)
        valid_ineq = constraint_by_name(P_minlp, "valid_ineq")
        if isnothing(valid_ineq)
            @constraint(P_minlp, -variable_by_name(P_minlp, "q_oil_total") >= C_relax)
        else
            set_normalized_rhs(valid_ineq, C_relax)
        end
        time_model_manipulation += time() - t

        t = time()
        C_fixed, vars_fixed = nlp_solver(P_minlp, time_limit=time_budget - (time() - start_time))
        time_nlp_solver += time() - t

        if C_fixed < C_minlp
            C_minlp = C_fixed
            q_liq_plat_index = findfirst(v -> name(v) == "q_liq_plat", vars_fixed)
            q_wat_plat_index = findfirst(v -> name(v) == "q_wat_plat", vars_fixed)
            q_oil_plat_index = findfirst(v -> name(v) == "q_oil_plat", vars_fixed)
            q_inj_plat_index = findfirst(v -> name(v) == "q_inj_plat", vars_fixed)
            q_liq_value = value(vars_fixed[q_liq_plat_index])
            q_wat_value = value(vars_fixed[q_wat_plat_index])
            q_oil_value = value(vars_fixed[q_oil_plat_index])
            q_inj_value = value(vars_fixed[q_inj_plat_index])
            println("The value of q_liq_plat is: ", q_liq_value)
            println("The value of q_wat_plat is: ", q_wat_value)
            println("The value of q_oil_plat is: ", q_oil_value)
            println("The value of q_inj_plat is: ", q_inj_value)
        end

        # Exclude and solve P_relax
        exclude!(P_relax, fixing_values)
        t = time()
        C_relax, vars_relax = milp_solver(P_relax, time_limit=time_budget - (time() - start_time))
        time_milp_solver += time() - t

        i += 1
    end

    final_time = time() - start_time
    C_relax = min(C_relax, C_minlp)

    # Calculate gap
    gap = if abs(C_minlp) > 0
        max(0, abs(C_minlp - C_relax) / abs(C_minlp))
    else
        abs(C_minlp - C_relax)
    end

    # Access platform-level variable values
    # q_liq_value = value(variable_by_name(P_minlp, "q_liq_plat"))
    # q_wat_value = value(variable_by_name(P_minlp, "q_wat_plat"))
    # q_oil_value = value(variable_by_name(P_minlp, "q_oil_plat"))
    # q_inj_value = value(variable_by_name(P_minlp, "q_inj_plat"))

    # Save results to a CSV file
    results = DataFrame(
        name = ["objective_value", "lower_bound", "gap", "runtime", "q_liq", "q_wat", "q_oil", "q_inj", "time_nlp_solver","time_milp_solver", "time_manipulating_model"],
        value = [C_minlp, C_relax, gap, final_time, q_liq_value, q_wat_value, q_oil_value, q_inj_value, time_nlp_solver, time_milp_solver, time_model_manipulation]
    )

    CSV.write(filename, results)

    println("\nResults saved to: $filename")
    println("Objective: ", C_minlp)
    println("Lower bound: ", C_relax)
    println("Gap: ", gap)
    println("Runtime: ", final_time)
end

# Function to test multiple constraints
function run_tests_with_constraints(
    platform::Platform, 
    q_inj_max_values::Union{Nothing, Vector{Float64}}, 
    q_liq_max_values::Union{Nothing, Vector{Float64}}, 
    time_budget::Float64
)
    if q_liq_max_values !== nothing
        for q_liq_max in q_liq_max_values
            println("Testing with q_liq_max = ", q_liq_max)
            new_platform = Platform(platform.p_sep, platform.satellite_wells, platform.manifolds, platform.q_inj_max, platform.q_water_max, platform.q_gas_max, q_liq_max)
            solve_minlp_with_constraints(new_platform, time_budget, "results_fixandexclude_q_liq_max_$(q_liq_max).csv")
        end
    end

    if q_inj_max_values !== nothing
        for q_inj_max in q_inj_max_values
            println("Testing with q_inj_max = ", q_inj_max)
            new_platform = Platform(platform.p_sep, platform.satellite_wells, platform.manifolds, q_inj_max, platform.q_water_max, platform.q_gas_max, platform.q_liq_max)
            solve_minlp_with_constraints(new_platform, time_budget, "results_fixandexclude_q_inj_max_$(q_inj_max).csv")
        end
    end
end

# Platform Initialization
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

# Parameters for tests
q_liq_max_values = [7500.0]
q_inj_max_values = [1600.0] .* 1e3
time_budget = 60.0 * 60  # 15 minutes

# Run tests
run_tests_with_constraints(platform, q_inj_max_values, q_liq_max_values, time_budget)

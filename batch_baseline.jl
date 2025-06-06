using Oil, CSV
using DataFrames
using FileIO
include("models.jl")
include("utils.jl")

# Function to solve the optimization problem with Gurobi
function solve_gurobi(P; time_limit = ∞, save_bounds=nothing)
    # set_time_limit_sec(P, time_limit)
    C, var, status = solve(P, time_limit, save_bounds)

    if solution_summary(P).result_count > 0
        C = objective_value(P)
    else
        C = ∞
    end
    best_bound = objective_bound(P)
    gap = MOI.get(P, MOI.RelativeGap())
    vars = all_variables(P)
    solving_time = solve_time(P)
    return C, vars, best_bound, gap, solving_time
end

# Function to update platform constraints
function update_platform_constraints(platform::Platform; 
    q_inj_max::Union{Nothing, Float64} = nothing, 
    q_liq_max::Union{Nothing, Float64} = nothing)
    return Platform(
        platform.p_sep,  
        platform.satellite_wells,
        platform.manifolds,
        q_inj_max !== nothing ? q_inj_max : platform.q_inj_max,
        platform.q_water_max,
        platform.q_gas_max,
        q_liq_max !== nothing ? q_liq_max : platform.q_liq_max
    )
end

# Function to solve and save results for a specific constraint scenario
function solve_and_save_results(platform, constraint_name::String, constraint_value::Union{Nothing, Float64}, time_budget::Float64)
    global times = Vector{Float64}()
    global uppers = Vector{Float64}()
    global lowers = Vector{Float64}()

    push!(times, 0.0)
    push!(uppers, ∞)
    push!(lowers, -∞)
    cb_calls = Cint[] 
    global start_time = time()

    P_minlp = get_minlp_problem(platform,sos2_with_binary=false)
    set_optimizer_attribute(P_minlp, "FeasibilityTol", 1e-9)
    set_optimizer_attribute(P_minlp, "OptimalityTol", 1e-9)
    set_optimizer_attribute(P_minlp, "IntFeasTol", 1e-8)
   


    set_optimizer(P_minlp, () -> Gurobi.Optimizer(GRB_ENV))
    best_obj, best_bound, gap, term_status, times, uppers, lowers = solve_minlp_gurobi(P_minlp, time_limit = time_budget)
    solving_time = time() - start_time

    q_liq_value = value(variable_by_name(P_minlp, "q_liq_plat"))
    q_wat_value = value(variable_by_name(P_minlp, "q_wat_plat"))
    q_oil_value = value(variable_by_name(P_minlp, "q_oil_plat"))
    q_inj_value = value(variable_by_name(P_minlp, "q_inj_plat"))

    println("Before Cleaning:")
    println("First 10 Uppers: ", first(uppers, 10))
    println("First 10 Lowers: ", first(lowers, 10))

    # Clean bounds
    clean_bounds!(uppers, lowers)

    # Validate cleaning
    println("After Cleaning:")
    println("First 10 Uppers: ", first(uppers, 10))
    println("First 10 Lowers: ", first(lowers, 10))

    println("Lengths -> Times: $(length(times)), Uppers: $(length(uppers)), Lowers: $(length(lowers))")

    # Create file name
    filename = "results_gurobi_$(constraint_name)_$(constraint_value).csv"

    # Prepare results as a DataFrame
    results = DataFrame(
        name = ["objective_value", "lower_bound", "gap", "runtime", "q_liq", "q_wat", "q_oil", "q_inj", "constraint_name", "constraint_value", "times", "lowers", "uppers"],
        value = [best_obj, best_bound, gap, solving_time, q_liq_value, q_wat_value, q_oil_value, q_inj_value, constraint_name, constraint_value, times, lowers, uppers]
    )

    # Write results to a CSV file
    CSV.write(filename, results, bufsize=500_000_000)

    println("Results saved to: $filename")
end

# Function to run tests with different constraints
function run_tests_with_constraints(
    platform::Platform,
    q_inj_max_values::Union{Nothing, Vector{Float64}}, 
    q_liq_max_values::Union{Nothing, Vector{Float64}},
    time_budget::Float64
)
    # Test q_liq_max values
    if q_liq_max_values !== nothing
        for q_liq_max in q_liq_max_values
            println("Running test with q_liq_max = ", q_liq_max)
            new_platform = update_platform_constraints(platform, q_liq_max=q_liq_max)
            solve_and_save_results(new_platform, "q_liq_max", q_liq_max, time_budget)
        end
    end

    # Test q_inj_max values
    if q_inj_max_values !== nothing
        for q_inj_max in q_inj_max_values
            println("Running test with q_inj_max = ", q_inj_max)
            new_platform = update_platform_constraints(platform, q_inj_max=q_inj_max)
            solve_and_save_results(new_platform, "q_inj_max", q_inj_max, time_budget)
        end
    end
end

# Define units and initialize platform
kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)

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

# Define test parameters
const GRB_ENV = Gurobi.Env()
const ∞  = Inf

q_liq_max_values = nothing  # Example values for q_liq_max
q_inj_max_values = [900.0] * 1e3  # Example values for q_inj_max
time_budget = 60.0 * 60  # 1 hour

# Run tests
run_tests_with_constraints(platform, q_inj_max_values, q_liq_max_values, time_budget)

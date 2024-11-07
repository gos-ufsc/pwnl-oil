using Oil, CSV
using DataFrames
using FileIO
include("models.jl")
include("utils.jl")


function solve_gurobi(P; time_limit = ∞, save_bounds=nothing)
    set_time_limit_sec(P, time_limit)
    optimize!(P)

    if solution_summary(P).result_count > 0
        C = objective_value(P)
    else
        C = ∞
    end
    best_bound = objective_bound(P)
    println("Best Bound: ", best_bound)
    gap = MOI.get(P, MOI.RelativeGap())
    println("Relative Gap: ", gap)
    vars = all_variables(P)
    return C, vars, best_bound, gap
end

function update_platform_constraints(platform::Platform; 
    q_inj_max::Union{Nothing, Float64} = nothing, 
    q_liq_max::Union{Nothing, Float64} = nothing)
    println("Atualizando valores:")
    println("q_inj_max: ", q_inj_max)
    println("q_liq_max: ", q_liq_max)

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

function solve_and_save_results(platform, constraint_name::String, constraint_value::Union{Nothing, Float64}, time_budget::Float64)
    P_minlp = get_minlp_problem(platform)

    set_optimizer(P_minlp, () -> Gurobi.Optimizer(GRB_ENV))
    C_minlp = ∞
    start_time = time()
    C_relax = ∞
    solution = nothing
    C_minlp, vars, bound, gap = solve_gurobi(P_minlp, time_limit = time_budget - (time() - start_time), save_bounds=nothing)
    solution = Dict(name(v) => value(v) for v in vars)
    final_time = time() - start_time
    filename = "results_gurobi_$(constraint_name)_$(constraint_value)_$(round(time(), digits=2)).csv"

    if C_minlp < ∞
        println("Objective: ", C_minlp)
        println("Lower bound: ", bound)
        println("Gap: ", gap)
        println("Runtime: ", final_time)

        CSV.write(filename, Dict(
            "constraint_name" => constraint_name,
            "constraint_value"=> constraint_value,
            "objective_value" => C_minlp,
            "lower_bound" => bound,
            "relative_gap" => gap,
            "runtime" => final_time
        ))

        println("Results saved to: $filename")
    else
        println("No feasible solution found.")
    end

    println("Objective: ", C_minlp)
    println("Lower bound: ", )
    println("Gap: ", gap)
    println("Runtime: ", final_time)
    println("Results saved to: $(filename).json")
end

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

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)

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

const GRB_ENV = Gurobi.Env()
const ∞  = Inf

q_liq_max_values = [4000.0]
q_inj_max_values = [650.0, 550.0] .*1e3 

time_budget = 60.0 * 60.0 * 1
run_tests_with_constraints(platform, q_inj_max_values, q_liq_max_values, time_budget)

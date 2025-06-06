using CSV, DataFrames, Gurobi
using JSON3
using Oil
using JuMP
include("models.jl")
include("utils.jl")

const GRB_ENV = Gurobi.Env()


function solve_and_store(platform, scenario_id, instance_id; 
                         results_dir="results", time_budget=30.0, sos2_with_binaries=false)
    # Create results directory
    scenario_results_dir = "$results_dir/scenario_$scenario_id"
    isdir(scenario_results_dir) || mkpath(scenario_results_dir)

    # Initialize tracking variables
    results = Dict(
        :scenario_id => scenario_id,
        :instance_id => instance_id,
        :status => "Error",
        :objective_value => NaN,
        :lower_bound => NaN,
        :gap => NaN,
        :solve_time => 0.0,
        :production => Dict()
    )

    try
        # Create and configure model
        model = get_minlp_problem(platform, sos2_with_binary=sos2_with_binaries)
        set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV))

        # Clean bounds
        clean_bounds!(uppers, lowers)

        best_obj, best_bound, gap, term_status, times, uppers, lowers = solve_minlp_gurobi(model, time_limit = time_budget)
        solving_time = time() - start_time

        println("Before Cleaning:")
        println("First 10 Uppers: ", first(uppers, 10))
        println("First 10 Lowers: ", first(lowers, 10))

        # Clean bounds
        clean_bounds!(uppers, lowers)

        # Collect production data safely
        production_data = Dict(
            :q_liq => safe_variable_value(model, "q_liq_plat"),
            :q_wat => safe_variable_value(model, "q_wat_plat"),
            :q_oil => safe_variable_value(model, "q_oil_plat"),
            :q_inj => safe_variable_value(model, "q_inj_plat")
        )

        # Update results dictionary
        results = merge(results, Dict(
            :status => string(term_status),
            :objective_value => best_obj,
            :lower_bound => best_bound,
            :gap => gap,
            :solve_time => solving_time,
            :production => production_data,
            :uppers => uppers,        # Fixed: Added => operator
            :lowers => lowers,        # Fixed: Added => operator
            :times => times           # Fixed: Added => operator
        ))
    catch e
        @error "Error solving scenario $scenario_id instance $instance_id" exception=(e, catch_backtrace())
        results = merge(results, Dict(
            :status => "Error: $(sprint(showerror, e))",
            :solve_time => time() - start_time
        ))
    end

    replace_nan!(results)

    # Save to JSON
    output_path = "$scenario_results_dir/instance_$instance_id.json"
    try
        open(output_path, "w") do f
            JSON3.pretty(f, results)
        end
    catch e
        @error "Failed to save results for scenario $scenario_id instance $instance_id" exception=(e, catch_backtrace())
    end
    
    return results
end

base_path = "scenarios"
results_dir = "results"
time_budget = 3600.0
sos2_with_binaries = false

# Create results directory structure
isdir(results_dir) || mkpath(results_dir)

# Process each scenario
for scenario_dir in readdir(base_path)
    scenario_path = joinpath(base_path, scenario_dir)
    isdir(scenario_path) || continue
    
    try
        # Extract scenario ID
        scenario_id = parse(Int, split(scenario_dir, "_")[2])

        if scenario_id > 6
            continue
        end
        
        # Process each instance
        for instance_file in readdir(scenario_path)
            instance_path = joinpath(scenario_path, instance_file)
            endswith(instance_path, ".json") || continue
            
            try
                instance_id = parse(Int, match(r"instance_(\d+)\.json", instance_file)[1])
                if instance_id > 1
                    continue  # only look at the first instance
                end
                println("Processing scenario $scenario_id instance $instance_id")
                
                platform = load_instance(scenario_id, instance_id; base_path)
                solve_and_store(platform, scenario_id, instance_id; results_dir, time_budget, sos2_with_binaries=sos2_with_binaries)
                println("Solved scenario $scenario_id instance $instance_id")
            
            catch e
                @error "Error processing instance $instance_file" exception=(e, catch_backtrace())
            end
        end
    
    catch e
        @error "Error processing scenario directory $scenario_dir" exception=(e, catch_backtrace())
    end
end
using CSV, DataFrames, SCIP
using JSON3
using Oil
using JuMP
include("models.jl")
include("utils.jl")


function solve_and_store(platform, scenario_id, instance_id; 
    results_dir="results", time_budget=3600.0)

    # Create results directory path
    scenario_results_dir = joinpath(results_dir, "scenario_$scenario_id")
    output_path = joinpath(scenario_results_dir, "instance_$instance_id.json")
    log_path = joinpath(scenario_results_dir, "instance_$instance_id.log")  # New log file path

    # Skip if file exists (regardless of content)
    if isfile(output_path) && isfile(log_path)
        println("Skipping scenario $scenario_id instance $instance_id (solution exists)")
        return nothing
    end

    # Initialize results structure
    results = Dict(
        :scenario_id => scenario_id,
        :instance_id => instance_id,
        :status => "Error",
        :objective_value => NaN,
        :solve_time => 0.0,
        :relative_gap => NaN,
        :production => Dict(
            :q_liq => NaN,
            :q_wat => NaN,
            :q_oil => NaN,
            :q_inj => NaN
        )
    )

    # Create results directory if needed
    isdir(scenario_results_dir) || mkpath(scenario_results_dir)
    
    # Open log file for writing
    log_file = open(log_path, "w")
    
    try

        # Redirect stdout and stderr to log file
        original_stdout = stdout
        original_stderr = stderr
        redirect_stdout(log_file)
        redirect_stderr(log_file)

        flush(log_file)

        # Create and configure model
        model = get_minlp_problem(platform, sos2_with_binary=true)
        set_optimizer(model, SCIP.Optimizer)
        set_optimizer_attribute(model, "limits/time", time_budget)
        version_info = MOI.get(model, MOI.SolverVersion())
        println("SCIP version: ", version_info)
        
        # Solve the model
        start_time = time()
        optimize!(model)
        solving_time = time() - start_time
        term_status = termination_status(model)

        # Determine if solved and calculate gap
        solved_statuses = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.TIME_LIMIT, MOI.SOLUTION_LIMIT]
        is_solved = term_status in solved_statuses && primal_status(model) == MOI.FEASIBLE_POINT
        dual_bound = objective_bound(model)
        obj_value = is_solved ? objective_value(model) : NaN
        relative_gap = NaN

        if is_solved
            sense = objective_sense(model)
            abs_gap = sense == MOI.MIN_SENSE ? max(0.0, obj_value - dual_bound) :
                      sense == MOI.MAX_SENSE ? max(0.0, dual_bound - obj_value) : 0.0
            denominator = max(abs(obj_value), 1e-6)
            relative_gap = abs_gap / denominator
        end

        # Collect production data
        production_data = Dict(
            :q_liq => safe_variable_value(model, "q_liq_plat"),
            :q_wat => safe_variable_value(model, "q_wat_plat"),
            :q_oil => safe_variable_value(model, "q_oil_plat"),
            :q_inj => safe_variable_value(model, "q_inj_plat")
        )

        # Update results
        results = merge(results, Dict(
            :status => "$term_status",
            :objective_value => obj_value,
            :solve_time => solving_time,
            :relative_gap => relative_gap,
            :production => production_data
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
    try
        open(output_path, "w") do f
            JSON3.pretty(f, results)
        end
        # redirect_stdout(original_stdout)
        # redirect_stderr(original_stderr)
        
        println("\n" * "="^80)
        println("FINISHED SCENARIO $scenario_id INSTANCE $instance_id")
        close(log_file)


    catch e
        @error "Failed to save results" exception=(e, catch_backtrace())
        println(log_file, "ERROR: ", sprint(showerror, e))
        println(log_file, "Stacktrace:")
        println(log_file, stacktrace(catch_backtrace()))
    end
    
    return results
end

base_path = "scenarios"
results_dir = "results"
time_budget = 3600.0

# Create results directory structure
isdir(results_dir) || mkpath(results_dir)

# Process each scenario
for scenario_dir in readdir(base_path)
    scenario_path = joinpath(base_path, scenario_dir)
    isdir(scenario_path) || continue
    try
        # Extract scenario ID
        scenario_id = parse(Int, split(scenario_dir, "_")[2])
        
        # Process each instance
        for instance_file in readdir(scenario_path)
            instance_path = joinpath(scenario_path, instance_file)
            endswith(instance_path, ".json") || continue
            
            try
                instance_id = parse(Int, match(r"instance_(\d+)\.json", instance_file)[1])
                println("Processing scenario $scenario_id instance $instance_id")
                
                platform = load_instance(scenario_id, instance_id; base_path)
                solve_and_store(platform, scenario_id, instance_id; results_dir, time_budget)
                println("Solved scenario $scenario_id instance $instance_id")
            
            catch e
                @error "Error processing instance $instance_file" exception=(e, catch_backtrace())
            end
        end
    
    catch e
        @error "Error processing scenario directory $scenario_dir" exception=(e, catch_backtrace())
    end
end

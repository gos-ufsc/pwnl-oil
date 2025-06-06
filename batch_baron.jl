using CSV, DataFrames
using JSON3
using Oil
using JuMP, AmplNLWriter
using MathOptInterface.FileFormats
include("models.jl")
include("utils.jl")

# ENV["BARON_EXEC"] = "/opt/ampl/baron"
# ENV["BARON_LICENSE"] = "/opt/ampl/ampl.lic"


function write_nl_file(platform; output_path="model.nl")
    model = get_minlp_problem(platform, sos2_with_binary=true)
    
    # Create NL format container
    nl_optimizer = MOI.FileFormats.Model(format=MOI.FileFormats.FORMAT_NL)
    
    # Copy model to NL container
    MOI.copy_to(nl_optimizer, backend(model))
    
    # Add metadata
    MOI.set(nl_optimizer, MOI.Name(), "oil_optimization_model")
    
    # Write file
    try
        MOI.write_to_file(nl_optimizer, output_path)
        @info "Successfully wrote NL file: $output_path"
        return output_path
    catch e
        @error "Failed to write NL file" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

function solve_and_store(platform, scenario_id, instance_id;
                         results_dir="results", time_budget=3600.0,
                         generate_nl_only=false)

    scenario_results_dir = joinpath(results_dir, "scenario_$scenario_id")
    output_path = joinpath(scenario_results_dir, "instance_$instance_id.json")
    nl_path = joinpath(scenario_results_dir, "instance_$instance_id.nl")

    # Generate NL file only
    if generate_nl_only
        if isfile(nl_path)
            println("Skipping NL generation: $nl_path exists")
            return nothing
        end
        mkpath(scenario_results_dir)
        write_nl_file(platform; output_path=nl_path)
        println("Generated NL file: $nl_path")
        return nothing
    end

    isdir(scenario_results_dir) || mkpath(scenario_results_dir)

    results = Dict(
        :scenario_id => scenario_id,
        :instance_id => instance_id,
        :status => "Error",
        :objective_value => NaN,
        :solve_time => 0.0
    )

    try
        # Cria o modelo com BARON via AMPL
        model = get_minlp_problem(platform, sos2_with_binary=true)
        println("MODEL CREATED")
        set_optimizer(model, () -> AmplNLWriter.Optimizer("/opt/ampl/baron"))
        set_optimizer_attribute(model, "maxtime", time_budget)

        start_time = time()
        optimize!(model)
        solving_time = time() - start_time

        term_status = termination_status(model)
        is_solved = term_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.TIME_LIMIT, MOI.SOLUTION_LIMIT] &&
                    primal_status(model) == MOI.FEASIBLE_POINT

        results[:status] = string(term_status)
        results[:solve_time] = solving_time
        results[:objective_value] = is_solved ? objective_value(model) : NaN

    catch e
        @error "Error solving scenario $scenario_id instance $instance_id" exception=(e, catch_backtrace())
        results[:status] = "Error: $(sprint(showerror, e))"
        results[:solve_time] = time() - start_time
    end

    replace_nan!(results)

    try
        open(output_path, "w") do f
            JSON3.pretty(f, results)
        end
    catch e
        @error "Failed to save results" exception=(e, catch_backtrace())
    end

    return results
end

base_path = "scenarios", 
results_dir = "results", 
time_budget = 3600.0,
generate_nl_only = true

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
                solve_and_store(platform, scenario_id, instance_id; results_dir, time_budget, generate_nl_only)
                println("Solved scenario $scenario_id instance $instance_id")
            
            catch e
                @error "Error processing instance $instance_file" exception=(e, catch_backtrace())
            end
        end
    
    catch e
        @error "Error processing scenario directory $scenario_dir" exception=(e, catch_backtrace())
    end
end
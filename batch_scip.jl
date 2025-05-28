using CSV, DataFrames, SCIP
using JSON3
using Oil
using JuMP
include("models.jl")
include("utils.jl")

function safe_variable_value(model, var_name)
    try
        var = variable_by_name(model, var_name)
        isnothing(var) ? NaN : value(var)
    catch e
        @warn "Error retrieving variable $var_name: $e"
        NaN
    end
end

function solve_and_store(platform, scenario_id, instance_id; 
                        results_dir="results", time_budget=3600.0)
    # Create results directory
    scenario_results_dir = "$results_dir/scenario_$scenario_id"
    isdir(scenario_results_dir) || mkpath(scenario_results_dir)

    # Initialize results structure
    results = Dict(
        :scenario_id => scenario_id,
        :instance_id => instance_id,
        :status => "Error",
        :objective_value => NaN,
        :solve_time => 0.0,
        :production => Dict()
    )

    try
        # Create and configure model
        model = get_minlp_problem(platform, sos2_with_binary=false)
        set_optimizer(model, SCIP.Optimizer)
        
        # Set SCIP parameters
        set_optimizer_attribute(model, "limits/time", time_budget)
        # set_optimizer_attribute(model, "display/verblevel", 0)

        # Solve the model
        global start_time = time()
        optimize!(model)
        solving_time = time() - start_time
        obj_value = NaN
        relative_gap = NaN

        # Collect results
        term_status = termination_status(model)
        is_solved = is_solved_and_feasible(model)
        dual_bound = objective_bound(model)

        if is_solved
            obj_value = objective_value(model)
            denominator = abs(obj_value) > 1e-6 ? abs(obj_value) : 1e-6
            relative_gap = (dual_bound - obj_value) / denominator
        end
        
        production_data = Dict(
            :q_liq => safe_variable_value(model, "q_liq_plat"),
            :q_wat => safe_variable_value(model, "q_wat_plat"),
            :q_oil => safe_variable_value(model, "q_oil_plat"),
            :q_inj => safe_variable_value(model, "q_inj_plat")
        )

        # Update results dictionary
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
    output_path = "$scenario_results_dir/instance_$instance_id.json"
    try
        open(output_path, "w") do f
            JSON3.pretty(f, results)
        end
    catch e
        @error "Failed to save results" exception=(e, catch_backtrace())
    end
    
    return results
end

function process_all_scenarios(; base_path="scenarios", results_dir="results", time_budget=3600.0)
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
end

function load_instance(scenario_id::Int, instance_id::Int; base_path="scenarios")
    path = "$base_path/scenario_$scenario_id/instance_$instance_id.json"
    try
        data = JSON3.read(read(path, String), Dict)
        
        # Load units
        kgf, g, m3, d = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day)
        
        # Rebuild satellite wells
        satellite_wells = [
            Well(
                w["name"],
                w["gor"],
                w["wct"],
                w["q_min"],
                w["q_max"],
                VLP("data/Well_$(w["vlp_type"])_VLP.Ecl"),
                IPR(w["p_res"] * kgf + g, w["IP"] * (m3/d)/kgf)
            ) for w in data["satellite_wells"]]
        
        # println(satellite_wells[1].name)
        
        # Rebuild manifolds
        manifolds = [
            Manifold(
                VLP("data/MSP_UEP_VFP.Ecl"),
                [
                    Well(
                        w["name"],
                        w["gor"],
                        w["wct"],
                        w["q_min"],
                        w["q_max"],
                        VLP("data/Well_SubseaManifold_VLP.Ecl"),
                        IPR(w["p_res"] * kgf + g, w["IP"] * (m3/d)/kgf)
                    ) for w in m["wells"]
                ],
                choke_enabled = false
            ) for m in data["manifolds"]]
        
        platform = Platform(
            data["metadata"]["p_sep"] * kgf + g,
            satellite_wells,
            manifolds,
            data["metadata"]["q_inj_max_plat"],   
            nothing,                              
            nothing,                              
            data["metadata"]["q_liq_max_plat"]   
        )
        # println("Created platform")

        return platform
    
    catch e
        @error "Error loading instance $instance_id for scenario $scenario_id" exception=(e, catch_backtrace())
        rethrow()
    end
end

function replace_nan!(data)
    if data isa AbstractDict
        for (k, v) in data
            if v isa AbstractDict || v isa AbstractArray
                replace_nan!(v)
            elseif v isa Number
                if isnan(v)
                    data[k] = nothing
                end
            end
        end
    elseif data isa AbstractArray
        for i in eachindex(data)
            item = data[i]
            if item isa AbstractDict || item isa AbstractArray
                replace_nan!(item)
            elseif item isa Number
                if isnan(item)
                    data[i] = nothing
                end
            end
        end
    end
end

# Run processing
process_all_scenarios()
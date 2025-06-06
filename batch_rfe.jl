using CSV, DataFrames
using JSON3
using Oil
using JuMP
include("models.jl")
include("utils.jl")
include("rfe_utils.jl")

const GRB_ENV = Gurobi.Env()

function sanitize_data!(data)
    if data isa AbstractDict
        for (k, v) in data
            if v isa Number
                if isnan(v) || isinf(v)
                    data[k] = nothing
                end
            elseif v isa AbstractDict || v isa AbstractArray
                sanitize_data!(v)
            end
        end
    elseif data isa AbstractArray
        for i in eachindex(data)
            v = data[i]
            if v isa Number
                if isnan(v) || isinf(v)
                    data[i] = nothing
                end
            elseif v isa AbstractDict || v isa AbstractArray
                sanitize_data!(v)
            end
        end
    end
end

function solve_and_store(platform, scenario_id, instance_id;
                         results_dir="results_rfe", time_budget=3600.0)
    # Setup directories
    scenario_dir = joinpath(results_dir, "scenario_$scenario_id")
    output_path = joinpath(scenario_dir, "instance_$instance_id.json")
    mkpath(scenario_dir)

    # Skip if already solved
    if isfile(output_path)
        println("Skipping scenario $scenario_id instance $instance_id (solution exists)")
        return nothing
    end

    println("\n" * "="^80)
    println("STARTING SCENARIO $scenario_id INSTANCE $instance_id")
    println("="^80)

    # Initialize results structure
    results = Dict{Symbol, Any}(
        :scenario_id => scenario_id,
        :instance_id => instance_id,
        :status => "Error",
        :objective_value => nothing,
        :lower_bound => nothing,
        :upper_bound => nothing,
        :gap => nothing,
        :solve_time => 0.0,
        :time_milp_solver => 0.0,
        :time_nlp_solver => 0.0,
        :time_model_manipulation => 0.0,
        :iterations => [],
        :production => Dict{Symbol, Any}(
        :q_liq => nothing,
        :q_wat => nothing,
        :q_oil => nothing,
        :q_inj => nothing
        )
    )

    # Timers and state variables
    start_time = time()
    total_time_milp = 0.0
    total_time_nlp = 0.0
    total_time_manip = 0.0
    iterations_data = []
    solution = nothing
    C_minlp = Inf   # Best upper bound (minimization)
    C_relax = Inf   # Best lower bound (relaxation)
    best_dual_bound = -Inf
    gap = 1.0

    # Print header for iteration table
    println("\nIteration Progress:")
    println("| Iter | Lower Bound  | Upper Bound  | Gap (%)      | MILP Time  | NLP Time   | Manip Time | Total Time |")

    try
        # Initialize models
        println("- Creating models...")
        P_minlp = get_minlp_problem(platform, sos2_with_binary=true)
        P_relax = get_milp_relaxation(platform, sos2_with_binary=true)
        set_optimizer(P_minlp, () -> Gurobi.Optimizer(GRB_ENV))
        set_optimizer(P_relax, () -> Gurobi.Optimizer(GRB_ENV))

        # Solve initial relaxation
        println("- Solving initial MILP relaxation...")
        t_milp = time()
        C_relax, vars_relax, status = milp_solver(P_relax, time_limit=time_budget)
        iter_milp_time = time() - t_milp
        total_time_milp += iter_milp_time

        if status == :OPTIMAL
            best_dual_bound = max(best_dual_bound, C_relax)
        end
    
        # Calculate and print gap
        gap = if C_minlp < Inf && best_dual_bound > -Inf
            abs(C_minlp - best_dual_bound) / max(abs(C_minlp), 1e-6) * 100  # Convert to percentage
        else
            100.0  # Default to 100% gap when no valid solution exists
        end
        elapsed = time() - start_time
        @printf("| %4d | %12.4f | %12.4f | %11.4f%% | %9.3fs | %10s | %10s | %9.3fs |\n",
        0, C_relax, C_minlp, gap, iter_milp_time, "N/A", "N/A", elapsed)

        # Record initial state
        push!(iterations_data, Dict{Symbol, Any}(
            :iteration => 0,
            :time => elapsed,
            :lower_bound => C_relax,
            :upper_bound => C_minlp,
            :milp_time => iter_milp_time,
            :nlp_time => nothing,
            :manipulation_time => nothing
        ))

        # Main RFE loop
        i = 0
        while (best_dual_bound < C_minlp) && (time() - start_time < time_budget) && (gap > 1e-3)
            i += 1
            iter_start = time()

            # Get fixing values and update model
            t_manip = time()
            fixing_values = get_fixing_values(vars_relax, platform)
            fix_model!(P_minlp, fixing_values)

            # Add valid inequality
            valid_ineq = constraint_by_name(P_minlp, "valid_ineq")
            if isnothing(valid_ineq)
                @constraint(P_minlp, -variable_by_name(P_minlp, "q_oil_total") >= C_relax)
            else
                set_normalized_rhs(valid_ineq, C_relax)
            end
            manip_time1 = time() - t_manip

            # Solve NLP
            t_nlp = time()
            C_fixed, vars_fixed = nlp_solver(P_minlp, time_limit=time_budget - (time() - start_time))
            iter_nlp_time = time() - t_nlp
            total_time_nlp += iter_nlp_time

            # Update best solution
            if C_fixed < C_minlp
                C_minlp = C_fixed
                solution = Dict(name(v) => value(v) for v in vars_fixed)
                # println("  ! New best solution found: $C_minlp")
            end

            # Update relaxation model
            t_manip = time()
            exclude!(P_relax, fixing_values)

            # Add valid inequality to relaxation
            valid_ineq = constraint_by_name(P_relax, "valid_ineq")
            if isnothing(valid_ineq)
                @constraint(P_relax, -variable_by_name(P_relax, "q_oil_total") >= C_relax)
            else
                set_normalized_rhs(valid_ineq, C_relax)
            end
            manip_time2 = time() - t_manip
            iter_manip_time = manip_time1 + manip_time2
            total_time_manip += iter_manip_time

            # Solve relaxation
            t_milp = time()
            C_relax, vars_relax = milp_solver(P_relax, time_limit=time_budget - (time() - start_time))
            iter_milp_time = time() - t_milp
            total_time_milp += iter_milp_time

            if status == MOI.OPTIMAL
                best_dual_bound = max(best_dual_bound, C_relax)
            end

            # Calculate gap and timing
            valid_lower_bound = best_dual_bound > -Inf ? best_dual_bound : -Inf
            gap = if C_minlp < Inf && best_dual_bound > -Inf
                abs(C_minlp - best_dual_bound) / max(abs(C_minlp), 1e-6) * 100  # Convert to percentage
            else
                100.0  # Default to 100% gap when no valid solution exists
            end

            elapsed = time() - start_time
            iter_total_time = time() - iter_start

            # Print iteration summary
            @printf("| %4d | %12.4f | %12.4f | %11.4f%% | %9.3fs | %9.3fs | %9.3fs | %9.3fs |\n",
            i, C_relax, C_minlp, gap, iter_milp_time, iter_nlp_time, 
            iter_manip_time, elapsed)

            # Record iteration
            push!(iterations_data, Dict{Symbol, Any}(
                :iteration => i,
                :time => elapsed,
                :lower_bound => valid_lower_bound,
                :upper_bound => C_minlp,
                :milp_time => iter_milp_time,
                :nlp_time => iter_nlp_time,
                :manipulation_time => iter_manip_time
            ))
        end

        # Final bounds adjustment
        C_relax = min(C_relax, C_minlp)
        final_time = time() - start_time
        valid_lower_bound = best_dual_bound > -Inf ? best_dual_bound : -Inf
        gap = if C_minlp < Inf && best_dual_bound > -Inf
            abs(C_minlp - best_dual_bound) / max(abs(C_minlp), 1e-6)
        else
            1.0  # 100% gap
        end

        # Update results
        status = time() - start_time >= time_budget ? "TIME_LIMIT" : "CONVERGED"
        results[:status] = status
        results[:objective_value] = -C_minlp  # Convert to positive oil production
        results[:lower_bound] = -C_minlp
        results[:upper_bound] = -C_relax
        results[:gap] = gap
        results[:solve_time] = final_time
        results[:time_milp_solver] = total_time_milp
        results[:time_nlp_solver] = total_time_nlp
        results[:time_model_manipulation] = total_time_manip
        results[:iterations] = iterations_data

        # Add production data if solution exists
        if solution !== nothing
            results[:production][:q_liq] = get(solution, "q_liq_plat", nothing)
            results[:production][:q_wat] = get(solution, "q_wat_plat", nothing)
            results[:production][:q_oil] = get(solution, "q_oil_plat", nothing)
            results[:production][:q_inj] = get(solution, "q_inj_plat", nothing)
        end

        # Print final summary
        println("\n" * "="^80)
        println("FINISHED SCENARIO $scenario_id INSTANCE $instance_id")
        println("Status: $status")
        @printf("Objective: %.4f\n", -C_minlp)
        @printf("Lower bound: %.4f\n", -C_minlp)
        @printf("Upper bound: %.4f\n", -C_relax)
        @printf("Gap: %.4f%%\n", gap * 100)
        @printf("Total time: %.2f seconds\n", final_time)
        @printf("Time MILP: %.2f seconds (%.1f%%)\n", total_time_milp, total_time_milp/final_time*100)
        @printf("Time NLP: %.2f seconds (%.1f%%)\n", total_time_nlp, total_time_nlp/final_time*100)
        @printf("Time Manip: %.2f seconds (%.1f%%)\n", total_time_manip, total_time_manip/final_time*100)
        println("="^80)

    catch e
    @error "Error solving scenario $scenario_id instance $instance_id" exception=(e, catch_backtrace())
    results[:status] = "Error: $(sprint(showerror, e))"
    results[:solve_time] = time() - start_time

    println("\n" * "!"^80)
    println("ERROR IN SCENARIO $scenario_id INSTANCE $instance_id")
    showerror(stdout, e)
    println("\n" * "!"^80)
    end

    # Clean data and save
    sanitize_data!(results)
    open(output_path, "w") do f
    JSON3.pretty(f, results)
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

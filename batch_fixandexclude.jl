using Oil, CSV
using DataFrames
using FileIO
include("models.jl")
include("utils.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)

# function validate_data_before_writing(times::Vector{Float64}, lowers::Vector{Float64}, uppers::Vector{Float64})
#     println("\nValidating Data Before Writing...")

#     # Check lengths
#     if length(times) != length(lowers) || length(times) != length(uppers)
#         println("Error: Mismatched lengths -> Times: $(length(times)), Lowers: $(length(lowers)), Uppers: $(length(uppers))")
#         return false
#     end

#     # Check for NaN or Inf in times
#     if any(isnan, times)
#         println("Warning: NaN found in 'times'.")
#     end
#     if any(x -> x == Inf || x == -Inf, times)
#         println("Warning: Inf/-Inf found in 'times'.")
#     end

#     # Check for NaN or Inf in lowers
#     if any(isnan, lowers)
#         println("Warning: NaN found in 'lowers'.")
#     end
#     if any(x -> x == Inf || x == -Inf, lowers)
#         println("Warning: Inf/-Inf found in 'lowers'.")
#     end

#     # Check for NaN or Inf in uppers
#     if any(isnan, uppers)
#         println("Warning: NaN found in 'uppers'.")
#     end
#     if any(x -> x == Inf || x == -Inf, uppers)
#         println("Warning: Inf/-Inf found in 'uppers'.")
#     end

#     println("Validation completed. Data looks good for writing.")
#     return true
# end


# function debug_written_results(filename::String)
#     # Read the CSV file into a DataFrame
#     results_df = CSV.read(filename, DataFrame)
#     println("\nLoaded Results from $filename:")
#     println(first(results_df, 10))  # Display the first 10 rows for validation

#     # Validate the structure of the DataFrame
#     println("Column names: ", names(results_df))
#     println("Column types: ", map(eltype, eachcol(results_df)))
#     println("Number of rows: ", size(results_df, 1))

#     # Helper function to parse array-like strings safely
#     function parse_array_string(array_str::String)
#         try
#             # Remove brackets and parse into a Float64 array
#             array_str = replace(array_str, r"^\[" => "", r"\]$" => "")  # Strip brackets
#             array_values = split(array_str, ",")
#             return parse.(Float64, array_values)
#         catch e
#             # println("Error parsing array: $e\nValue: $array_str")
#             return Float64[]  # Return an empty array if parsing fails
#         end
#     end

#     # Extract and parse the arrays from the DataFrame
#     if "value" in names(results_df)
#         println("Data loaded correctly. Validating arrays...")
    
#         times = parse_array_string(results_df[results_df.name .== "times", :value][1])
#         lowers = parse_array_string(results_df[results_df.name .== "lowers", :value][1])
#         uppers = parse_array_string(results_df[results_df.name .== "uppers", :value][1])
    
#         println("Lengths after reading -> Times: $(length(times)), Lowers: $(length(lowers)), Uppers: $(length(uppers))")
#     end

#         # Debugging information
#         println("Lengths -> Times: $(length(times)), Lowers: $(length(lowers)), Uppers: $(length(uppers))")
#         println("First 10 Times: ", first(times, min(length(times), 10)))
#         println("First 10 Lowers: ", first(lowers, min(length(lowers), 10)))
#         println("First 10 Uppers: ", first(uppers, min(length(uppers), 10)))
#     else
#         println("No 'value' column found in the DataFrame. Check CSV file format.")
#     end

#     return results_df
# end


function solve_minlp_with_constraints(platform::Platform, time_budget::Float64, filename::String)
    ## Profiling
    time_milp_solver = 0
    time_nlp_solver = 0
    time_model_manipulation = 0

    global times = Vector{Float64}()
    global uppers = Vector{Float64}()
    global lowers = Vector{Float64}()

    push!(times, 0.0)
    push!(uppers, 1e5)
    push!(lowers, -1e5)
    global cb_calls = Cint[]  # just for debugging
    

    ## 1. Build MINLP problem P
    P_minlp = get_minlp_problem(platform, sos2_with_binary=true)
    C_minlp = ∞  # Following the minimization standard

    ## 2. Build the MILP relaxation \tilde{P}
    P_relax = get_milp_relaxation(platform, sos2_with_binary=true)
    # set_optimizer_attribute(P_relax, "MIPGap", 1e-6)
    C_relax = ∞
    best_lower_bound = -10e5

    ## Start the optimization process
    global start_time = time()
    println("global def of start time: ", start_time)
    t = time()
    C_relax, vars_relax, milp_status = milp_solver(P_relax, time_limit=time_budget - (time() - start_time))
    push!(times, (time()-start_time))
    push!(uppers, uppers[end])
    push!(lowers, C_relax)
    time_milp_solver += time() - t

    @printf("| %-5s | %-15s | %-15s | %-15s | %-10s |\n", "It.", "Lower Bound", "Upper Bound", "Gap", "Time")

    q_liq_value = 0
    q_wat_value = 0
    q_oil_value = 0
    q_inj_value = 0

    global i = 0
    while (C_relax < C_minlp) & (time() - start_time < time_budget) & (milp_status == MOI.OPTIMAL)
        t = time()
        fixing_values = get_fixing_values(vars_relax, platform)
        # println("Time being printed: ", time()-start_time)
        @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100 * abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)

        # Fix and solve P_minlp
        fix_model!(P_minlp, fixing_values)
        # valid_ineq = constraint_by_name(P_minlp, "valid_ineq")
        # if isnothing(valid_ineq)
        #     @constraint(P_minlp, -variable_by_name(P_minlp, "q_oil_total") >= C_relax)
        # else
        #     set_normalized_rhs(valid_ineq, C_relax)
        # end
        time_model_manipulation += time() - t

        t = time()
        C_fixed, vars_fixed, nlp_status = nlp_solver(P_minlp, time_limit=time_budget - (time() - start_time), save_bounds=false)
        time_nlp_solver += time() - t

        if C_fixed < C_minlp
            push!(times, (time()-start_time))
            push!(uppers, C_fixed)
            push!(lowers, lowers[end])
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
        # println("Current time: ", t)
        # println("Current time - start_time:", (time() - start_time))
        # println("Current time limit: ", time_budget - (time() - start_time))

        if (time_budget - (time() - start_time)) < 60
            println("Breaking before solving MILP doe to time limit")
            i += 1
            break 
        end
        
        C_relax, vars_relax, milp_status = milp_solver(P_relax, time_limit=time_budget - (time() - start_time),save_bounds=false)
        println("Termination status: ", milp_status)
        time_milp_solver += time() - t

        if (C_relax > best_lower_bound) & (milp_status == MOI.OPTIMAL)
            best_lower_bound = C_relax
            push!(times, (time()-start_time))
            push!(uppers, uppers[end])
            push!(lowers, C_relax)
        end

        i += 1
    end
    println("VALUE O MILP FINAL ITERATION: ", C_relax)
    final_time = time() - start_time
    C_relax = min(best_lower_bound, C_minlp)

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

    # Validate results
    println("Before Cleaning:")
    println("First 10 Uppers: ", first(uppers, 10))
    println("First 10 Lowers: ", first(lowers, 10))

    # Clean bounds
    # clean_bounds!(uppers, lowers)

    # Validate cleaning
    # println("After Cleaning:")
    # println("First 10 Uppers: ", first(uppers, 10))
    # println("First 10 Lowers: ", first(lowers, 10))

    println("Lengths before writing -> Times: $(length(times)), Lowers: $(length(lowers)), Uppers: $(length(uppers))")

    # Save results to a CSV file
    # is_valid = validate_data_before_writing(times, lowers, uppers)

    # CSV.write("debug_times.csv", DataFrame(times=times))
    # CSV.write("debug_lowers.csv", DataFrame(lowers=lowers))
    # CSV.write("debug_uppers.csv", DataFrame(uppers=uppers))


    # println("Data is valid. Proceeding to write...")
    # Write data to CSV
    results = DataFrame(
        name = ["objective_value", "lower_bound", "gap", "runtime", "iterations", "q_liq", "q_wat", "q_oil", "q_inj", "time_nlp_solver", "time_milp_solver", "time_manipulating_model", "times", "lowers", "uppers"],
        value = [C_minlp, C_relax, gap, final_time, i, q_liq_value, q_wat_value, q_oil_value, q_inj_value, time_nlp_solver, time_milp_solver, time_model_manipulation, times, lowers, uppers]
    )
    CSV.write(filename, results, bufsize=500_000_000)

    println("\nResults saved to: $filename")
    println("Objective: ", C_minlp)
    println("Lower bound: ", C_relax)
    println("Gap: ", gap)
    println("Runtime: ", final_time)
    println("Iterations", i)
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


# SCENARIO 3
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
q_liq_max_values = nothing # Example values for q_liq_max
q_inj_max_values = [900.0] * 1e3  # Example values for q_inj_max
time_budget = 60.0 * 60 * 24  # 15 minutes

# Run tests
run_tests_with_constraints(platform, q_inj_max_values, q_liq_max_values, time_budget)

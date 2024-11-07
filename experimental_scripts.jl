using Oil
using CSV
using Printf
using Gurobi

include("models.jl")
include("utils.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)

## SCENARIO 1
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
function print_non_zero_breakpoints(vars)
    println("Non-zero breakpoints:")
    for v in vars
        if startswith(name(v), "ξ") && value(v) > 1e-6  # Filter only ξ variables with non-zero values
            @printf("Variable: %-20s Value: %f\n", name(v), value(v))
        end
    end
    println()
end

function find_infeasible_constraints_gurobi(model::Model)
    # Check if the backend optimizer is Gurobi
    if !(JuMP.optimizer_with_attributes(model) isa Gurobi.Optimizer)
        println("This function is specific to Gurobi models.")
        return
    end

    # Access the underlying Gurobi model
    gurobi_model = JuMP.backend(model).optimizer.model  # Adjusted access for Gurobi model

    # Compute IIS
    Gurobi.computeIIS(gurobi_model)

    # Print infeasible constraints in the IIS
    println("Infeasible Constraints (IIS):")
    for c in JuMP.all_constraints(model)
        constraint_name = JuMP.name(c)
        # Check if constraint is marked as part of the IIS
        if Gurobi.get_constr_attr(gurobi_model, "IISConstr", constraint_name) == 1
            println("Constraint in IIS: ", constraint_name)
        end
    end
end


function fix_and_exclude(platform::Platform, filename::String, constraint_name::String, constraint_value::Float64, time_budget::Float64)
    # Atualizar a restrição especificada na plataforma
    if constraint_name == "q_inj_max"
        platform = Platform(platform.p_sep; satellite_wells=platform.satellite_wells, 
                            manifolds=platform.manifolds,
                            q_inj_max=constraint_value)
    elseif constraint_name == "q_liq_max"
        platform = Platform(platform.p_sep; satellite_wells=platform.satellite_wells, 
                            manifolds=platform.manifolds,
                            q_liq_max=constraint_value)
    else
        error("Invalid constraint name. Use 'q_inj_max' or 'q_liq_max'")
    end

    P_minlp = get_minlp_problem(platform)
    # set_optimizer_attribute(P_minlp, "OptimalityTol", 1e-8)     # Optimality tolerance
    set_optimizer_attribute(P_minlp, "FeasibilityTol", 1e-9)    # Feasibility tolerance
    set_optimizer_attribute(P_minlp, "IntFeasTol", 1e-9)        # Integrality tolerance
    # set_optimizer_attribute(P_minlp, "BarConvTol", 1e-8)        # Barrier convergence tolerance

    C_minlp = ∞
    solution = nothing
    times = Vector{Float64}()
    uppers = Vector{Float64}()
    lowers = Vector{Float64}()

    push!(times, 0.0)
    push!(uppers, ∞)
    push!(lowers, -∞) 

    
    global start_time = time()

    ## 2. Build the MILP relaxation \tilde{P}
    P_relax = get_milp_relaxation(platform)
    # set_optimizer_attribute(P_relax, "OptimalityTol", 1e-8)     # Optimality tolerance
    set_optimizer_attribute(P_relax, "FeasibilityTol", 1e-9)    # Feasibility tolerance
    set_optimizer_attribute(P_relax, "IntFeasTol", 1e-9)        # Integrality tolerance
    # set_optimizer_attribute(P_relax, "BarConvTol", 1e-8)        # Barrier convergence tolerance
    C_relax = ∞

    ## 3. Solve the MILP relaxation \tilde{P}, get solution and cost
    C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time), save_bounds = false)
    if termination_status(P_relax) == MOI.OPTIMAL
        obj_value = objective_value(P_relax)
        push!(times, time() - start_time)
        push!(lowers, lowers[end])
        push!(uppers, obj_value)
        # print_non_zero_breakpoints(vars_relax)  # Print non-zero breakpoints here
    end

    @printf("| %-5s | %-15s | %-15s | %-15s | %-10s |\n", "It.", "Lower Bound", "Upper Bound", "Gap", "Time")
    

    i = 0
    while (C_relax < C_minlp) & (time() - start_time < time_budget)
        @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)

        # 4. Build P_fixed by constraining P to x_relax
        fix_model!(P_minlp, vars_relax)
        C_fixed = ∞
        # set_optimizer_attribute(P_minlp, "InfUnbdInfo", 1)
        # set_optimizer_attribute(P_minlp, "DualReductions", 0)


        # 5. Solve P_fixed, update x and C
        C_fixed, vars_fixed = nlp_solver(P_minlp, time_limit = time_budget - (time() - start_time), save_bounds = false)

        if termination_status(P_minlp) == MOI.OPTIMAL
            obj_value = objective_value(P_minlp)
            push!(times, time() - start_time)
            push!(lowers, obj_value)
            push!(uppers, uppers[end])
            # print_non_zero_breakpoints(vars_fixed) 
        end

        status = termination_status(P_minlp)
        # println("Model status: ", status)
    
        # if status == MOI.INFEASIBLE
        #     println("Model is infeasible. Identifying infeasible constraints:")
        #     find_infeasible_constraints_gurobi(P_minlp)
        # end

        if C_fixed < C_minlp
            C_minlp = C_fixed;
            solution = Dict(name(v) => value(v) for v in vars_fixed)
        end

        if time() - start_time >= time_budget
            break
        end

        # 6. Exclude x_relax from P_relax
        exclude!(P_relax, vars_relax)

        # 7. Solve P_relax, get solution x_relax and C_relax
        C_relax, vars_relax = milp_solver(P_relax, time_limit = time_budget - (time() - start_time), save_bounds = false)

        if termination_status(P_relax) == MOI.OPTIMAL
            obj_value = objective_value(P_relax)
            push!(times, time() - start_time)
            push!(lowers, lowers[end])
            push!(uppers, obj_value)
        end
        

        i += 1
    end
    final_time = time() - start_time

    C_relax = min(C_relax, C_minlp)  # in case the last relax solution was worse than the candidate
    @printf("| %5d | %15f | %15f | %14f%% | %9fs |\n", i, C_relax, C_minlp, 100*abs(C_minlp - C_relax) / abs(C_minlp), time()-start_time)

    # Store final results
    push!(times, final_time)
    push!(lowers, -C_minlp)
    push!(uppers, -C_relax)
    q_liq_plat = solution["q_liq_platform"]
    q_inj_plat = solution["q_inj_platform"]

    CSV.write(filename, Dict(
        "times"=>times,
        "lowers"=>lowers,
        "uppers"=>uppers,
        "q_inj_plat"=>q_inj_plat,
        "q_liq_plat"=>q_liq_plat,
        "iterations"=>i
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
end

function run_scenarios(platform::Platform, constraint_name::String, values::Vector{Float64}, time_budget::Float64)
    for value in values
        # Definir o nome do arquivo de saída com base na restrição e valor
        filename = "results_$(constraint_name)_$(value).csv"
        
        println("Running scenario for $(constraint_name) = $(value)")
        fix_and_exclude(platform, filename, constraint_name, value, time_budget)
    end
end

# q_inj_max_values = [700.0, 600.0] .*10^3 
q_liq_max_values = [4000.0]
q_inj_max_values = [6500.0] .*1e3 
time_budget = 60.0 * 60.0 * 24

# run_scenarios(platform, "q_inj_max", q_inj_max_values, time_budget)

run_scenarios(platform, "q_liq_max", q_liq_max_values, time_budget)
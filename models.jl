using Oil
using JuMP, Gurobi

const GRB_ENV = Gurobi.Env()

M_Pressure = 1000.0  # Big-M


function add_custom_SOS2_constraint(model::GenericModel, ξ; name = "sos2")
    z = @variable(model, [1:length(ξ)-1], Bin, base_name="z_$name")

    @constraint(model, sum(z) == 1, base_name="single_interval_$name")
    @constraint(model, sum(ξ) == 1, base_name="convex_combination_ξ_$name")

    # Constraints for IGLR
    for (i, key_tuple) in enumerate(keys(ξ))
        if i == 1
            # Edge constraint for first value
            @constraint(model, ξ[key_tuple] <= z[i], base_name="link_ξ_$(i)_$name")
        elseif i == length(keys(ξ))
            # Edge constraint for last value
            @constraint(model, ξ[key_tuple] <= z[i-1], base_name="link_ξ_$(i)_$name")
        else
            # Middle constraint for intermediate values
            @constraint(model, ξ[key_tuple] <= z[i-1] + z[i], base_name="link_ξ_$(i)_$name")
        end
    end

    return z
end

# TODO: convert these functions into macros
function add_nonlinear_well(model::GenericModel, well::Oil.AbstractWell; sos2_with_binary = false)
    # preprocess VLP curve for a piecewise formulation
    well = Oil.PiecewiseLinearWell(well, true)

    n = well.name

    # Variables

    ## Flow rates
    q_liq_n = @variable(model, base_name="q_liq_$n", lower_bound = 0)
    q_water_n = @variable(model, base_name="q_water_$n", lower_bound = 0)
    q_oil_n = @variable(model, base_name="q_oil_$n", lower_bound = 0)
    q_gas_n = @variable(model, base_name="q_gas_$n", lower_bound = 0)
    q_inj_n = @variable(model, base_name="q_inj_$n", lower_bound = 0)

    ## Pressures
    whp_n = @variable(model, base_name="whp_$n", lower_bound = 0)  # wellhead
    bhp_n = @variable(model, base_name="bhp_$n", lower_bound = 0)  # bottomhole
    wfp_n = @variable(model, base_name="wfp_$n", lower_bound = 0)  # well-flowing
    dsp_n = @variable(model, base_name="dsp_$n", lower_bound = 0)  # downstream

    ## Parameters
    iglr_n = @variable(model, base_name="iglr_$n", lower_bound = 0)  # injected-gas to liquid ratio
    wct_n = well.wct
    gor_n = well.gor

    ## Decisions
    y_n = @variable(model, base_name="y_$n", binary = true)        # (62h) well activation
    t_gl_n = @variable(model, base_name="t_gl_$n", binary = true)  # (62u) well lifting

    ## VLP curve variables
    λ_vlp_n = @variable(
        model,
        [well.IGLR, well.WHP, well.Q_liq_vlp],
        lower_bound=0.0,  # (62i)
        base_name="λ_vlp_$n"
    )
    ξ_iglr_n = @variable(model, [well.IGLR], lower_bound=0.0, base_name="ξ_iglr_$n")
    ξ_whp_n = @variable(model, [well.WHP], lower_bound=0.0, base_name="ξ_whp_$n")
    ξ_qliq_n = @variable(model, [well.Q_liq_vlp], lower_bound=0.0, base_name="ξ_qliq_vlp_$n")

    # Constraints

    ## VLP

    ### (62a)
    @constraint(model, iglr_n == sum(
        iglr_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name="pwl_iglr_$n")

    ### (62b)
    @constraint(model, whp_n == sum(
        whp_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_whp_$n")

    ### (62c)
    @constraint(model, q_liq_n == sum(
        q_liq_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_q_liq_vlp_$n")

    ### (62d)
    # @NonlinearConstraint
    @constraint(model, q_inj_n == q_liq_n * iglr_n, base_name = "bilinear_q_inj_$n")

    ### (62e)
    @constraint(model, wfp_n == sum(
        well.WFP_vlp[iglr_bp, whp_bp, q_liq_bp] * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_bhp_vlp_lb_$n")

    ### (62f)
    @constraint(model, whp_n >= dsp_n - M_Pressure * (1 - y_n), base_name = "well_choke_dP_$n")

    ### (62g)
    @constraint(model, y_n == sum(  # SELECTION
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_vlp_curve_$n")

    ### (62j)
    @constraint(model, [iglr_bp in well.IGLR], ξ_iglr_n[iglr_bp] == sum(
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "eta_iglr_$n")

    ### (62k)
    @constraint(model, [whp_bp in well.WHP], ξ_whp_n[whp_bp] == sum(
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "eta_whp_$n")

    ### (62l)
    @constraint(model, [q_liq_bp in well.Q_liq_vlp], ξ_qliq_n[q_liq_bp] == sum(
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
    ), base_name = "eta_q_liq_vlp_$n")

    if !sos2_with_binary
        ### (62m)
        @constraint(model, ξ_iglr_n[well.IGLR] in SOS2(well.IGLR), base_name = "sos2_iglr_$n")

        ### (62n)
        @constraint(model, ξ_whp_n[well.WHP] in SOS2(well.WHP), base_name = "sos2_whp_$n")

        ### (62o)
        @constraint(model, ξ_qliq_n[well.Q_liq_vlp] in SOS2(well.Q_liq_vlp), base_name = "sos2_qliq_vlp_$n")
    else
        ### (62m)
        z_iglr_n = add_custom_SOS2_constraint(model, ξ_iglr_n, name="iglr_$n")

        ### (62n)
        z_whp_n = add_custom_SOS2_constraint(model, ξ_whp_n, name="whp_$n")

        ### (62o)
        z_qliq_n = add_custom_SOS2_constraint(model, ξ_qliq_n, name="qliq_$n")
    end

    ### (62p)
    @constraint(model, q_oil_n == q_liq_n * (1 - wct_n), base_name = "ratio_q_oil_$n")

    ### (62q)
    @constraint(model, q_water_n == q_liq_n * wct_n, base_name = "ratio_q_water_$n")

    ### (62r)
    @constraint(model, q_gas_n == q_oil_n * gor_n, base_name = "ratio_q_gas_$n")

    ### (62s)
    if ~isnothing(well.min_q_inj)
        @constraint(model, q_inj_n >= well.min_q_inj * t_gl_n, base_name="min_q_inj_$n")
    end
    if ~isnothing(well.max_q_inj)
        @constraint(model, q_inj_n <= well.max_q_inj * t_gl_n, base_name="max_q_inj_$n")
    else
        @constraint(model, q_inj_n <= M_qinj * t_gl_n, base_name="max_q_inj_$n")
    end

    ### (62t)
    @constraint(model, t_gl_n <= y_n, base_name = "gas_lifted_$n")  # closed wells are "surgent"

    ### Multilinear Interpolation Constraints
    bilinear_left = @variable(model, [well.Q_liq_vlp, well.IGLR], base_name="bilinear_left_$n")
    # @NonlinearConstraint
    @constraint(model, [q_liq_bp in well.Q_liq_vlp, iglr_bp in well.IGLR],
        bilinear_left[q_liq_bp, iglr_bp] == ξ_qliq_n[q_liq_bp] * ξ_iglr_n[iglr_bp],
        base_name = "bilinear_left_constraint_$n"
    )
    # @NonlinearConstraint
    @constraint(model, [q_liq_bp in well.Q_liq_vlp, iglr_bp in well.IGLR, whp_bp in well.WHP],
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp] == bilinear_left[q_liq_bp, iglr_bp] * ξ_whp_n[whp_bp],
        base_name = "well_bilinear_constraint_$n"
    )

    ## IPR

    ### LinearIPR
    @constraint(model, q_liq_n == well.ipr.IP * (well.ipr.p_res * y_n - bhp_n),
                base_name = "linear_ipr_$n")

    ## Infeasible Points

    for iglr_bp in well.IGLR
        for whp_bp in well.WHP
            for q_liq_bp in well.Q_liq_vlp
                if well.WFP_vlp[iglr_bp, whp_bp, q_liq_bp] == -1
                    @constraint(
                        model,
                        # infeasible_bps_vlp[n,iglr_bp,whp_bp,q_liq_bp],
                        0 == λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
                    )
                end
            end
        end
    end

    return model
end

function add_nonlinear_manifold(model::GenericModel, manifold::Manifold, p_sep::Float64; sos2_with_binary = false)
    N_man = names(manifold.wells)  # names of all wells in manifold

    m = manifold.name

    for well in manifold.wells
        model = add_nonlinear_well(model, well, sos2_with_binary = sos2_with_binary)
    end

    manifold = Oil.PiecewiseLinearManifold(manifold, p_sep)

    # Variables

    ## Flow rates
    q_liq_m = @variable(model, base_name="q_liq_$m", lower_bound = 0)
    q_oil_m = @variable(model, base_name="q_oil_$m", lower_bound = 0)
    q_water_m = @variable(model, base_name="q_water_$m", lower_bound = 0)
    q_gas_m = @variable(model, base_name="q_gas_$m", lower_bound = 0)
    q_inj_m = @variable(model, base_name="q_inj_$m", lower_bound = 0)

    ## Pressures
    dsp_m = @variable(model, base_name="dsp_$m", lower_bound = 0)  # manifold downstream pressure
    p_sup_m = @variable(model, base_name="p_sup_$m", lower_bound = 0)  # manifold choke pressure
    ΔP = @variable(model, base_name="ΔP_$m", lower_bound = 0)  # riser pressure drop

    ## Parameters
    gor_m = @variable(model, base_name="gor_$m", lower_bound=0.0)
    wct_m = @variable(model, base_name="wct_$m", lower_bound=0.0)
    iglr_m = @variable(model, base_name="iglr_$m", lower_bound=0.0)

    ## Discrete decisions
    y_m = @variable(model, base_name="y_$m", binary = true)  # (65d)

    ## VLP curve variables
    λ_riser_m = @variable(model, [manifold.Q_liq, manifold.GOR, manifold.WCT, manifold.IGLR],
                          lower_bound=0.0, base_name="λ_riser_$m")  # (65d)
    ξ_qliq_m = @variable(model, [manifold.Q_liq], lower_bound=0.0, base_name="ξ_qliq_$m")
    ξ_gor_m = @variable(model, [manifold.GOR], lower_bound=0.0, base_name="ξ_gor_$m")
    ξ_wct_m = @variable(model, [manifold.WCT], lower_bound=0.0, base_name="ξ_wct_$m")
    ξ_iglr_m = @variable(model, [manifold.IGLR], lower_bound=0.0, base_name="ξ_iglr_$m")

    # Constraints

    ## Manifold input

    ### (65a)
    @constraint(model, q_liq_m == sum(variable_by_name(model, "q_liq_$n") for n in N_man), base_name="liquid_rate_$m")
    @constraint(model, q_oil_m == sum(variable_by_name(model, "q_oil_$n") for n in N_man), base_name="oil_rate_$m")
    @constraint(model, q_water_m == sum(variable_by_name(model, "q_water_$n") for n in N_man), base_name="water_rate_$m")
    @constraint(model, q_gas_m == sum(variable_by_name(model, "q_gas_$n") for n in N_man), base_name="gas_rate_$m")
    @constraint(model, q_inj_m == sum(variable_by_name(model, "q_inj_$n") for n in N_man), base_name="inj_gas_rate_$m")

    ## VLP

    ### (65b)
    @constraint(model, q_liq_m == sum(
        q_liq_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_qliq_riser_$m")
    @constraint(model, gor_m == sum(
        gor_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_gor_riser_$m")
    @constraint(model, wct_m == sum(
        wct_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_wct_riser_$m")
    @constraint(model, iglr_m == sum(
        iglr_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_iglr_riser_$m")
    @constraint(model, p_sup_m <= dsp_m + M_Pressure * (1 - y_m), base_name="p_sup_upper_$m")
    @constraint(model, p_sup_m >= dsp_m - M_Pressure * (1 - y_m), base_name="riser_choke_dP_$m")

    ### (65c)
    @constraint(model, q_oil_m == q_liq_m - q_water_m, base_name = "bilinear_q_oil_riser_$m")
    # @NonlinearConstraint
    @constraint(model, q_water_m == q_liq_m * wct_m, base_name = "bilinear_q_water_riser_$m")
    # @NonlinearConstraint
    @constraint(model, q_gas_m == q_oil_m * gor_m, base_name = "bilinear_q_gas_riser_$m")
    # @NonlinearConstraint
    @constraint(model, q_inj_m == q_liq_m * iglr_m, base_name = "bilinear_q_inj_riser_$m")

    @constraint(model, ΔP == sum(
        manifold.ΔP[q_liq_bp, gor_bp, wct_bp, iglr_bp] * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_ΔP_riser_$m")

    ### (65d)
    @constraint(model, y_m == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_riser_curve_$m")
    @constraint(model, [n in N_man], y_m >= variable_by_name(model, "y_$n"), base_name="forced_activation_riser_$m")
    @constraint(model, y_m <= sum(variable_by_name(model, "y_$n") for n in N_man), base_name="forced_deactivation_riser_$m")

    ### (65e)
    @constraint(model, [q_liq_bp in manifold.Q_liq], ξ_qliq_m[q_liq_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="eta_qliq_riser_$m")
    @constraint(model, [gor_bp in manifold.GOR], ξ_gor_m[gor_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="eta_gor_riser_$m")
    @constraint(model, [wct_bp in manifold.WCT], ξ_wct_m[wct_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for iglr_bp in manifold.IGLR
    ), base_name="eta_wct_riser_$m")
    @constraint(model, [iglr_bp in manifold.IGLR], ξ_iglr_m[iglr_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
    ), base_name="eta_iglr_riser_$m")

    if !sos2_with_binary
        @constraint(model, ξ_qliq_m[manifold.Q_liq] in SOS2(manifold.Q_liq), base_name="sos2_qliq_riser_$m")
        @constraint(model, ξ_gor_m[manifold.GOR] in SOS2(manifold.GOR), base_name="sos2_gor_riser_$m")
        @constraint(model, ξ_wct_m[manifold.WCT] in SOS2(manifold.WCT), base_name="sos2_wct_riser_$m")
        @constraint(model, ξ_iglr_m[manifold.IGLR] in SOS2(manifold.IGLR), base_name="sos2_iglr_riser_$m")
    else
        z_qliq_m = add_custom_SOS2_constraint(model, ξ_qliq_m, name="qliq_$m")
        z_gor_m = add_custom_SOS2_constraint(model, ξ_gor_m, name="gor_$m")
        z_wct_m = add_custom_SOS2_constraint(model, ξ_wct_m, name="wct_$m")
        z_iglr_m = add_custom_SOS2_constraint(model, ξ_iglr_m, name="iglr_$m")
    end

    ## Multilinear Interpolation

    bilinear_left = @variable(model, [manifold.Q_liq, manifold.GOR], base_name="bilinear_left_$m")
    bilinear_right = @variable(model, [manifold.WCT, manifold.IGLR], base_name="bilinear_right_$m")

    # @NonlinearConstraint
    @constraint(model, [q_liq_bp in manifold.Q_liq, gor_bp in manifold.GOR],
        bilinear_left[q_liq_bp, gor_bp] == ξ_qliq_m[q_liq_bp] * ξ_gor_m[gor_bp],
    base_name="bilinear_left_riser_$m")
    # @NonlinearConstraint
    @constraint(model, [wct_bp in manifold.WCT, iglr_bp in manifold.IGLR],
        bilinear_right[wct_bp, iglr_bp] == ξ_wct_m[wct_bp] * ξ_iglr_m[iglr_bp],
    base_name="riser_bilinear_right_$m")
    # @NonlinearConstraint
    @constraint(model, [q_liq_bp in manifold.Q_liq, gor_bp in manifold.GOR, wct_bp in manifold.WCT, iglr_bp in manifold.IGLR],
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp] == bilinear_left[q_liq_bp,gor_bp] * bilinear_right[wct_bp,iglr_bp],
    base_name="riser_bilinear_$m")

    ## Infeasible points
    for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
            for wct_bp in manifold.WCT
                for iglr_bp in manifold.IGLR
                    if manifold.ΔP[(q_liq_bp, gor_bp, wct_bp, iglr_bp)] < 0
                        @constraint(
                            model,
                            # infeasible_bps_vlp[n,iglr_bp,whp_bp,q_liq_bp],
                            0 == λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
                        )
                    end
                end
            end
        end
    end

    return model
end

function add_platform(model::GenericModel, platform::Platform)
    N_sat = names(platform.satellite_wells)

    N = copy(N_sat)
    for manifold in platform.manifolds
        N = [N ; names(manifold.wells)]
    end

    # Constraints

    ### (64b)/(66b)
    @constraint(
        model, vlp_ipr_match[n in N],
        variable_by_name(model, "bhp_$n") == variable_by_name(model, "wfp_$n")
    )

    ### (64c)/(66c)
    @constraint(
        model, satellite_wells_p_sep[n in N_sat],
        variable_by_name(model, "dsp_$n") == platform.p_sep
    )

    # the following constraints suit the add_manifold_riser function better
    for manifold in platform.manifolds
        m = manifold.name

        ### (66d)
        @constraint(
            model, [n in names(manifold.wells)],
            variable_by_name(model, "dsp_$n") == variable_by_name(model, "p_sup_$m")
                                                 + variable_by_name(model, "ΔP_$m"),
            base_name="manifold_wells_downstream_p_$m"
        )

        ### (66e)
        @constraint(
            model, base_name = "downstream_pressure_manifold_$m",
            variable_by_name(model, "dsp_$m") == platform.p_sep
        )
    end

    ### (64d)/(66f)
    if ~isnothing(platform.q_inj_max)
        @constraint(
            model, gl_injection_limit,
            platform.q_inj_max >= sum(variable_by_name(model, "q_inj_$n") for n in N)
        )
    end

    ### (64e)/(66g)
    if ~isnothing(platform.q_water_max)
        @constraint(
            model, water_production_limit,
            platform.q_water_max >= sum(variable_by_name(model, "q_water_$n") for n in N)
        )
    end

    ### (64f)/(66h)
    if ~isnothing(platform.q_gas_max)
        @constraint(
            model, gas_capacity,
            platform.q_gas_max >= sum(
                variable_by_name(model, "q_gas_$n") + variable_by_name(model, "q_inj_$n")
                for n in N
            )
        )
    end

    ### (64g)/(66i)
    if ~isnothing(platform.q_liq_max)
        @constraint(
            model, liquid_production_limit,
            platform.q_liq_max >= sum(variable_by_name(model, "q_liq_$n") for n in N)
        )
    end

    return model
end

function get_minlp_problem(platform::Platform; sos2_with_binary = false)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))

    # Objective
    @variable(model, q_oil_total >= 0)
    q_oil_total_expr = 0  # modeled like this to facilitate LBD valid inequality

    # Wells
    for well in platform.satellite_wells
        model = add_nonlinear_well(model, well, sos2_with_binary = sos2_with_binary)
        q_oil_total_expr = q_oil_total_expr + variable_by_name(model, "q_oil_$(well.name)")
    end

    # Manifolds
    for manifold in platform.manifolds
        model = add_nonlinear_manifold(model, manifold, platform.p_sep, sos2_with_binary = sos2_with_binary)  # (65)
        q_oil_total_expr = q_oil_total_expr + variable_by_name(model, "q_oil_$(manifold.name)")
    end

    # Platform
    model = add_platform(model, platform)

    # Objective
    @constraint(model, q_oil_total == q_oil_total_expr)
    @objective(model, Max, q_oil_total)  # (64a)/(66a)

    return model
end

function add_linear_well(model::GenericModel, well::Oil.AbstractWell; sos2_with_binary = false)
    # preprocess VLP curve for a piecewise formulation
    well = Oil.PiecewiseLinearWell(well, true)

    n = well.name

    # Variables

    ## Flow rates
    q_liq_n = @variable(model, base_name="q_liq_$n", lower_bound = 0)
    q_water_n = @variable(model, base_name="q_water_$n", lower_bound = 0)
    q_oil_n = @variable(model, base_name="q_oil_$n", lower_bound = 0)
    q_gas_n = @variable(model, base_name="q_gas_$n", lower_bound = 0)
    q_inj_n = @variable(model, base_name="q_inj_$n", lower_bound = 0)

    ## Pressures
    whp_n = @variable(model, base_name="whp_$n", lower_bound = 0)  # wellhead
    bhp_n = @variable(model, base_name="bhp_$n", lower_bound = 0)  # bottomhole
    wfp_n = @variable(model, base_name="wfp_$n", lower_bound = 0)  # well-flowing
    dsp_n = @variable(model, base_name="dsp_$n", lower_bound = 0)  # downstream

    ## Parameters
    iglr_n = @variable(model, base_name="iglr_$n", lower_bound = 0)  # injected-gas to liquid ratio
    wct_n = well.wct
    gor_n = well.gor

    ## Decisions
    y_n = @variable(model, base_name="y_$n", binary = true)        # (62h) well activation
    t_gl_n = @variable(model, base_name="t_gl_$n", binary = true)  # (62u) well lifting

    ## VLP curve variables
    λ_vlp_n = @variable(
        model,
        [well.IGLR, well.WHP, well.Q_liq_vlp],
        lower_bound=0.0,  # (62i)
        base_name="λ_vlp_$n"
    )
    ξ_iglr_n = @variable(model, [well.IGLR], lower_bound=0.0, base_name="ξ_iglr_$n")
    ξ_whp_n = @variable(model, [well.WHP], lower_bound=0.0, base_name="ξ_whp_$n")
    ξ_qliq_n = @variable(model, [well.Q_liq_vlp], lower_bound=0.0, base_name="ξ_qliq_vlp_$n")

    # Constraints

    ## VLP

    ### (62a)
    @constraint(model, iglr_n == sum(
        iglr_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name="pwl_iglr_$n")

    ### (62b)
    @constraint(model, whp_n == sum(
        whp_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_whp_$n")

    ### (62c)
    @constraint(model, q_liq_n == sum(
        q_liq_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_q_liq_vlp_$n")

    ### (62d)
    # @PiecewiseRelaxation
    @constraint(model, q_inj_n == sum(
        q_liq_bp * iglr_bp * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_q_inj_$n")

    ### (62e)
    @constraint(model, wfp_n == sum(
        well.WFP_vlp[iglr_bp, whp_bp, q_liq_bp] * λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_bhp_vlp_lb_$n")

    ### (62f)
    @constraint(model, whp_n >= dsp_n - M_Pressure * (1 - y_n), base_name = "well_choke_dP_$n")

    ### (62g)
    @constraint(model, y_n == sum(  # SELECTION
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "pwl_vlp_curve_$n")

    ### (62j)
    @constraint(model, [iglr_bp in well.IGLR], ξ_iglr_n[iglr_bp] == sum(
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for whp_bp in well.WHP
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "eta_iglr_$n")

    ### (62k)
    @constraint(model, [whp_bp in well.WHP], ξ_whp_n[whp_bp] == sum(
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for q_liq_bp in well.Q_liq_vlp
    ), base_name = "eta_whp_$n")

    ### (62l)
    @constraint(model, [q_liq_bp in well.Q_liq_vlp], ξ_qliq_n[q_liq_bp] == sum(
        λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
        for iglr_bp in well.IGLR
        for whp_bp in well.WHP
    ), base_name = "eta_q_liq_vlp_$n")

    if !sos2_with_binary
        ### (62m)
        @constraint(model, ξ_iglr_n[well.IGLR] in SOS2(well.IGLR), base_name = "sos2_iglr_$n")

        ### (62n)
        @constraint(model, ξ_whp_n[well.WHP] in SOS2(well.WHP), base_name = "sos2_whp_$n")

        ### (62o)
        @constraint(model, ξ_qliq_n[well.Q_liq_vlp] in SOS2(well.Q_liq_vlp), base_name = "sos2_qliq_vlp_$n")
    else
        ### (62m)
        z_iglr_n = add_custom_SOS2_constraint(model, ξ_iglr_n, name="iglr_$n")

        ### (62n)
        z_whp_n = add_custom_SOS2_constraint(model, ξ_whp_n, name="whp_$n")

        ### (62o)
        z_qliq_n = add_custom_SOS2_constraint(model, ξ_qliq_n, name="qliq_$n")
    end

    ### (62p)
    @constraint(model, q_oil_n == q_liq_n * (1 - wct_n), base_name = "ratio_q_oil_$n")

    ### (62q)
    @constraint(model, q_water_n == q_liq_n * wct_n, base_name = "ratio_q_water_$n")

    ### (62r)
    @constraint(model, q_gas_n == q_oil_n * gor_n, base_name = "ratio_q_gas_$n")

    ### (62s)
    if ~isnothing(well.min_q_inj)
        @constraint(model, q_inj_n >= well.min_q_inj * t_gl_n, base_name="min_q_inj_$n")
    end
    if ~isnothing(well.max_q_inj)
        @constraint(model, q_inj_n <= well.max_q_inj * t_gl_n, base_name="max_q_inj_$n")
    else
        @constraint(model, q_inj_n <= M_qinj * t_gl_n, base_name="max_q_inj_$n")
    end

    ### (62t)
    @constraint(model, t_gl_n <= y_n, base_name = "gas_lifted_$n")  # closed wells are "surgent"

    ## IPR

    ### LinearIPR
    @constraint(model, q_liq_n == well.ipr.IP * (well.ipr.p_res * y_n - bhp_n),
                base_name = "linear_ipr_$n")

    ## Infeasible Points

    for iglr_bp in well.IGLR
        for whp_bp in well.WHP
            for q_liq_bp in well.Q_liq_vlp
                if well.WFP_vlp[iglr_bp, whp_bp, q_liq_bp] == -1
                    @constraint(
                        model,
                        # infeasible_bps_vlp[n,iglr_bp,whp_bp,q_liq_bp],
                        0 == λ_vlp_n[iglr_bp, whp_bp, q_liq_bp]
                    )
                end
            end
        end
    end

    return model
end

function add_linear_manifold(model::GenericModel, manifold::Manifold, p_sep::Float64; sos2_with_binary = false)
    N_man = names(manifold.wells)  # names of all wells in manifold

    m = manifold.name

    for well in manifold.wells
        model = add_linear_well(model, well, sos2_with_binary = sos2_with_binary)
    end

    manifold = Oil.PiecewiseLinearManifold(manifold, p_sep)

    # Variables

    ## Flow rates
    q_liq_m = @variable(model, base_name="q_liq_$m", lower_bound = 0)
    q_oil_m = @variable(model, base_name="q_oil_$m", lower_bound = 0)
    q_water_m = @variable(model, base_name="q_water_$m", lower_bound = 0)
    q_gas_m = @variable(model, base_name="q_gas_$m", lower_bound = 0)
    q_inj_m = @variable(model, base_name="q_inj_$m", lower_bound = 0)

    ## Pressures
    dsp_m = @variable(model, base_name="dsp_$m", lower_bound = 0)  # manifold downstream pressure
    p_sup_m = @variable(model, base_name="p_sup_$m", lower_bound = 0)  # manifold choke pressure
    ΔP = @variable(model, base_name="ΔP_$m", lower_bound = 0)  # riser pressure drop

    ## Parameters
    gor_m = @variable(model, base_name="gor_$m")
    wct_m = @variable(model, base_name="wct_$m")
    iglr_m = @variable(model, base_name="iglr_$m")

    ## Discrete decisions
    y_m = @variable(model, base_name="y_$m", binary = true)  # (65d)

    ## VLP curve variables
    λ_riser_m = @variable(model, [manifold.Q_liq, manifold.GOR, manifold.WCT, manifold.IGLR],
                          lower_bound=0.0, base_name="λ_riser_$m")  # (65d)
    ξ_qliq_m = @variable(model, [manifold.Q_liq], lower_bound=0.0, base_name="ξ_qliq_$m")
    ξ_gor_m = @variable(model, [manifold.GOR], lower_bound=0.0, base_name="ξ_gor_$m")
    ξ_wct_m = @variable(model, [manifold.WCT], lower_bound=0.0, base_name="ξ_wct_$m")
    ξ_iglr_m = @variable(model, [manifold.IGLR], lower_bound=0.0, base_name="ξ_iglr_$m")

    # Constraints

    ## Manifold input

    ### (65a)
    @constraint(model, q_liq_m == sum(variable_by_name(model, "q_liq_$n") for n in N_man), base_name="manifold_liquid_rate_$m")
    @constraint(model, q_oil_m == sum(variable_by_name(model, "q_oil_$n") for n in N_man), base_name="manifold_oil_rate_$m")
    @constraint(model, q_water_m == sum(variable_by_name(model, "q_water_$n") for n in N_man), base_name="manifold_water_rate_$m")
    @constraint(model, q_gas_m == sum(variable_by_name(model, "q_gas_$n") for n in N_man), base_name="manifold_gas_rate_$m")
    @constraint(model, q_inj_m == sum(variable_by_name(model, "q_inj_$n") for n in N_man), base_name="manifold_inj_gas_rate_$m")

    ## VLP

    ### (65b)
    @constraint(model, q_liq_m == sum(
        q_liq_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_qliq_riser_$m")
    @constraint(model, gor_m == sum(
        gor_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_gor_riser_$m")
    @constraint(model, wct_m == sum(
        wct_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_wct_riser_$m")
    @constraint(model, iglr_m == sum(
        iglr_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_iglr_riser_$m")
    @constraint(model, p_sup_m <= dsp_m + M_Pressure * (1 - y_m), base_name="p_sup_upper_$m")
    @constraint(model, p_sup_m >= dsp_m - M_Pressure * (1 - y_m), base_name="riser_choke_dP_$m")

    ### (65c)
    @constraint(model, q_oil_m == q_liq_m - q_water_m, base_name = "bilinear_q_oil_riser_$m")
    # @PiecewiseRelaxation
    @constraint(model, q_oil_m == sum(
        q_liq_bp * (1 - wct_bp) * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_qoil_riser_$m")
    # @PiecewiseRelaxation
    @constraint(model, q_water_m == sum(
        q_liq_bp * wct_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_qwater_riser_$m")
    # @PiecewiseRelaxation
    @constraint(model, q_gas_m == sum(
        q_liq_bp * (1 - wct_bp) * gor_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_qgas_riser_$m")
    # @PiecewiseRelaxation
    @constraint(model, q_inj_m == sum(
        q_liq_bp * iglr_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_qinj_riser_$m")

    @constraint(model, ΔP == sum(
        manifold.ΔP[q_liq_bp, gor_bp, wct_bp, iglr_bp] * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_ΔP_riser_$m")

    ### (65d)
    @constraint(model, y_m == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="pwl_riser_curve_$m")
    @constraint(model, [n in N_man], y_m >= variable_by_name(model, "y_$n"), base_name="forced_activation_riser_$m")
    @constraint(model, y_m <= sum(variable_by_name(model, "y_$n") for n in N_man), base_name="forced_deactivation_riser_$m")

    ### (65e)
    @constraint(model, [q_liq_bp in manifold.Q_liq], ξ_qliq_m[q_liq_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="eta_qliq_riser_$m")
    @constraint(model, [gor_bp in manifold.GOR], ξ_gor_m[gor_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for wct_bp in manifold.WCT
        for iglr_bp in manifold.IGLR
    ), base_name="eta_gor_riser_$m")
    @constraint(model, [wct_bp in manifold.WCT], ξ_wct_m[wct_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for iglr_bp in manifold.IGLR
    ), base_name="eta_wct_riser_$m")
    @constraint(model, [iglr_bp in manifold.IGLR], ξ_iglr_m[iglr_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
        for wct_bp in manifold.WCT
    ), base_name="eta_iglr_riser_$m")

    if !sos2_with_binary
        @constraint(model, ξ_qliq_m[manifold.Q_liq] in SOS2(manifold.Q_liq), base_name="sos2_qliq_riser_$m")
        @constraint(model, ξ_gor_m[manifold.GOR] in SOS2(manifold.GOR), base_name="sos2_gor_riser_$m")
        @constraint(model, ξ_wct_m[manifold.WCT] in SOS2(manifold.WCT), base_name="sos2_wct_riser_$m")
        @constraint(model, ξ_iglr_m[manifold.IGLR] in SOS2(manifold.IGLR), base_name="sos2_iglr_riser_$m")
    else
        z_qliq_m = add_custom_SOS2_constraint(model, ξ_qliq_m, name="qliq_$m")
        z_gor_m = add_custom_SOS2_constraint(model, ξ_gor_m, name="gor_$m")
        z_wct_m = add_custom_SOS2_constraint(model, ξ_wct_m, name="wct_$m")
        z_iglr_m = add_custom_SOS2_constraint(model, ξ_iglr_m, name="iglr_$m")
    end

    ## Infeasible points
    for q_liq_bp in manifold.Q_liq
        for gor_bp in manifold.GOR
            for wct_bp in manifold.WCT
                for iglr_bp in manifold.IGLR
                    if manifold.ΔP[(q_liq_bp, gor_bp, wct_bp, iglr_bp)] < 0
                        @constraint(
                            model,
                            # infeasible_bps_vlp[n,iglr_bp,whp_bp,q_liq_bp],
                            0 == λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
                        )
                    end
                end
            end
        end
    end

    return model
end

function get_milp_relaxation(platform::Platform; sos2_with_binary = false)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))

    # Objective
    @variable(model, q_oil_total >= 0)
    q_oil_total_expr = 0  # modeled like this to facilitate LBD valid inequality

    # Wells
    for well in platform.satellite_wells
        model = add_linear_well(model, well, sos2_with_binary = sos2_with_binary)
        q_oil_total_expr = q_oil_total_expr + variable_by_name(model, "q_oil_$(well.name)")
    end

    # Manifolds
    for manifold in platform.manifolds
        model = add_linear_manifold(model, manifold, platform.p_sep, sos2_with_binary = sos2_with_binary)
        q_oil_total_expr = q_oil_total_expr + variable_by_name(model, "q_oil_$(manifold.name)")
    end

    # Platform
    model = add_platform(model, platform)

    # Objective
    @constraint(model, q_oil_total == q_oil_total_expr)
    @objective(model, Max, q_oil_total)  # (64a)/(66a)

    return model
end

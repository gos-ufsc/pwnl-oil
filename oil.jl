using Oil
using JuMP, SCIP

M_Pressure = 1000.0  # Big-M

# TODO: convert these functions into macros
function add_nonlinear_well(model::GenericModel, well::Oil.AbstractWell)
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

    ### (62m)
    @constraint(model, ξ_iglr_n[well.IGLR] in SOS2(well.IGLR), base_name = "sos2_iglr_$n")

    ### (62n)
    @constraint(model, ξ_whp_n[well.WHP] in SOS2(well.WHP), base_name = "sos2_whp_$n")

    ### (62o)
    @constraint(model, ξ_qliq_n[well.Q_liq_vlp] in SOS2(well.Q_liq_vlp), base_name = "sos2_qliq_vlp_$n")

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

function add_nonlinear_riser(model::GenericModel, riser::Riser, p_sep::Float64)
    N_man = names(riser.manifold_wells)  # names of all wells in manifold

    for well in riser.manifold_wells
        model = add_nonlinear_well(model, well)
    end

    riser = Oil.PiecewiseLinearRiser(riser, p_sep, false)

    # Variables

    ## Flow rates
    @variables(model, begin
        q_liq_m >= 0.0
        q_oil_m >= 0.0
        q_water_m >= 0.0
        q_gas_m >= 0.0
        q_inj_m >= 0.0
    end)

    ## Pressures
    @variables(model, begin
        dsp_m   >= 0.0  # manifold downstream pressure
        p_sup_m >= 0.0  # manifold choke pressure
        ΔP      >= 0.0  # riser pressure drop
    end)

    ## Parameters
    @variables(model, begin
        gor_m
        wct_m
        iglr_m
    end)

    ## Discrete decisions
    @variable(model, y_m, Bin)  # (65d)

    ## VLP curve variables
    λ_riser_m = @variable(model, [riser.Q_liq, riser.GOR, riser.WCT, riser.IGLR],
                          lower_bound=0.0, base_name="λ_riser_m")  # (65d)
    ξ_qliq_m = @variable(model, [riser.Q_liq], lower_bound=0.0, base_name="ξ_qliq_m")
    ξ_gor_m = @variable(model, [riser.GOR], lower_bound=0.0, base_name="ξ_gor_m")
    ξ_wct_m = @variable(model, [riser.WCT], lower_bound=0.0, base_name="ξ_wct_m")
    ξ_iglr_m = @variable(model, [riser.IGLR], lower_bound=0.0, base_name="ξ_iglr_m")

    # Constraints

    ## Manifold input

    ### (65a)
    @constraint(model, manifold_liquid_rate, q_liq_m == sum(variable_by_name(model, "q_liq_$n") for n in N_man))
    @constraint(model, manifold_oil_rate, q_oil_m == sum(variable_by_name(model, "q_oil_$n") for n in N_man))
    @constraint(model, manifold_water_rate, q_water_m == sum(variable_by_name(model, "q_water_$n") for n in N_man))
    @constraint(model, manifold_gas_rate, q_gas_m == sum(variable_by_name(model, "q_gas_$n") for n in N_man))
    @constraint(model, manifold_inj_gas_rate, q_inj_m == sum(variable_by_name(model, "q_inj_$n") for n in N_man))

    ## VLP

    ### (65b)
    @constraint(model, pwl_qliq_riser_m, q_liq_m == sum(
        q_liq_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, pwl_gor_riser_m, gor_m == sum(
        gor_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, pwl_wct_riser_m, wct_m == sum(
        wct_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, pwl_iglr_riser_m, iglr_m == sum(
        iglr_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, p_sup_upper, p_sup_m <= dsp_m + M_Pressure * (1 - y_m))
    @constraint(model, riser_choke_dP, p_sup_m >= dsp_m - M_Pressure * (1 - y_m))

    ### (65c)
    @constraint(model, q_oil_m == q_liq_m - q_water_m, base_name = "bilinear_q_oil_riser_m")
    # @NonlinearConstraint
    @constraint(model, q_water_m == q_liq_m * wct_m, base_name = "bilinear_q_water_riser_m")
    # @NonlinearConstraint
    @constraint(model, q_gas_m == q_oil_m * gor_m, base_name = "bilinear_q_gas_riser_m")
    # @NonlinearConstraint
    @constraint(model, q_inj_m == q_liq_m * iglr_m, base_name = "bilinear_q_inj_riser_m")

    @constraint(model, pwl_ΔP_riser_m, ΔP == sum(
        riser.ΔP[q_liq_bp, gor_bp, wct_bp, iglr_bp] * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))

    ### (65d)
    @constraint(model, pwl_riser_curve, y_m == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, forced_activation_riser_m[n in N_man], y_m >= variable_by_name(model, "y_$n"))
    @constraint(model, forced_deactivation_riser_m, y_m <= sum(variable_by_name(model, "y_$n") for n in N_man))

    ### (65e)
    @constraint(model, eta_qliq_riser_m[q_liq_bp in riser.Q_liq], ξ_qliq_m[q_liq_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, eta_gor_riser_m[gor_bp in riser.GOR], ξ_gor_m[gor_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, eta_wct_riser_m[wct_bp in riser.WCT], ξ_wct_m[wct_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, eta_iglr_riser_m[iglr_bp in riser.IGLR], ξ_iglr_m[iglr_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
    ))
    @constraint(model, sos2_qliq_riser_m, ξ_qliq_m[riser.Q_liq] in SOS2(riser.Q_liq))
    @constraint(model, sos2_gor_riser_m, ξ_gor_m[riser.GOR] in SOS2(riser.GOR))
    @constraint(model, sos2_wct_riser_m, ξ_wct_m[riser.WCT] in SOS2(riser.WCT))
    @constraint(model, sos2_iglr_riser_m, ξ_iglr_m[riser.IGLR] in SOS2(riser.IGLR))

    ## Multilinear Interpolation

    @variable(model, bilinear_left[riser.Q_liq, riser.GOR])
    @variable(model, bilinear_right[riser.WCT, riser.IGLR])

    # @NonlinearConstraint
    @constraint(model, riser_bilinear_left[q_liq_bp in riser.Q_liq, gor_bp in riser.GOR],
        bilinear_left[q_liq_bp, gor_bp] == ξ_qliq_m[q_liq_bp] * ξ_gor_m[gor_bp]
    )
    # @NonlinearConstraint
    @constraint(model, riser_bilinear_right[wct_bp in riser.WCT, iglr_bp in riser.IGLR],
        bilinear_right[wct_bp, iglr_bp] == ξ_wct_m[wct_bp] * ξ_iglr_m[iglr_bp]
    )
    # @NonlinearConstraint
    @constraint(model, riser_bilinear[q_liq_bp in riser.Q_liq, gor_bp in riser.GOR, wct_bp in riser.WCT, iglr_bp in riser.IGLR],
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp] == bilinear_left[q_liq_bp,gor_bp] * bilinear_right[wct_bp,iglr_bp]
    )

    ## Infeasible points
    for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
            for wct_bp in riser.WCT
                for iglr_bp in riser.IGLR
                    if riser.ΔP[(q_liq_bp, gor_bp, wct_bp, iglr_bp)] < 0
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

    if ~isnothing(platform.riser)
        N_man = names(platform.riser.manifold_wells)

        N = [N_sat ; N_man]
    else
        N = N_sat
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
    if ~isnothing(platform.riser)
        ### (66d)
        @constraint(
            model, manifold_wells_downstream_p[n in N_man],
            variable_by_name(model, "dsp_$n") == variable_by_name(model, "p_sup_m")
                                                 + variable_by_name(model, "ΔP")
        )

        ### (66e)
        @constraint(
            model, downstream_pressure_manifold_m,
            variable_by_name(model, "dsp_m") == platform.p_sep
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

function get_minlp_problem(platform::Platform)
    model = Model(SCIP.Optimizer)

    # Objective
    q_oil_total = 0

    # Wells
    for well in platform.satellite_wells
        model = add_nonlinear_well(model, well)
        q_oil_total = q_oil_total + variable_by_name(model, "q_oil_$(well.name)")
    end

    # Manifold
    if ~isnothing(platform.riser)
        model = add_nonlinear_riser(model, platform.riser, platform.p_sep)  # (65)
        q_oil_total = q_oil_total + variable_by_name(model, "q_oil_m")
    end

    # Platform
    model = add_platform(model, platform)

    # Objective
    @objective(model, Max, q_oil_total)  # (64a)/(66a)

    return model
end

function add_linear_well(model::GenericModel, well::Oil.AbstractWell)
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
    # @McCormickEnvelope
    q_liq_L = min(well.Q_liq_vlp...)
    q_liq_U = max(well.Q_liq_vlp...)
    # TODO: these can probably be tightened with well parameters and other constraints
    iglr_L = min(well.IGLR...)
    iglr_U = max(well.IGLR...)
    @constraint(model, q_inj_n >= q_liq_L * iglr_n + q_liq_n * iglr_L - q_liq_L * iglr_L, base_name = "q_inj_n_nccormick_under_L")
    @constraint(model, q_inj_n >= q_liq_U * iglr_n + q_liq_n * iglr_U - q_liq_U * iglr_U, base_name = "q_inj_n_nccormick_under_U")
    @constraint(model, q_inj_n <= q_liq_U * iglr_n + q_liq_n * iglr_L - q_liq_U * iglr_L, base_name = "q_inj_n_nccormick_over_UL")
    @constraint(model, q_inj_n <= q_liq_L * iglr_n + q_liq_n * iglr_U - q_liq_L * iglr_U, base_name = "q_inj_n_nccormick_over_LU")

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

    ### (62m)
    @constraint(model, ξ_iglr_n[well.IGLR] in SOS2(well.IGLR), base_name = "sos2_iglr_$n")

    ### (62n)
    @constraint(model, ξ_whp_n[well.WHP] in SOS2(well.WHP), base_name = "sos2_whp_$n")

    ### (62o)
    @constraint(model, ξ_qliq_n[well.Q_liq_vlp] in SOS2(well.Q_liq_vlp), base_name = "sos2_qliq_vlp_$n")

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

function add_linear_riser(model::GenericModel, riser::Riser, p_sep::Float64)
    N_man = names(riser.manifold_wells)  # names of all wells in manifold

    # Bounds for McCormick Envelope
    q_liq_lowers = Vector{Float64}()
    q_liq_uppers = Vector{Float64}()
    iglr_lowers = Vector{Float64}()
    iglr_uppers = Vector{Float64}()
    wcts = [well.wct for well in riser.manifold_wells]
    gors = [well.gor for well in riser.manifold_wells]

    for well in riser.manifold_wells
        model = add_nonlinear_well(model, well)

        pwl_well = Oil.PiecewiseLinearWell(well, true)
        push!(q_liq_lowers, min(pwl_well.Q_liq_vlp...))
        push!(q_liq_uppers, max(pwl_well.Q_liq_vlp...))
        # TODO: I believe these can be tightened through min/max_q_inj
        push!(iglr_lowers, min(pwl_well.IGLR...))
        push!(iglr_uppers, max(pwl_well.IGLR...))
    end

    q_liq_L = sum(q_liq_lowers...)
    q_liq_U = sum(q_liq_uppers...)
    iglr_L = min(iglr_lowers...)
    iglr_U = max(iglr_uppers...)
    wct_L = min(wcts...)
    wct_U = max(wcts...)
    gor_L = min(gors...)
    gor_U = max(gors...)

    riser = Oil.PiecewiseLinearRiser(riser, p_sep, false)

    # Variables

    ## Flow rates
    @variables(model, begin
        q_liq_m >= 0.0
        q_oil_m >= 0.0
        q_water_m >= 0.0
        q_gas_m >= 0.0
        q_inj_m >= 0.0
    end)

    ## Pressures
    @variables(model, begin
        dsp_m   >= 0.0  # manifold downstream pressure
        p_sup_m >= 0.0  # manifold choke pressure
        ΔP      >= 0.0  # riser pressure drop
    end)

    ## Parameters
    @variables(model, begin
        gor_m
        wct_m
        iglr_m
    end)

    ## Discrete decisions
    @variable(model, y_m, Bin)  # (65d)

    ## VLP curve variables
    λ_riser_m = @variable(model, [riser.Q_liq, riser.GOR, riser.WCT, riser.IGLR],
                          lower_bound=0.0, base_name="λ_riser_m")  # (65d)
    ξ_qliq_m = @variable(model, [riser.Q_liq], lower_bound=0.0, base_name="ξ_qliq_m")
    ξ_gor_m = @variable(model, [riser.GOR], lower_bound=0.0, base_name="ξ_gor_m")
    ξ_wct_m = @variable(model, [riser.WCT], lower_bound=0.0, base_name="ξ_wct_m")
    ξ_iglr_m = @variable(model, [riser.IGLR], lower_bound=0.0, base_name="ξ_iglr_m")

    # Constraints

    ## Manifold input

    ### (65a)
    @constraint(model, manifold_liquid_rate, q_liq_m == sum(variable_by_name(model, "q_liq_$n") for n in N_man))
    @constraint(model, manifold_oil_rate, q_oil_m == sum(variable_by_name(model, "q_oil_$n") for n in N_man))
    @constraint(model, manifold_water_rate, q_water_m == sum(variable_by_name(model, "q_water_$n") for n in N_man))
    @constraint(model, manifold_gas_rate, q_gas_m == sum(variable_by_name(model, "q_gas_$n") for n in N_man))
    @constraint(model, manifold_inj_gas_rate, q_inj_m == sum(variable_by_name(model, "q_inj_$n") for n in N_man))

    ## VLP

    ### (65b)
    @constraint(model, pwl_qliq_riser_m, q_liq_m == sum(
        q_liq_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, pwl_gor_riser_m, gor_m == sum(
        gor_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, pwl_wct_riser_m, wct_m == sum(
        wct_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, pwl_iglr_riser_m, iglr_m == sum(
        iglr_bp * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, p_sup_upper, p_sup_m <= dsp_m + M_Pressure * (1 - y_m))
    @constraint(model, riser_choke_dP, p_sup_m >= dsp_m - M_Pressure * (1 - y_m))

    ### (65c)
    @constraint(model, q_oil_m == q_liq_m - q_water_m, base_name = "bilinear_q_oil_riser_m")
    # @McCormickEnvelope
    @constraint(model, q_water_m >= q_liq_L * wct_m + q_liq_m * wct_L - q_liq_L * wct_L, base_name = "q_water_m_mccormick_under_L")
    @constraint(model, q_water_m >= q_liq_U * wct_m + q_liq_m * wct_U - q_liq_U * wct_U, base_name = "q_water_m_mccormick_under_U")
    @constraint(model, q_water_m <= q_liq_U * wct_m + q_liq_m * wct_L - q_liq_U * wct_L, base_name = "q_water_m_mccormick_over_UL")
    @constraint(model, q_water_m <= q_liq_L * wct_m + q_liq_m * wct_U - q_liq_L * wct_U, base_name = "q_water_m_mccormick_over_LU")
    # @McCormickEnvelope
    @constraint(model, q_gas_m >= q_oil_L * gor_m + q_oil_m * gor_L - q_oil_L * gor_L, base_name = "q_gas_m_mccormick_under_L")
    @constraint(model, q_gas_m >= q_oil_U * gor_m + q_oil_m * gor_U - q_oil_U * gor_U, base_name = "q_gas_m_mccormick_under_U")
    @constraint(model, q_gas_m <= q_oil_U * gor_m + q_oil_m * gor_L - q_oil_U * gor_L, base_name = "q_gas_m_mccormick_over_UL")
    @constraint(model, q_gas_m <= q_oil_L * gor_m + q_oil_m * gor_U - q_oil_L * gor_U, base_name = "q_gas_m_mccormick_over_LU")
    # @McCormickEnvelope
    @constraint(model, q_inj_m >= q_liq_L * iglr_m + q_liq_m * iglr_L - q_liq_L * iglr_L, base_name = "q_inj_m_mccormick_under_L")
    @constraint(model, q_inj_m >= q_liq_U * iglr_m + q_liq_m * iglr_U - q_liq_U * iglr_U, base_name = "q_inj_m_mccormick_under_U")
    @constraint(model, q_inj_m <= q_liq_U * iglr_m + q_liq_m * iglr_L - q_liq_U * iglr_L, base_name = "q_inj_m_mccormick_over_UL")
    @constraint(model, q_inj_m <= q_liq_L * iglr_m + q_liq_m * iglr_U - q_liq_L * iglr_U, base_name = "q_inj_m_mccormick_over_LU")

    @constraint(model, pwl_ΔP_riser_m, ΔP == sum(
        riser.ΔP[q_liq_bp, gor_bp, wct_bp, iglr_bp] * λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))

    ### (65d)
    @constraint(model, pwl_riser_curve, y_m == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, forced_activation_riser_m[n in N_man], y_m >= variable_by_name(model, "y_$n"))
    @constraint(model, forced_deactivation_riser_m, y_m <= sum(variable_by_name(model, "y_$n") for n in N_man))

    ### (65e)
    @constraint(model, eta_qliq_riser_m[q_liq_bp in riser.Q_liq], ξ_qliq_m[q_liq_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, eta_gor_riser_m[gor_bp in riser.GOR], ξ_gor_m[gor_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for wct_bp in riser.WCT
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, eta_wct_riser_m[wct_bp in riser.WCT], ξ_wct_m[wct_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for iglr_bp in riser.IGLR
    ))
    @constraint(model, eta_iglr_riser_m[iglr_bp in riser.IGLR], ξ_iglr_m[iglr_bp] == sum(
        λ_riser_m[q_liq_bp, gor_bp, wct_bp, iglr_bp]
        for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
        for wct_bp in riser.WCT
    ))
    @constraint(model, sos2_qliq_riser_m, ξ_qliq_m[riser.Q_liq] in SOS2(riser.Q_liq))
    @constraint(model, sos2_gor_riser_m, ξ_gor_m[riser.GOR] in SOS2(riser.GOR))
    @constraint(model, sos2_wct_riser_m, ξ_wct_m[riser.WCT] in SOS2(riser.WCT))
    @constraint(model, sos2_iglr_riser_m, ξ_iglr_m[riser.IGLR] in SOS2(riser.IGLR))

    ## Infeasible points
    for q_liq_bp in riser.Q_liq
        for gor_bp in riser.GOR
            for wct_bp in riser.WCT
                for iglr_bp in riser.IGLR
                    if riser.ΔP[(q_liq_bp, gor_bp, wct_bp, iglr_bp)] < 0
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

function get_milp_relaxation(platform::Platform)
    model = Model(SCIP.Optimizer)

    # Objective
    q_oil_total = 0

    # Wells
    for well in platform.satellite_wells
        model = add_nonlinear_well(model, well)
        q_oil_total = q_oil_total + variable_by_name(model, "q_oil_$(well.name)")
    end

    # Manifold
    if ~isnothing(platform.riser)
        model = add_nonlinear_riser(model, platform.riser, platform.p_sep)  # (65)
        q_oil_total = q_oil_total + variable_by_name(model, "q_oil_m")
    end

    # Platform
    model = add_platform(model, platform)

    # Objective
    @objective(model, Max, q_oil_total)  # (64a)/(66a)

    return model
end

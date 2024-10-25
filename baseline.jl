using Oil

include("models.jl")

kgf, g, m3, d, kPa = latin_si(:kgf), latin_si(:gauge), latin_si(:m3), latin_si(:day), latin_si(:kPa)

platform = Platform(
    10.001 * kgf + g,
    satellite_wells = [
        Well("P05", 65.0, 0.55, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(175.0 * kgf + g, 44.123932 * (m3 / d) / kgf)),
        # Well("P13", 70.0, 0.25, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(220.0 * kgf + g, 117.64259 * (m3 / d) / kgf)),
    ],
    # riser = Riser(
    #     VLP("data/MSP_UEP_VFP.Ecl"),
    #     [
    #         Well("P01", 70.0, 0.10, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(200*kgf + g, 98.11 * (m3 / d) / kgf)),
    #         # Well("P02", 120.0, 0.50, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(180.0 * kgf + g, 48.33 * (m3 / d) / kgf)),
    #         # Well("P03", 100.0, 0.05, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(170.0 * kgf + g, 28.08 * (m3 / d) / kgf)),
    #         # Well("P04", 150.0, 0.75, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(220.0 * kgf + g, 73.01 * (m3 / d) / kgf)),
    #     ],
    #     choke_enabled = false
    # ),
    q_inj_max=10*1e3
)

model = get_minlp_problem(platform)
# model = get_milp_relaxation(platform)

optimize!(model)

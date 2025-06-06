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

## SCENARIO 2
platform = Platform(
    10.001 * kgf + g,
    # satellite_wells = Vector{Well}(),
    satellite_wells = [
        Well("P05", 65.0, 0.55, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(175.0 * kgf + g, 44.123932 * (m3 / d) / kgf)),
        Well("P13", 70.0, 0.25, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_UEP_VLP.Ecl"), IPR(220.0 * kgf + g, 117.64259 * (m3 / d) / kgf)),
    ],
    manifolds = [
        Manifold(
            VLP("data/MSP_UEP_VFP.Ecl"),
            [
                Well("P01", 70.0, 0.10, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(200*kgf + g, 98.11 * (m3 / d) / kgf)),
                Well("P02", 120.0, 0.50, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(180.0 * kgf + g, 48.33 * (m3 / d) / kgf)),
                Well("P03", 100.0, 0.05, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(170.0 * kgf + g, 28.08 * (m3 / d) / kgf)),
                Well("P04", 150.0, 0.75, 100 * 1e3*m3/d, 200 * 1e3*m3/d, VLP("data/Well_SubseaManifold_VLP.Ecl"), IPR(220.0 * kgf + g, 73.01 * (m3 / d) / kgf)),
            ],
            choke_enabled = false
        ),
    ],
)

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
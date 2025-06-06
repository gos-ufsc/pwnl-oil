using Random
using JSON3
using Oil

include("models.jl")

function generate_random_well_params(name::String;
    gor_range=(45.0, 500.0),
    wct_range=(0.1, 0.9),
    p_res_range=(180.0, 270.0),
    IP_range=(30.0, 150.0),
    vlp_type="UEP"
)
    params = Dict(
        :name => name,
        :gor => rand() * (gor_range[2] - gor_range[1]) + gor_range[1],
        :wct => rand() * (wct_range[2] - wct_range[1]) + wct_range[1],
        :p_res => rand() * (p_res_range[2] - p_res_range[1]) + p_res_range[1],
        :IP => rand() * (IP_range[2] - IP_range[1]) + IP_range[1],
        :q_min => 100.0e3,
        :q_max => 200.0e3,
        :vlp_type => vlp_type
    )
    # Ensure q_max > q_min
    params[:q_max] = max(params[:q_max], params[:q_min] + 1e3)
    return params
end

function create_scenario_instance(scenario_id::Int, instance_id::Int;
    n_satellite=4,
    n_manifolds=2,
    wells_per_manifold=3,
    base_path="scenarios_teste",
    param_ranges=Dict()
)
    scenario_dir = "$base_path/scenario_$scenario_id"
    isdir(scenario_dir) || mkpath(scenario_dir)
    total_wells = n_satellite + (n_manifolds * wells_per_manifold)

    platform_constraints = Dict(
        :q_liq_max_plat => (rand() * (1200.0 - 700.0) + 700.0) * total_wells,
        :q_inj_max_plat => (rand() * (150.0e3 - 100.0e3) + 100.0e3) * total_wells
    )
    
    # Generate compact data structure
    instance_data = Dict(
        :metadata => merge(Dict(
            :scenario_id => scenario_id,
            :instance_id => instance_id,
            :p_sep => 10.001
        ),platform_constraints),
        :satellite_wells => [
            generate_random_well_params("SAT_$(instance_id)_$i"; 
                param_ranges..., vlp_type="UEP")
            for i in 1:n_satellite
        ],
        :manifolds => [
            Dict(
                :wells => [
                    generate_random_well_params("MAN_$(instance_id)_$(m)_$i"; 
                        param_ranges..., vlp_type="Subsea")
                    for i in 1:wells_per_manifold
                ],
                :choke_enabled => false
            ) for m in 1:n_manifolds
        ]
    )

    # Save pretty-printed JSON
    output_path = "$scenario_dir/instance_$instance_id.json"
    open(output_path, "w") do f
        JSON3.pretty(f, instance_data)
    end
end

# Reconstruction function (for later use)
function load_instance(scenario_id, instance_id; base_path="scenarios")
    data = JSON3.read("$base_path/scenario_$scenario_id/instance_$instance_id.json")
    
    # Rebuild objects from parameters
    satellite_wells = [
        Well(w.name, w.gor, w.wct, w.q_min, w.q_max, 
            VLP("data/Well_$(w.vlp_type)_VLP.Ecl"), 
            IPR(w.p_res, w.IP)) for w in data.satellite_wells
    ]
    
    manifolds = [
        Manifold(
            VLP("data/MSP_UEP_VFP.Ecl"),
            [Well(w.name, w.gor, w.wct, w.q_min, w.q_max, 
                VLP("data/Well_SubseaManifold_VLP.Ecl"), 
                IPR(w.p_res, w.IP)) for w in m.wells],
            m.choke_enabled
        ) for m in data.manifolds
    ]
    
    Platform(
        data.metadata.p_sep,
        satellite_wells,
        manifolds,
        nothing, nothing, nothing, nothing
    )
end

function generate_all_scenarios(;
    scenarios=[
        (id=1, n_sat=1, n_man=0, wells_per_man=0),
        (id=2, n_sat=2, n_man=0, wells_per_man=0),
        (id=3, n_sat=3, n_man=0, wells_per_man=0),
        (id=4, n_sat=2, n_man=1, wells_per_man=2),
        (id=5, n_sat=2, n_man=1, wells_per_man=3),
        (id=6, n_sat=3, n_man=1, wells_per_man=3),
        (id=7, n_sat=3, n_man=2, wells_per_man=2),
        (id=8, n_sat=2, n_man=2, wells_per_man=3),
        (id=9, n_sat=3, n_man=2, wells_per_man=3),
    ],
    instances_per_scenario=8,
    param_ranges=Dict(
        :gor_range => (45.0, 200.0),
        :wct_range => (0.1, 0.9),
        :p_res_range => (180.0, 270.0),
        :IP_range => (30.0, 150.0),
    ),
    base_path = "scenarios"
)
    for (i, scenario) in enumerate(scenarios)
        for instance_id in 1:instances_per_scenario
            create_scenario_instance(
                scenario.id,
                instance_id;
                n_satellite=scenario.n_sat,
                n_manifolds=scenario.n_man,
                wells_per_manifold=scenario.wells_per_man,
                base_path=base_path,
                param_ranges=param_ranges
            )
            println("Created scenario $i out of $(length(scenarios)) (instance $instance_id)")
        end
    end
end

# Generate all scenarios and instances
generate_all_scenarios()
import glob
import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'vibrant'])

fig, ax = plt.subplots(figsize=(6, 4))

INF = float('inf')

for result_dir in glob.glob("results_*"):
    method = os.path.basename(result_dir)[len("results_"):]

    print(f"Processing results for method: {method}")

    feas_times = [0.0]
    opt_times = [0.0]
    for i in range(1, 7):  # scenario id 1 to 6
        print(f"Scenario {i}:")
        optimal = 0
        feasible = 0
        min_lb = INF
        avg_gap = 0.0
        max_gap = 0.0
        for j in range(1, 9):  # instance id 1 to 8
            path = os.path.join(result_dir, f"scenario_{i}", f"instance_{j}.json")
            if os.path.isfile(path):
                with open(path, "r") as f:
                    results = json.load(f)

                # Process results as needed
                if results["objective_value"] > 0.0:
                    feasible += 1
                    # Find first time where lower > 1e-3
                    lowers = results["lowers"]
                    times = results["times"]
                    first_feas_t = next((t for l, t in zip(lowers, times) if l > 1e-3), None)
                    if first_feas_t is not None:
                        feas_times.append(first_feas_t)
                if results["status"] == "OPTIMAL":
                    optimal += 1
                    opt_times.append(results["solve_time"])

                gap = results["gap"]
                avg_gap += gap / 8
                max_gap = max(max_gap, gap)
                min_lb = min(min_lb, results["objective_value"])

            else:
                print(f"ERROR: No results found for scenario {i}, instance {j} at path: {path}")

        print(f"\t Optimal = \t {optimal} / 8")
        print(f"\t Feasible = \t {feasible} / 8")
        print(f"\t Average Gap = \t {100*avg_gap} %")
        print(f"\t Max Gap = \t {100*max_gap} %")
        print(f"\t Min LB = \t {min_lb}")
    
    feas_times = np.sort(feas_times)
    opt_times = np.sort(opt_times)

    # For step plots, we need the "n" and "t" arrays as in Julia
    opt_times_n = np.arange(len(opt_times) + 1)
    opt_times_t = np.append(opt_times, 3600.0)
    feas_times_n = np.arange(len(feas_times) + 1)
    feas_times_t = np.append(feas_times, 3600.0)

    l = ax.step(feas_times_t, feas_times_n, where='post', linestyle='--')
    ax.step(opt_times_t, opt_times_n, where='post', linestyle='-', label=method, color=l[0].get_color())

plt.xlim(0, 3600)
plt.ylim(0, 6*8)
plt.xlabel("Time [s]")
plt.ylabel("Number of solved instances")
plt.legend()
fig.savefig("results.pdf", bbox_inches='tight')
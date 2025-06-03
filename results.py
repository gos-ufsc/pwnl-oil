import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

from glob import glob

plt.style.use(['science', 'vibrant'])

def add_step_curve(x1,y1,x2,y2):
    """
    Add two step curves.
    """
    # assert y1[0] == 0.0
    # assert y2[0] == 0.0

    x = np.union1d(x1, x2)
    x = np.sort(x)
    assert x[0] == 0.0

    delta1 = y1 - np.roll(y1, 1)
    delta2 = y2 - np.roll(y2, 1)
    delta1[0] = y1[0]
    delta2[0] = y2[0]

    y = np.zeros(x.shape)
    for i in np.arange(1,len(x)):
        y[i] = y[i-1]
        if x[i] in x1:
            y[i] += delta1[np.where(x1 == x[i])[0][0]]
        if x[i] in x2:
            y[i] += delta2[np.where(x2 == x[i])[0][0]]

    return x, y

fig_instances, ax_instances = plt.subplots(figsize=(6, 4))
fig_gaps, ax_gaps = plt.subplots(figsize=(6, 4))
fig_gaps_scatter, ax_gaps_scatter = plt.subplots(figsize=(6, 4))
fig_lbs, ax_lbs = plt.subplots(figsize=(6, 4))

n_methods = len(list(glob("results_*")))
width = 0.4 / n_methods

INF = float('inf')

for k, result_dir in enumerate(glob("results_*")):
    method = os.path.basename(result_dir)[len("results_"):]

    print(f"Processing results for method: {method}")

    t = [0.0]
    lb = [0.0]
    feas_times = [0.0]
    opt_times = [0.0]
    max_gaps = []
    avg_gaps = []
    min_gaps = []
    for i in range(1, 7):  # scenario id 1 to 6
        print(f"Scenario {i}:")
        optimal = 0
        feasible = 0
        min_lb = INF
        avg_gap = 0.0
        max_gap = 0.0
        min_gap = INF
        for j in range(1, 9):  # instance id 1 to 8
            path = os.path.join(result_dir, f"scenario_{i}", f"instance_{j}.json")
            if os.path.isfile(path):
                with open(path, "r") as f:
                    results = json.load(f)

                lb_ref = 0.0
                if method.startswith("gurobi"):
                    # Process results as needed
                    if results["objective_value"] > 0.0:
                        feasible += 1
                        # Find first time where lower > 1e-3
                        lowers = np.array(results["lowers"])
                        times = np.array(results["times"])
                        first_feas_t = next((t for l, t in zip(lowers, times) if l > 1e-3), None)
                        if first_feas_t is not None:
                            feas_times.append(first_feas_t)
                    if results["status"] == "OPTIMAL":
                        optimal += 1
                        opt_times.append(results["solve_time"])

                        lb_ref = results["lower_bound"]
                    
                    # ax_lbs.step(
                    #     times, np.array(lowers) / lb_ref,
                    #     alpha=0.25,
                    #     color=f"C{k}",
                    #     where="post",
                    # )

                    gap = results["gap"]
                    max_gap = max(max_gap, gap)
                    avg_gap += gap / 8
                    min_gap = min(min_gap, gap)
                else:
                    # Process results as needed
                    if results["objective_value"] > 0.0:
                        feasible += 1
                        # Find first time where lower > 1e-3
                        for iteration in results["iterations"]:
                            if iteration["upper_bound"] is not None:
                                if -iteration["upper_bound"] > 1e-3:  # UB in minimization
                                    feas_times.append(iteration["time"])
                                    break

                        times = [0.0]
                        lowers = [0.0]
                        for iteration in results["iterations"]:
                            if iteration["upper_bound"] is not None:
                                times.append(iteration["time"])
                                lowers.append(-iteration["upper_bound"])
                        times = np.array(times)
                        lowers = np.array(lowers)

                    if results["status"] == "CONVERGED":
                        optimal += 1
                        opt_times.append(results["solve_time"])

                        lb_ref = results["lower_bound"]

                    # ax_lbs.step(
                    #     times, np.array(lowers) / lb_ref,
                    #     alpha=0.25,
                    #     color=f"C{k}",
                    #     where="post",
                    # )

                    gap = results["gap"]
                    max_gap = max(max_gap, gap)
                    avg_gap += gap / 8
                    min_gap = min(min_gap, gap)

                if lb_ref == 0.0:
                    for result_path in glob(f"results_*/scenario_{i}/instance_{j}.json"):
                        with open(result_path, "r") as f:
                            result = json.load(f)
                            if result["lower_bound"] > lb_ref:
                                lb_ref = result["lower_bound"]
                lowers[lowers < 0] = 0.0
                lowers = lowers / lb_ref

                ax_gaps_scatter.scatter(
                    i + k * width + width / 2, 100 * gap,
                    marker='_',
                    s=50,
                    label=method if i == 1 and j == 1 else None,
                    color=f"C{k}",
                )

                t, lb = add_step_curve(t, lb, times, lowers / (6 * 8))
            else:
                print(f"ERROR: No results found for scenario {i}, instance {j} at path: {path}")

        print(f"\t Optimal = \t {optimal} / 8")
        print(f"\t Feasible = \t {feasible} / 8")
        print(f"\t Max Gap = \t {100*max_gap} %")
        print(f"\t Average Gap = \t {100*avg_gap} %")
        print(f"\t Min Gap = \t {100*min_gap} %")

        max_gaps.append(max_gap)
        avg_gaps.append(avg_gap)
        min_gaps.append(min_gap)

    feas_times = np.sort(feas_times)
    opt_times = np.sort(opt_times)

    # For step plots, we need the "n" and "t" arrays as in Julia
    opt_times_n = np.arange(len(opt_times) + 1)
    opt_times_t = np.append(opt_times, 3600.0)
    feas_times_n = np.arange(len(feas_times) + 1)
    feas_times_t = np.append(feas_times, 3600.0)

    # plot instances solved over time
    l = ax_instances.step(feas_times_t, feas_times_n, where='post', linestyle='--')
    ax_instances.step(opt_times_t, opt_times_n, where='post', linestyle='-', label=method, color=l[0].get_color())

    # plot gaps
    max_gaps = np.array(max_gaps)
    avg_gaps = np.array(avg_gaps)
    min_gaps = np.array(min_gaps)
    ax_gaps.errorbar(
        np.arange(1,7) + k * width + width / 2, 100 * avg_gaps,
        yerr=100 * np.vstack([avg_gaps - min_gaps, max_gaps - avg_gaps]),
        linewidth=0,
        elinewidth=1, marker='_', markersize=8,
        # fmt='o', capsize=5, label=method if i == 1 else None,
        label=method, color=l[0].get_color(),
        # ecolor=f"C{k}", elinewidth=2, markersize=8
    )

    # eb = ax_gaps_scatter.errorbar(
    #     np.arange(1,7) + k * width + width / 2, 100 * avg_gaps,
    #     yerr=100 * np.vstack([avg_gaps - min_gaps, max_gaps - avg_gaps]),
    #     linewidth=0,
    #     elinewidth=1, marker='^', markersize=3,
    #     # fmt='o', capsize=5, label=method if i == 1 else None,
    #     label=method, color=l[0].get_color(),
    #     # ecolor=f"C{k}", elinewidth=2, markersize=8
    # )
    # eb[-1][0].set_linestyle('--')
    ax_gaps_scatter.scatter(
        np.arange(1,7) + k * width + width / 2, 100 * avg_gaps,
        marker='*',
        s=25,
        color=f"C{k}",
    )


    ax_lbs.step(
        t, np.array(lb),
        color=f"C{k}",
        where="post",
        label=method
    )
    

ax_instances.set_xlim(0, 3600)
ax_instances.set_ylim(0, 6*8)
ax_instances.set_xlabel("Time [s]")
ax_instances.set_ylabel("Number of solved instances")
ax_instances.legend()
fig_instances.savefig("results.pdf", bbox_inches='tight')

# ax_gaps.set_xticks(np.arange(1, 7) + width * (len(list(glob.glob("results_*"))) - 1) / 2)
ax_gaps.set_xticks(np.arange(1, 7) + width * n_methods / 2)
ax_gaps.set_xticklabels([f"S{i}" for i in range(1, 7)])
ax_gaps.set_xlim(3.5, 7)
ax_gaps.set_ylim(1e-4, 1000.0)
ax_gaps.set_yscale("log")
ax_gaps.set_ylabel("Gap [\%]")
ax_gaps.legend(loc='upper left')
fig_gaps.savefig("gaps.pdf", bbox_inches='tight')

ax_gaps_scatter.set_xticks(np.arange(1, 7) + width * n_methods / 2)
ax_gaps_scatter.set_xticklabels([f"S{i}" for i in range(1, 7)])
# ax_gaps_scatter.set_xlim(3.5, 7)
# ax_gaps_scatter.set_ylim(1e-4, 1000.0)
ax_gaps_scatter.set_yscale("log")
ax_gaps_scatter.set_ylabel("Gap [\%]")
ax_gaps_scatter.legend(loc='upper left')
fig_gaps_scatter.savefig("gaps_scatter.pdf", bbox_inches='tight')

ax_lbs.set_ylim(0, 1)
ax_lbs.set_xlim(0, 3600)
ax_lbs.set_xlabel("Time [s]")
ax_lbs.set_ylabel("Average Normalized Lower Bound")
ax_lbs.legend()
fig_lbs.savefig("lbs.pdf", bbox_inches='tight')
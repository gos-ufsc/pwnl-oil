from copy import copy
import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

from glob import glob

plt.style.use(['science', 'vibrant'])

np.random.seed(33)

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

eg_i, eg_j = np.random.randint(1, 7), np.random.randint(1, 9)

n_methods = len(list(glob("results_*")))
width = 0.4 / n_methods

INF = float('inf')
eps_opt = 0.0001

# methods = ["rfe", "gurobi_sos1", "gurobi_binary", "scip", "scip_bin", "baron"]
methods = ["rfe", "gurobi_sos2", "gurobi_binary", "scip", "baron"]

labels = {
    "rfe": "RFE",
    "gurobi_sos2": "Gurobi (SOS2)",
    "gurobi_binary": "Gurobi (Bin.)",
    "baron": "BARON",
    "scip": "SCIP",
    # "scip": "SCIP (SOS2)",
    # "scip_bin": "SCIP (Bin.)",
}

for k, method in enumerate(methods):
    print(f"Processing results for method: {method}")

    result_dir = f"results_{method}"

    t = [0.0]
    lb = [0.0]
    feas_times = [0.0]
    opt_times = [0.0]
    max_gaps = []
    avg_gaps = []
    min_gaps = []
    for i in range(1, 10):  # scenario id 1 to 6
        print(f"Scenario {i}:")
        optimal = 0
        feasible = 0
        min_lb = INF
        avg_gap = 0.0
        max_gap = 0.0
        min_gap = INF
        max_time = 0.0
        avg_time = 0.0

        if len(glob(f"{result_dir}/scenario_{i}/*.json")) < 8:
            print(f"\t Skipping (incomplete results).")
            continue  # Skip scenarios with less than 8 instances

        for j in range(1, 9):  # instance id 1 to 8
            path = os.path.join(result_dir, f"scenario_{i}", f"instance_{j}.json")
            if os.path.isfile(path):
                with open(path, "r") as f:
                    results = json.load(f)

                solve_time = results.get("solve_time", None)
                if solve_time is None:
                    solve_time = 3600.0  # Default to 3600 seconds if not present
                solve_time = min(solve_time, 3600.0)  # Cap at 3600 seconds
                avg_time += solve_time / 8
                max_time = max(max_time, solve_time)

                gap = INF
                lb_ref = 0.0
                if method.startswith("rfe"):
                    # Process results as needed
                    if (results["objective_value"] is not None) and (results["objective_value"] > 0.0):
                        feasible += 1
                        # Find first time where lower > 1e-3
                        for iteration in results["iterations"]:
                            if iteration["upper_bound"] is not None:
                                if -iteration["upper_bound"] > 1e-3:  # UB in minimization
                                    feas_times.append(iteration["time"])
                                    break

                        times = [0.0]
                        lowers = [0.0]
                        uppers = [INF]
                        for iteration in results["iterations"]:
                            times.append(iteration["time"])
                            if iteration["upper_bound"] is not None:
                                lowers.append(-iteration["upper_bound"])
                            else:
                                lowers.append(lowers[-1])

                            if iteration["lower_bound"] is not None:
                                uppers.append(-iteration["lower_bound"])
                            else:
                                uppers.append(uppers[-1])
                        times = np.array(times)
                        lowers = np.array(lowers)
                        uppers = np.array(uppers)

                        g_times = times[(uppers < 100000) & (lowers > 0)]
                        g_uppers = uppers[(uppers < 100000) & (lowers > 0)]
                        g_lowers = lowers[(uppers < 100000) & (lowers > 0)]
                        gaps = (g_uppers - g_lowers) / g_lowers

                        u_times = times[uppers < 100000]
                        uppers = uppers[uppers < 100000]

                        times = times[lowers > 0]
                        lowers = lowers[lowers > 0]

                        if np.any(gaps < eps_opt):
                            gap = min(gaps[-1], results["gap"])
                            gap = max(0.0, gap)
                            optimal += 1
                            opt_times.append(g_times[gaps < eps_opt][0])

                            lb_ref = results["objective_value"]
                        elif results["status"] == "CONVERGED":
                            gap = max(0.0, results["gap"])
                            optimal += 1
                            opt_times.append(results["solve_time"])

                            lb_ref = results["objective_value"]

                    if (gap == INF) and (results["gap"] is not None):
                        gap = results["gap"]

                    # ax_lbs.step(
                    #     times, np.array(lowers) / lb_ref,
                    #     alpha=0.25,
                    #     color=f"C{k}",
                    #     where="post",
                    # )
                else:
                    # Process results as needed
                    if (results["objective_value"] is not None) and (results["objective_value"] >= 0.1):
                        feasible += 1
                        # Find first time where lower > 1e-3
                        if method.startswith("scip"):
                            lowers = np.array(results["uppers"])
                            uppers = np.array(results["lowers"])
                        else:
                            lowers = np.array(results["lowers"])
                            uppers = np.array(results["uppers"])
                        times = np.array(results["times"])

                        g_times = times[(uppers < 100000) & (lowers > 0)]
                        g_uppers = uppers[(uppers < 100000) & (lowers > 0)]
                        g_lowers = lowers[(uppers < 100000) & (lowers > 0)]
                        gaps = (g_uppers - g_lowers) / g_lowers

                        u_times = times[uppers < 100000]
                        uppers = uppers[uppers < 100000]

                        times = np.append(times, results["solve_time"])
                        lowers = np.append(lowers, results["objective_value"])

                        times = times[lowers > 0]
                        lowers = lowers[lowers > 0]

                        feas_times.append(times[lowers > 1e-3][0])

                        if feas_times[-1] >= results["solve_time"]:
                            print(path)
                        # np.append(u_times, results["solve_time"])
                        # np.append(uppers, results["upper_bound"])

                        if ("gap" in results.keys()) and (results["gap"] is not None):
                            gap = results["gap"]

                        if np.any(gaps < eps_opt):
                            gap = min(gaps[-1], gap)
                            optimal += 1
                            opt_times.append(g_times[gaps < eps_opt][0])

                            lb_ref = results["objective_value"]
                        elif gap < eps_opt:
                            optimal += 1
                            opt_times.append(results["solve_time"])

                            # if results["solve_time"] >= 3600.0:
                            #     print(f"WARNING: Scenario {i}, Instance {j} took too long: {results['solve_time']} s. Path = {path}")

                            lb_ref = results["objective_value"]

                max_gap = max(max_gap, gap)
                avg_gap += gap / 8
                min_gap = min(min_gap, gap)

                if (method == "rfe") or (method.startswith("gurobi")):
                    lowers[lowers < 0] = 0.0

                    if lb_ref == 0.0:
                        for result_path in glob(f"results_*/scenario_{i}/instance_{j}.json"):
                            with open(result_path, "r") as f:
                                result = json.load(f)
                                if (result["objective_value"] is not None) and (result["objective_value"] > lb_ref):
                                    lb_ref = result["objective_value"]
                    lowers = lowers / lb_ref

                    ax_gaps_scatter.scatter(
                        i + k * width + width / 2, 100 * gap,
                        marker='_',
                        s=50,
                        label=labels[method] if i == 1 and j == 1 else None,
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
        print(f"\t Max Time = \t {max_time} s")
        print(f"\t Avg Time = \t {avg_time} s")

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
    ax_instances.step(opt_times_t, opt_times_n, where='post', linestyle='-', label=labels[method], color=f"C{k}")
    # l = ax_instances.step(feas_times_t, feas_times_n, where='post', linestyle='--')
    # ax_instances.step(opt_times_t, opt_times_n, where='post', linestyle='-', label=labels[method], color=l[0].get_color())

    # plot gaps
    if (method == "rfe") or (method.startswith("gurobi")):
        max_gaps = np.array(max_gaps)
        avg_gaps = np.array(avg_gaps)
        min_gaps = np.array(min_gaps)
        ax_gaps.errorbar(
            np.arange(1,len(avg_gaps)+1) + k * width + width / 2, 100 * avg_gaps,
            yerr=100 * np.vstack([avg_gaps - min_gaps, max_gaps - avg_gaps]),
            linewidth=0,
            elinewidth=1, marker='_', markersize=8,
            # fmt='o', capsize=5, label=method if i == 1 else None,
            label=labels[method], color=f"C{k}",
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
            np.arange(1,len(avg_gaps)+1) + k * width + width / 2, 100 * avg_gaps,
            marker='*',
            s=25,
            color=f"C{k}",
        )


        ax_lbs.step(
            t, np.array(lb),
            color=f"C{k}",
            where="post",
            label=labels[method]
        )


ax_instances.set_xlim(1., 3600)
ax_instances.set_xscale("log")
ax_instances.set_ylim(0, 9*8)
ax_instances.set_xlabel("Time [s]")
ax_instances.set_ylabel("Number of instances solved")
ax_instances.legend(loc='upper left')
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
ax_lbs.set_xlim(1, 3600)
ax_lbs.set_xscale("log")
ax_lbs.set_xlabel("Time [s]")
ax_lbs.set_ylabel("Average Normalized Lower Bound")
ax_lbs.legend()
fig_lbs.savefig("lbs.pdf", bbox_inches='tight')

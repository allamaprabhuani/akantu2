#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("TKAgg")
import matplotlib.pyplot as plt

# Same font as JOSS
plt.rcParams['font.sans-serif'] = 'cmss10'

print("Using:", matplotlib.get_backend())
plots = {
    "elastic_gcc": {"prefix": "timmings_", "suffix": "_jed"},
    "cohesive_gcc": {"prefix": "timmings_", "suffix": "_jed"},
}


axes = {}
fig, (axes["elastic"], axes["cohesive"]) = plt.subplots(2, 1, sharex=True)
handles = []
for plot_name, data in plots.items():
    mat, compiler = plot_name.split("_")
    data["df"] = pd.read_csv(
        f"""{data["prefix"]}{plot_name}{data["suffix"]}.csv""",
        sep=",",
        skipinitialspace=True,
    )

plotting = "TTS"

psize = np.array(sorted(set(plots[list(plots.keys())[0]]["df"]["psize"])))

for plot_name, data in plots.items():
    mat, compiler = plot_name.split("_")
    df = data["df"]

    measures = {}
    if mat == "elastic":
        step = df["solve_step"]
    elif mat == "cohesive":
        step = df["check_cohesive_stress"] + df["solve_step"]
    measures["TTS"] = step * df["solve_step nb_rep"]
    measures["speedup"] = step[0] / step
    measures["mumps"] = df["static_solve"]

    for name, measure in measures.items():
        data[name] = {
            "values": measure,
            "median": [],
            "mean": [],
            "std": [],
            "max": [],
            "min": [],
            "grouped": [],
            "nb_measure": [],
        }
        for s in psize:
            values = measure[df["psize"] == s]
            data[name]["std"].append(values.std())
            data[name]["mean"].append(values.mean())
            data[name]["median"].append(values.median())
            data[name]["max"].append(values.max())
            data[name]["min"].append(values.min())
            data[name]["grouped"].append(values)
            data[name]["nb_measure"].append(len(values))

    # axes[mat].boxplot(
    #     data[plotting]["grouped"], positions=psize, widths=[0.1 * s for s in psize]
    # )

    print(f"""{mat}: {data[plotting]["nb_measure"]}\n {psize}""")
    axes[mat].plot(
        psize,
        data[plotting]["median"],
        "-o",
        label=f"""Median {plotting}.""",
    )

    axes[mat].fill_between(
        psize,
        data[plotting]["min"],
        data[plotting]["max"],
        alpha=0.2,
        label=f"Min/Max {plotting}.",
    )

    axes[mat].plot(
        psize,
        data[plotting]["median"][0] / psize if plotting == "TTS" else psize,
        "-",
        label=f"Ideal {plotting}.",
    )

labels = [psize[0]]
for label in psize:
    if label >= labels[-1] * 2:
        labels.append(label)

del labels[-1]
labels.append(psize[-1])
print(labels)

for name, ax in axes.items():
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_ylabel(f"""{plotting} [{"s" if plotting == "TTS" else "-"}]""")
    ax.set_xticks(ticks=labels, labels=[str(label) for label in labels])

axes["cohesive"].set_title("With cohesive insertion")
axes["elastic"].set_title("Without cohesive insertion")

plt.xlabel("Nb Cores [-]")
plt.legend(bbox_to_anchor=(1.04, 2), loc="center left", borderaxespad=0)

plt.subplots_adjust(right=0.7)
plt.savefig(f"{plotting}.svg", bbox_inches="tight")
plt.show()

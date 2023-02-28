#!/usr/bin/env python3

import pandas as pd
import numpy as np

# needed to generate the plots on jed
import matplotlib

# matplotlib.use("TKAgg")
import matplotlib.pyplot as plt

# Same font as JOSS
plt.rcParams["font.sans-serif"] = "cmss10"

# Loading data
plots = {
    "elastic_gcc": {"prefix": "timmings_", "suffix": "_jed"},
    "cohesive_gcc": {"prefix": "timmings_", "suffix": "_jed"},
}

fig, ax = plt.subplots(figsize=(4.5, 4))
# fig, ax = plt.subplots(1, 1)

plotting = "TTS"

handles = []
for plot_name, data in plots.items():
    material, compiler = plot_name.split("_")
    data["df"] = pd.read_csv(
        f"""{data["prefix"]}{plot_name}{data["suffix"]}.csv""",
        sep=",",
        skipinitialspace=True,
    )
    data["material"] = material
    data["compiler"] = compiler

    df = data["df"]
    step = df["solve_step"] * df["solve_step nb_rep"]
    if material == "cohesive":
        step = step + df["check_cohesive_stress"] * df["check_cohesive_stress nb_rep"]

    df["TTS"] = step
    df["speedup"] = step[0] / step
    df["mumps"] = df["static_solve"] * df["static_solve nb_rep"]


def plot_measure(ax, df, plotting, label, **kwargs):
    """Plot a given measure."""
    grouped = df.groupby("psize")  # compute stats grouped by number of procs
    med = grouped.median()
    min = grouped.min()
    max = grouped.max()
    print(grouped)

    (l,) = ax.plot(med.index, med[plotting], label=f"{label} (median)", **kwargs)

    ax.fill_between(
        med.index, min[plotting], max[plotting], color=l.get_color(), alpha=0.2
    )

    ax.plot(med.index, med[plotting][1] / med.index, ls="--", color=l.get_color())
    # ax.boxplot(
    #     data[plotting]["grouped"], positions=psize, widths=[0.1 * s for s in psize]
    # )


plot_measure(
    ax,
    plots["cohesive_gcc"]["df"],
    plotting,
    "insertion",
    marker="o",
)

plot_measure(
    ax,
    plots["elastic_gcc"]["df"],
    plotting,
    "no insertion",
    marker="o",
)


# Selecting appropriate tick values
psize = np.array(np.unique(plots[list(plots.keys())[0]]["df"]["psize"]))
labels = np.concatenate(
    [[psize[0]], psize[1:][psize[1:] >= 2 * psize[:-1]], [psize[-1]]]
)
# for name, ax in axes.items():
ax.set_xscale("log", base=2)
ax.set_yscale("log")

ylabel = plotting if plotting != "TTS" else "Time to solution"
yunit = "s" if plotting != "speedup" else "-"

ax.set_xlabel("Nb Cores [-]")
ax.set_ylabel(f"""{ylabel} [{yunit}]""")

ax.set_xticks(ticks=labels, labels=map(str, labels))

# Constructing legend with min/max and ideal labels
handles, labels = ax.get_legend_handles_labels()
handles += [
    matplotlib.lines.Line2D([], [], linestyle="--", color="k"),
    matplotlib.patches.Patch(color="k", alpha=0.2),
]
labels += ["ideal", "min/max"]
ax.legend(handles=handles, labels=labels)

fig.tight_layout()
fig.savefig(f"{plotting}.svg", transparent=True, bbox_inches="tight", pad_inches=0.1)
# plt.show()

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Same font as JOSS
plt.rcParams["font.sans-serif"] = "cmss10"

# Loading data
elastic = pd.read_csv("timmings_elastic_gcc_jed.csv", sep=",",
                      skipinitialspace=True)
cohesive = pd.read_csv("timmings_cohesive_gcc_jed.csv", sep=",",
                       skipinitialspace=True)

# Solve time
elastic["step"] = elastic["solve_step"]
cohesive["step"] = cohesive["check_cohesive_stress"] + cohesive["solve_step"]

fig, ax = plt.subplots(figsize=(4, 3.5))


def plot_tts(ax, df, **kwargs):
    g = df.groupby("psize")  # compute stats grouped by number of procs
    med = g.median()
    min = g.min()
    max = g.max()
    l, = ax.plot(med.index, med["step"], **kwargs)
    ax.fill_between(med.index,
                    min["step"],
                    max["step"],
                    color=l.get_color(),
                    alpha=.2)
    ax.plot(med.index, med["step"][1] / med.index, ls='--', color=l.get_color())


plot_tts(ax, cohesive, marker="s", label="insertion (median)")
plot_tts(ax, elastic, marker="o", label="no insertion (median)")

ax.set_xscale("log", base=2)
ax.set_yscale("log")

ax.set_xlabel("Nb Cores [-]")
ax.set_ylabel("Time to Solution [s]")

# Selecting appropriate tick values
psize = np.sort(np.unique(elastic["psize"]))
psize = np.concatenate([[psize[0]],
                        psize[1:][psize[1:] >= 2 * psize[:-1]],
                        [psize[-1]]])
ax.set_xticks(ticks=psize, labels=map(str, psize))

# Constructing legend with min/max and ideal labels
handles, labels = ax.get_legend_handles_labels()
handles += [
    Line2D([], [], linestyle='--', color='k'),
    Patch(color='k', alpha=.2),
]
labels += ["ideal", "min/max"]
ax.legend(handles=handles, labels=labels)

fig.tight_layout()
fig.savefig("TTS.svg", transparent=True, bbox_inches="tight", pad_inches=0.1)

plt.show()

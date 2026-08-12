#!/usr/bin/env python3
# Real-data demo: build a protein-protein interaction (PPI) network from the
# STRING database and rank hub genes by degree centrality (the CytoHubba idea).
# Downloads a live network from STRING when possible; falls back to a bundled
# real subset (the TP53 neighbourhood) so the demo always runs.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from string_utils import (fetch_network, load_bundled, nodes_and_degree,
                          spring_layout)

SEED_GENES = ["TP53"]
BUNDLED = "data/string_network.tsv"


def load_edges():
    if os.environ.get("STRING_FORCE_BUNDLED") != "1":
        try:
            print("Fetching real PPI network from STRING ...")
            edges = fetch_network(SEED_GENES, required_score=700, limit=40)
            if len(edges) >= 20:
                print(f"Live download: {len(edges)} interactions")
                return edges, "live STRING download"
        except Exception as e:
            print(f"Live download failed ({e}); using bundled network.")
    else:
        print("Forced offline; using bundled real network.")
    return load_bundled(BUNDLED), "bundled real subset"


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    edges, source = load_edges()
    nodes, deg = nodes_and_degree(edges)
    pos, idx = spring_layout(nodes, edges)
    degrees = np.array([deg[n] for n in nodes])
    hubs = nodes[:10]

    fig, ax = plt.subplots(2, 2, figsize=(13, 10))

    # Panel 1: the network. Node size/colour = degree; top hubs labelled.
    a = ax[0, 0]
    for u, v, s in edges:
        xa, ya = pos[idx[u]]; xb, yb = pos[idx[v]]
        a.plot([xa, xb], [ya, yb], color="#CCCCCC", lw=0.4 + s * 0.8, zorder=1)
    sc = a.scatter(pos[:, 0], pos[:, 1], s=40 + degrees ** 1.8 * 6,
                   c=degrees, cmap="viridis", zorder=2, edgecolor="white", linewidth=0.5)
    for n in hubs[:8]:
        a.text(pos[idx[n], 0], pos[idx[n], 1], n, fontsize=8, fontweight="bold",
               ha="center", va="center", zorder=3)
    a.set_title(f"PPI network ({len(nodes)} genes, {len(edges)} interactions)")
    a.set_xticks([]); a.set_yticks([])
    fig.colorbar(sc, ax=a, fraction=0.046, label="degree")

    # Panel 2: hub genes by degree centrality (the CytoHubba ranking).
    a = ax[0, 1]
    top = hubs[::-1]
    a.barh(top, [deg[n] for n in top], color="#4C72B0")
    a.set_xlabel("degree (number of interactions)")
    a.set_title("Top hub genes by degree centrality")

    # Panel 3: degree distribution — PPI networks are scale-free-ish.
    a = ax[1, 0]
    a.hist(degrees, bins=range(1, degrees.max() + 2), color="#55A868", edgecolor="white", align="left")
    a.set_xlabel("degree"); a.set_ylabel("number of genes")
    a.set_title("Degree distribution (few hubs, many leaves)")

    # Panel 4: edge confidence-score distribution.
    a = ax[1, 1]
    a.hist([s for _, _, s in edges], bins=20, color="#8172B3", edgecolor="white")
    a.set_xlabel("STRING combined confidence score"); a.set_ylabel("edges")
    a.set_title("Interaction confidence scores")

    fig.suptitle("STRING PPI network & hub-gene analysis (real data)", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig("figures/demo.png", dpi=120)

    pd.DataFrame({
        "gene": nodes, "degree": degrees,
        "rank": np.arange(1, len(nodes) + 1),
        "is_hub_top10": [n in hubs for n in nodes],
    }).to_csv("results/summary.csv", index=False)

    print(f"Source: {source}")
    print(f"Network: {len(nodes)} genes, {len(edges)} interactions")
    print(f"Top hub genes: {', '.join(f'{n}({deg[n]})' for n in hubs[:5])}")
    print("Wrote figures/demo.png and results/summary.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Helpers for a STRING protein-protein interaction (PPI) network: live download,
# degree centrality (hub genes), and a numpy-only force-directed layout.
import json
import sys
import urllib.request
import urllib.parse
import numpy as np

API = "https://string-db.org/api/tsv/network"


def fetch_network(genes, species=9606, required_score=700, limit=40, timeout=30):
    # Download a real PPI network from STRING. Returns list of (a, b, score).
    params = urllib.parse.urlencode({
        "identifiers": "%0d".join(genes), "species": species,
        "required_score": required_score, "limit": limit})
    with urllib.request.urlopen(f"{API}?{params}", timeout=timeout) as r:
        text = r.read().decode("utf-8")
    edges = []
    for line in text.splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 6:
            edges.append((f[2], f[3], float(f[5])))
    return edges


def load_bundled(path):
    edges = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            a, b, s = line.rstrip("\n").split("\t")
            edges.append((a, b, float(s)))
    return edges


def nodes_and_degree(edges):
    deg = {}
    for a, b, _ in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    nodes = sorted(deg, key=lambda n: deg[n], reverse=True)
    return nodes, deg


def spring_layout(nodes, edges, iterations=200, seed=0):
    # Fruchterman-Reingold force-directed layout in pure numpy.
    rng = np.random.default_rng(seed)
    n = len(nodes)
    idx = {node: i for i, node in enumerate(nodes)}
    pos = rng.uniform(-1, 1, size=(n, 2))
    E = np.array([(idx[a], idx[b]) for a, b, _ in edges])
    k = 1.0 / np.sqrt(n)
    t = 0.1
    for _ in range(iterations):
        disp = np.zeros((n, 2))
        diff = pos[:, None, :] - pos[None, :, :]
        dist = np.sqrt((diff ** 2).sum(-1)) + 1e-9
        rep = (k * k / dist)[..., None] * (diff / dist[..., None])
        disp += rep.sum(axis=1)
        d = pos[E[:, 0]] - pos[E[:, 1]]
        dd = np.sqrt((d ** 2).sum(1)) + 1e-9
        att = (dd * dd / k)[:, None] * (d / dd[:, None])
        np.add.at(disp, E[:, 0], -att)
        np.add.at(disp, E[:, 1], att)
        length = np.sqrt((disp ** 2).sum(1)) + 1e-9
        pos += (disp / length[:, None]) * np.minimum(length, t)[:, None]
        t *= 0.99
    return pos, idx


if __name__ == "__main__":
    e = load_bundled("data/string_network.tsv")
    nodes, deg = nodes_and_degree(e)
    print("edges:", len(e), "nodes:", len(nodes))
    print("top hubs:", [(n, deg[n]) for n in nodes[:5]])

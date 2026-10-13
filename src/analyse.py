"""
Interpret the screen: compare observed matches against what chance alone would
give, normalise per megabase of non-target sequence, pick a construct region
for each gene, and draw the match map.
"""
import json, math, os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
MIN_CONSTRUCT = 200          # nt; a workable hairpin arm
SHORT = {"Apis mellifera (honey bee)": "honey bee",
         "Danaus plexippus (monarch butterfly)": "monarch",
         "Homo sapiens (human)": "human"}

d = json.load(open(os.path.join(RES, "screen_results.json")))
sizes = d["non_target_sizes_nt"]
total_nt = sum(sizes.values())


def expected_flagged(n_windows, k, nt):
    """Windows expected to match by chance: both strands, uniform base usage."""
    p_site = 2 * nt / (4 ** k)
    return n_windows * (1 - math.exp(-p_site))


def gc(s):
    return 100.0 * sum(c in "GC" for c in s) / len(s) if s else 0.0


rows = []
print(f"{'gene':<16}{'k':>4}{'obs':>7}{'chance':>9}{'ratio':>8}   per-species per Mnt")
for g in d["genes"]:
    label = g["header"].split("(")[0].split("PREDICTED: ")[-1].strip()
    for k in d["k_values"]:
        t = g["tiers"][str(k)]
        n = t["windows_tested"]
        obs = t["windows_flagged"]
        exp = expected_flagged(n, k, total_nt)
        ratio = obs / exp if exp > 0 else float("nan")
        rates = {SHORT[s]: round(len(p) / (sizes[s] / 1e6), 3)
                 for s, p in t["flagged_positions_by_species"].items()}
        print(f"{g['accession']:<16}{k:>4}{obs:>7}{exp:>9.1f}{ratio:>8.1f}   {rates}")
        rows.append({"gene": label, "accession": g["accession"], "length_nt": g["length"],
                     "k": k, "windows_tested": n, "windows_flagged": obs,
                     "expected_by_chance": round(exp, 2), "observed_over_chance": round(ratio, 2),
                     **{f"flagged_per_Mnt_{v}": rates.get(v, 0.0) for v in SHORT.values()}})

with open(os.path.join(RES, "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ---- pick a construct per gene at the regulatory k=21 criterion ----
print("\nRecommended dsRNA construct regions (no 21 nt match in any non-target):")
fasta = []
picks = []
for g in d["genes"]:
    t = g["tiers"]["21"]
    regions = [r for r in t["longest_safe_regions"] if r["length"] >= MIN_CONSTRUCT]
    if not regions:
        print(f"  {g['accession']}: no clean region >= {MIN_CONSTRUCT} nt")
        continue
    r = regions[0]
    seq = g["sequence"][r["start_1based"] - 1:r["end_1based"]]
    picks.append((g, r, seq))
    print(f"  {g['accession']}: {r['length']} nt at {r['start_1based']}-{r['end_1based']}, GC {gc(seq):.1f}%")
    fasta.append(f">{g['accession']}_safe_{r['start_1based']}_{r['end_1based']} "
                 f"len={r['length']} GC={gc(seq):.1f}%\n" +
                 "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)))

with open(os.path.join(RES, "candidate_dsRNA_regions.fasta"), "w") as fh:
    fh.write("\n".join(fasta) + "\n")

# ---- figure ----
fig, axes = plt.subplots(len(d["genes"]), 1, figsize=(11, 2.6 * len(d["genes"])))
colours = {"honey bee": "#E6A700", "monarch": "#C2410C", "human": "#1F6FEB"}
for ax, g in zip(axes, d["genes"]):
    t = g["tiers"]["21"]
    name = g["header"].split("PREDICTED: ")[-1].split(" (LOC")[0].split(" mRNA")[0].replace(g["accession"] + " ", "")
    ax.set_xlim(0, g["length"]); ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.plot([0, g["length"]], [0.22, 0.22], color="#D0D7DE", lw=7, solid_capstyle="butt")
    for i, (sp, short) in enumerate(SHORT.items()):
        pos = t["flagged_positions_by_species"].get(sp, [])
        y = 0.5 + i * 0.16
        if pos:
            ax.vlines(pos, y - 0.055, y + 0.055, color=colours[short], lw=0.8)
        ax.text(g["length"] * 1.005, y, f"{short} ({len(pos)})", va="center",
                fontsize=8, color=colours[short])
    for (gg, r, _s) in picks:
        if gg["accession"] == g["accession"]:
            ax.plot([r["start_1based"], r["end_1based"]], [0.22, 0.22],
                    color="#1A7F37", lw=7, solid_capstyle="butt")
            ax.text((r["start_1based"] + r["end_1based"]) / 2, 0.06,
                    f"construct region {r['length']} nt", ha="center", fontsize=8,
                    color="#1A7F37")
    ax.set_title(f"{name}  ({g['accession']}, {g['length']} nt)", fontsize=9, loc="left")
    ax.set_xlabel("position in transcript (nt)", fontsize=8)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
fig.suptitle("Where a 21 nt stretch of each locust target also occurs in a non-target transcriptome",
             fontsize=11, y=0.995)
fig.tight_layout(rect=[0, 0, 0.88, 0.97])
fig.savefig(os.path.join(RES, "offtarget_map.png"), dpi=170)
print("\nwrote summary.csv, candidate_dsRNA_regions.fasta, offtarget_map.png")

"""
Off-target screen for RNAi (dsRNA) pest-control designs.

Question: if we express a dsRNA against a locust gene, does any short stretch of
it also match a transcript of a species we do not want to harm?

Method: RNAi specificity is governed by short exact matches, so we take every
k-nucleotide window of each target transcript (both strands) and ask whether it
occurs anywhere in a non-target transcriptome. Positions with no match in any
non-target are candidate sites for a dsRNA construct.

k = 21 is the standard regulatory criterion (a shared >=21 nt exact match).
k = 16 is reported as a stricter tier; at that length chance matches are
expected in a large transcriptome, which is itself worth showing.
"""
import gzip, json, os, sys
from collections import defaultdict
import ahocorasick

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
KS = (21, 16)
COMP = str.maketrans("ACGT", "TGCA")

NON_TARGETS = [
    ("Apis mellifera (honey bee)", "bee.fna.gz"),
    ("Danaus plexippus (monarch butterfly)", "monarch.fna.gz"),
    ("Homo sapiens (human)", "human.fna.gz"),
]


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name:
        yield name, "".join(chunks)


def main():
    targets = [(n.split()[0], n, s.upper()) for n, s in
               read_fasta(os.path.join(DATA, "targets.fa"))]
    print(f"targets: {len(targets)}")

    # pattern -> list of (target_accession, start_position, k)
    index = defaultdict(list)
    for acc, _hdr, seq in targets:
        for k in KS:
            for i in range(len(seq) - k + 1):
                w = seq[i:i + k]
                if set(w) - set("ACGT"):
                    continue
                index[w].append((acc, i, k))
                index[rc(w)].append((acc, i, k))

    A = ahocorasick.Automaton()
    for pat, locs in index.items():
        A.add_word(pat, (pat, locs))
    A.make_automaton()
    print(f"patterns in automaton: {len(index)}")

    # hits[acc][k][start_position] = {species: count}
    hits = {acc: {k: defaultdict(lambda: defaultdict(int)) for k in KS}
            for acc, _h, _s in targets}

    sizes = {}
    for species, fname in NON_TARGETS:
        path = os.path.join(DATA, fname)
        n_rec = n_nt = n_hit = 0
        for _name, seq in read_fasta(path):
            seq = seq.upper()
            n_rec += 1
            n_nt += len(seq)
            for _end, (_pat, locs) in A.iter(seq):
                n_hit += 1
                for acc, pos, k in locs:
                    hits[acc][k][pos][species] += 1
        print(f"{species}: {n_rec:,} transcripts, {n_nt:,} nt, {n_hit:,} raw matches")
        sizes[species] = n_nt

    out = {"k_values": list(KS),
           "non_targets": [s for s, _ in NON_TARGETS],
           "non_target_sizes_nt": sizes,
           "genes": []}

    for acc, hdr, seq in targets:
        gene = {"accession": acc, "header": hdr, "sequence": seq, "length": len(seq), "tiers": {}}
        for k in KS:
            hit_pos = hits[acc][k]
            n_windows = len(seq) - k + 1
            flagged = sorted(hit_pos)
            per_species = defaultdict(int)
            for pos in flagged:
                for sp in hit_pos[pos]:
                    per_species[sp] += 1
            # safe runs: maximal stretches of consecutive unflagged window starts
            runs, start = [], None
            flagged_set = set(flagged)
            for i in range(n_windows):
                if i not in flagged_set:
                    if start is None:
                        start = i
                elif start is not None:
                    runs.append((start, i - 1))
                    start = None
            if start is not None:
                runs.append((start, n_windows - 1))
            # a construct spanning [s, e+k) contains only unflagged windows
            constructs = sorted(((e + k - s), s, e + k) for s, e in runs)[::-1]
            by_species_positions = defaultdict(list)
            for pos in flagged:
                for sp in hit_pos[pos]:
                    by_species_positions[sp].append(pos)
            gene["tiers"][str(k)] = {
                "flagged_positions": flagged,
                "flagged_positions_by_species": {s: v for s, v in by_species_positions.items()},
                "windows_tested": n_windows,
                "windows_flagged": len(flagged),
                "flagged_by_species": dict(per_species),
                "longest_safe_regions": [
                    {"length": L, "start_1based": s + 1, "end_1based": e}
                    for L, s, e in constructs[:5]
                ],
            }
        out["genes"].append(gene)

    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "screen_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    for g in out["genes"]:
        print(f"\n{g['accession']}  {g['length']} nt")
        for k in KS:
            t = g["tiers"][str(k)]
            print(f"  k={k}: {t['windows_flagged']}/{t['windows_tested']} windows flagged "
                  f"{t['flagged_by_species']}")
            best = t["longest_safe_regions"][0] if t["longest_safe_regions"] else None
            if best:
                print(f"        longest clean region: {best['length']} nt "
                      f"({best['start_1based']}-{best['end_1based']})")


if __name__ == "__main__":
    main()

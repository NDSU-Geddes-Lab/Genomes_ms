#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import SeqIO

# size thresholds (bp)
PSYMB_MIN = 1_400_000
PSYMB_MAX =   3_000_000
ACC_MIN   =      30_000
ACC_MAX   =     750_000

def find_gene(rec, name):
    """Return True if any feature has gene==name or product contains name (case-insensitive)."""
    name = name.lower()
    for f in rec.features:
        for tag in ("gene","product"):
            for val in f.qualifiers.get(tag, []):
                if name in val.lower():
                    return True
    return False

def has_all(rec, genes):
    return all(find_gene(rec, g) for g in genes)

def classify_one(path):
    path = Path(path)
    strain = path.stem
    recs = list(SeqIO.parse(path, "genbank"))
    assigned = {}

    # 1) chromosome = the single largest contig
    chrom = max(recs, key=lambda r: len(r.seq))
    assigned[chrom.id] = f"{strain}_chromosome"

    # 2) pSymA = any contig with nifA
    for r in recs:
        if r.id in assigned: continue
        if find_gene(r, "nifA"):
            assigned[r.id] = f"{strain}_pSymA"
            break

    # 3) pSymB = 1.4–3 Mb AND repA/B/C
    for r in recs:
        if r.id in assigned: continue
        L = len(r.seq)
        if PSYMB_MIN <= L <= PSYMB_MAX and has_all(r, ["repA","repB","repC"]):
            assigned[r.id] = f"{strain}_pSymB"
            break

    # 4) accessory = 30 kb–750 kb AND repA/B/C
    accs = [r for r in recs
            if r.id not in assigned
            and ACC_MIN <= len(r.seq) <= ACC_MAX
            and has_all(r, ["repA","repB","repC"])]
    for i, r in enumerate(accs, 1):
        assigned[r.id] = f"{strain}_accessory{i}"

    # 5) everything else → contigN
    leftover = [r for r in recs if r.id not in assigned]
    for i, r in enumerate(leftover, 1):
        assigned[r.id] = f"{strain}_contig{i}"

    # write out one GBK per replicon
    for r in recs:
        new_id = assigned[r.id]
        r.id = new_id
        r.name = new_id
        r.description = ""
        out = path.parent / f"{new_id}.gbk"
        SeqIO.write(r, out, "genbank")
        print(f"Wrote {out}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: annotate_gbk_contigs.py <GBK-file-or-dir> [more...]", file=sys.stderr)
        sys.exit(1)

    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.is_dir():
            # recurse and find every .gbk under this directory
            for gbk in p.rglob("*.gbk"):
                classify_one(gbk)
        elif p.is_file() and p.suffix.lower() == ".gbk":
            classify_one(p)
        else:
            print(f"Skipping {p}: not a .gbk file or directory", file=sys.stderr)

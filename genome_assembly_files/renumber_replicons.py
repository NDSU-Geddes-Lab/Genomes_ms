#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

def find_repA(rec):
    """Return the first feature whose gene/product qualifier contains 'repA'."""
    for f in rec.features:
        for tag in ("gene","product"):
            for val in f.qualifiers.get(tag, []):
                if "repa" in val.lower():
                    return f
    return None

def rotate_and_renumber(rec, origin):
    """
    Rotate sequence so that `origin` (0‑based) becomes new position 0,
    shift all features by that offset (with wrap), and return the modified rec.
    """
    seq = rec.seq
    L   = len(seq)
    # rotate sequence: tail then head
    new_seq = seq[origin:] + seq[:origin]

    new_feats = []
    for f in rec.features:
        s0 = int(f.location.start)
        e0 = int(f.location.end)
        # compute new coords
        if s0 >= origin:
            ns = s0 - origin
            ne = e0 - origin
        else:
            ns = s0 + L - origin
            ne = e0 + L - origin
        new_loc = FeatureLocation(ns, ne, strand=f.location.strand)
        new_feats.append(SeqFeature(new_loc, type=f.type, qualifiers=f.qualifiers))

    rec.seq      = new_seq
    rec.features = new_feats
    rec.annotations["comment"] = f"renumbered: original RepA start was {origin+1}"
    return rec

def process_path(inp, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    if inp.is_dir():
        paths = list(inp.rglob("*.gbk"))
    elif inp.is_file() and inp.suffix.lower()==".gbk":
        paths = [inp]
    else:
        print(f"[ERROR] {inp} is not a .gbk file or directory", file=sys.stderr)
        sys.exit(1)

    for gbk in paths:
        for rec in SeqIO.parse(gbk, "genbank"):
            # only target pSymA/pSymB/accessory records
            if not any(x in rec.id for x in ("pSymA","pSymB","accessory")):
                continue

            repA = find_repA(rec)
            if repA is None:
                print(f"[WARN] no repA in {rec.id}", file=sys.stderr)
                continue

            origin = int(repA.location.start)
            rec2   = rotate_and_renumber(rec, origin)

            outp = outdir / f"{rec2.id}.gbk"
            SeqIO.write(rec2, outp, "genbank")
            print(f"Wrote renumbered → {outp}")

if __name__=="__main__":
    if len(sys.argv)!=3:
        print("Usage: renumber_replicons.py <gbk-file-or-dir> <output-dir>", file=sys.stderr)
        sys.exit(1)

    inp    = Path(sys.argv[1])
    outdir = Path(sys.argv[2])
    process_path(inp, outdir)

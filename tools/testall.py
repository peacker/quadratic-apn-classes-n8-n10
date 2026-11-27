#!/usr/bin/env python3
"""Point 1: test quadratic APN functions for CCZ-equivalence to a permutation, and record class invariants.

Input: text files with one lookup table per line (any separators, e.g. "{0, 1, 8, ...}").
Output (TSV, appended, resumable): source  line  apn  inv  M  status  seconds
  inv     hash of the ortho-derivative invariant (differential + extended Walsh spectrum of the
          ortho-derivative): a CCZ-invariant of quadratic APN functions, used to key the database
  M       largest dimension of a subspace of the non-bent components together with 0
  status  'none' (not CCZ-equivalent to a permutation) or 'permutation' (from the paper's test:
          counting criterion if M < n/2, otherwise the dual test of Theorem 8)
Sharding: --shard i --nshards k processes lines with (line % k) == i.

usage: testall.py --n 8 --source NAME --input FILE --out OUT.tsv [--shard i --nshards k]
"""
import argparse, hashlib, os, re, sys, time
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "lib"))
from apnsearch import ccz, fast  # noqa: E402
from apnsearch.invariants import ortho_invariant  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, required=True)
ap.add_argument("--source", required=True)
ap.add_argument("--input", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--shard", type=int, default=0)
ap.add_argument("--nshards", type=int, default=1)
ap.add_argument("--limit", type=int, default=0, help="stop after this many functions (testing)")
ap.add_argument("--keep-lut", action="store_true", help="append the lookup table (first occurrence of each invariant)")
a = ap.parse_args()
N = 1 << a.n

done = set()
if os.path.exists(a.out):
    for line in open(a.out):
        t = line.split("\t")
        if len(t) >= 7 and t[0] == a.source:
            done.add(int(t[1]))

num = re.compile(r"\d+")
seen_inv = set()
if a.keep_lut and os.path.exists(a.out):
    seen_inv = {l.split("\t")[3] for l in open(a.out) if l.count("\t") >= 6}
out = open(a.out, "a")
count = 0
with open(a.input) as f:
    for i, line in enumerate(f):
        if i % a.nshards != a.shard or i in done:
            continue
        vals = num.findall(line)
        if len(vals) != N:
            continue
        F = np.array(list(map(int, vals)), dtype=np.int64)
        t0 = time.time()
        apn = fast.is_apn(F)
        inv = M = st = "-"
        if apn:
            try:
                inv = hashlib.sha1(repr(ortho_invariant(F, a.n)).encode()).hexdigest()[:16]
                quad = True
            except ValueError:  # APN but not quadratic: no ortho-derivative invariant
                inv, quad = "nonquadratic", False
            if quad:
                st, info = ccz.ccz_to_permutation_dual(F, a.n, quadratic=True)
                M = info.get("max_dim_NB")
            else:  # counting criterion only; a full decision needs kasami/run_wzs.py
                st, info = "untested", {}
                try:
                    st, info = ccz.ccz_to_permutation_dual(F, a.n, quadratic=False)
                except Exception as e:  # non-affine supports: the dual test does not apply
                    st = "needs_wzs"
                M = info.get("max_dim_NB", "-")
        lut = ""
        if a.keep_lut and apn and inv not in seen_inv:
            seen_inv.add(inv)
            lut = "\t" + ",".join(map(str, F.tolist()))
        out.write(f"{a.source}\t{i}\t{int(apn)}\t{inv}\t{M}\t{st}\t{time.time() - t0:.3f}{lut}\n")
        if st == "permutation":
            open(a.out + ".PERMUTATION", "a").write(f"{a.source}\t{i}\t{line}")
        count += 1
        if count % 200 == 0:
            out.flush()
        if a.limit and count >= a.limit:
            break
out.close()

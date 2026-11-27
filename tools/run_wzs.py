#!/usr/bin/env python3
"""Decide whether an (n,n)-function is CCZ-equivalent to a permutation, for Walsh supports that need not be
affine, by listing the Walsh zero spaces U with n - M <= dim pi(U) <= M (Proposition 7 of our paper) with
wzs.c and looking for a pair with trivial intersection (Canteaut-Perrin criterion, Prop. 3 of the paper).

functions: kim6 (n = 6, must give a permutation: Dillon), gold10 (x^3, n = 10, must give none),
           kasami10 (x^57, n = 10, open), dobbertin10 (x^339; not plateaued), or a LUT file.
usage: run_wzs.py FUNC [--jobs J]
"""
import argparse, itertools, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "lib"))
from apnsearch import fast  # noqa: E402
from apnsearch.gf import GF  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("func")
ap.add_argument("--jobs", type=int, default=1)
ap.add_argument("--alldims", action="store_true", help="list Walsh zero spaces of every projection dimension (validation)")
a = ap.parse_args()

if a.func == "kim6":
    from apnsearch import families
    n, F = 6, families.kim6()[1]
elif a.func in ("gold10", "kasami10", "dobbertin10"):
    n = 10
    F = GF(10).power_table({"gold10": 3, "kasami10": 57, "dobbertin10": 339}[a.func])
else:
    F = np.array(list(map(int, open(a.func).read().replace(",", " ").split())))
    n = int(len(F)).bit_length() - 1
F = np.asarray(F, dtype=np.int64)
N = 1 << n
assert fast.is_apn(F), "not APN"

# Walsh transform: W[a, b]
par = np.array([bin(v).count("1") & 1 for v in range(N)])
M_ = 1 - 2 * par[np.bitwise_and.outer(F, np.arange(N))]   # (x, b)
h = 1
W = M_.astype(np.int64)
while h < N:
    W = W.reshape(N // (2 * h), 2, h, -1)
    W = np.concatenate([W[:, 0] + W[:, 1], W[:, 0] - W[:, 1]], axis=1).reshape(N, -1)
    h *= 2
nb = np.zeros(N, dtype=np.uint8)
for b in range(1, N):
    if np.any(W[:, b] == 0):
        nb[b] = 1
amps = {b: sorted(set(np.abs(W[:, b]).tolist()) - {0}) for b in range(1, N) if nb[b]}
plateaued = all(len(v) == 1 for v in amps.values())
print(f"{a.func}: n = {n}, non-bent components {int(nb.sum())}, plateaued non-bent components: {plateaued}")

# subspaces of NB u {0}
M = 0
Bs = []
for d in range(n, 0, -1):
    S = fast.find_subspaces(nb, d)
    if S is not None and len(S):
        if not M:
            M = d
        if d < n - M and not a.alldims:
            break
        Bs += [list(map(int, s)) for s in S]
print(f"M = {M}; subspaces B with {n - M} <= dim <= {M}: {len(Bs)}")
if M < n / 2 and not a.alldims:
    print("RESULT: not CCZ-equivalent to a permutation (counting criterion)")
    sys.exit(0)

tmp = os.path.join(here, f"_{a.func.replace('/', '_')}")
os.makedirs(tmp, exist_ok=True)
open(os.path.join(tmp, "lut.txt"), "w").write(",".join(map(str, F.tolist())))
open(os.path.join(tmp, "B.txt"), "w").write("".join(" ".join(map(str, B)) + "\n" for B in Bs))
exe = os.path.join(here, "wzs")
t0 = time.time()


def shard(i):
    r = subprocess.run([exe, str(n), os.path.join(tmp, "lut.txt"), os.path.join(tmp, "B.txt"), str(i), str(a.jobs)],
                       capture_output=True, text=True)
    return r.stdout


with ThreadPoolExecutor(a.jobs) as ex:
    outs = list(ex.map(shard, range(a.jobs)))
U = [list(map(int, l.split()[2:])) for o in outs for l in o.splitlines() if l.startswith("U ")]
if a.alldims:
    U.append(list(range(0, 0)) + [1 << i for i in range(n)])  # pi(U) = {0}: U = F_2^n x {0}
print(f"Walsh zero spaces with admissible projection: {len(U)} ({time.time() - t0:.1f}s)")


def rank(vs):
    basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


assert all(rank(u) == n for u in U), "a listed space is not n-dimensional"
# every nonzero vector (a, b) of every listed space must be a Walsh zero: W(a, b) = 0
Wz = (W == 0)
Wz[0, 0] = True
for u in U:
    span = [0]
    for v in u:
        span += [x ^ v for x in span]
    pts = np.array(span[1:])
    assert Wz[pts & (N - 1), pts >> n].all(), "a listed space contains a nonzero Walsh value"
print("check: every listed space consists of Walsh zeros")
pair = None
for i, j in itertools.combinations_with_replacement(range(len(U)), 2):
    if rank(U[i] + U[j]) == 2 * n:
        pair = (i, j)
        break
if pair:
    print(f"RESULT: CCZ-equivalent to a permutation (Walsh zero spaces {pair[0]} and {pair[1]} intersect trivially)")
    open(os.path.join(tmp, "pair.txt"), "w").write(f"{U[pair[0]]}\n{U[pair[1]]}\n")
else:
    print("RESULT: not CCZ-equivalent to a permutation (no pair of Walsh zero spaces with trivial intersection)")

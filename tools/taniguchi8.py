#!/usr/bin/env python3
"""All members of the Taniguchi family at n = 8 (bivariate on F_16 x F_16, family T of Göloğlu-Kölsch Table 2.1):
F(x, y) = (x^(q+1) + x y^q + d y^(q+1), x y^(q^2)), q = 2^k, gcd(k, 4) = 1, x^(q+1) + x + d != 0 on F_16.
Writes the APN members, one lookup table per line, (x, y) <-> x + 16 y."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from apnsearch import fast
from apnsearch.gf import GF
G = GF(4)
X = np.arange(256); xs, ys = X & 15, X >> 4
p = lambda V, e: G.pow_vec(V, e % 15 or 15) if e else np.ones_like(V)
m = lambda A, B: G.mul_vec(A, B)
c = lambda a, V: G.mul_vec(np.full(len(V), a), V)
out = []
for k in (1, 3):
    q = 1 << k
    for d in range(16):
        V = np.arange(16)
        if np.any((p(V, q + 1) ^ V ^ d) == 0):
            continue
        f = p(xs, q + 1) ^ m(xs, p(ys, q)) ^ c(d, p(ys, q + 1))
        g = m(xs, p(ys, q * q))
        F = (f | (g << 4)).astype(np.int64)
        if fast.is_apn(F):
            out.append(F)
open(sys.argv[1], "w").write("".join(",".join(map(str, F.tolist())) + "\n" for F in out))
print(len(out), "APN members")

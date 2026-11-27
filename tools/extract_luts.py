#!/usr/bin/env python3
"""Turn the sboxU 8-bit APN data (Python modules with lists of lookup tables) and the Edel-Pott polynomial
list into plain files with one lookup table per line, for testall.py.

usage: extract_luts.py SBOXU_APNDB_DIR OUTDIR
  writes OUTDIR/{sboxu_firstQAM,sboxu_secondQAM,sboxu_WenTanGon,sboxu_EdelPott}.txt
Edel-Pott polynomials (sboxU eightBitAPN.all_quadratic_polynomials, Table 9 of Edel-Pott) are evaluated in
GF(2^8) = GF(2)[a]/(a^8+a^4+a^3+a^2+1), the field Sage uses for GF(2**8) by default, with g = a.
"""
import os, re, sys

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
lut_re = re.compile(r"\[\s*(\d+(?:\s*,\s*\d+){255})\s*,?\s*\]")
hex_re = re.compile(r"\[\s*(0x[0-9a-fA-F]+(?:\s*,\s*0x[0-9a-fA-F]+){255})\s*,?\s*\]")


def extract(path, dst):
    text = open(path).read()
    luts = [m.group(1) for m in lut_re.finditer(text)]
    luts += [",".join(str(int(t, 16)) for t in m.group(1).split(",")) for m in hex_re.finditer(text)]
    with open(dst, "w") as f:
        for l in luts:
            f.write(re.sub(r"\s+", "", l) + "\n")
    return len(luts)


print("firstQAM", extract(os.path.join(src, "firstQAM.py"), os.path.join(out, "sboxu_firstQAM.txt")))
print("secondQAM", extract(os.path.join(src, "secondQAM.py"), os.path.join(out, "sboxu_secondQAM.txt")))
# Wen-Tan-Gong: the LUT lists inside all_WenTanGon()
txt = open(os.path.join(src, "eightBitAPN.py")).read()
wtg = txt[txt.index("def all_WenTanGon"):txt.index("def all_quadratics")]
luts = [m.group(1) for m in lut_re.finditer(wtg)]
open(os.path.join(out, "sboxu_WenTanGon.txt"), "w").write("".join(re.sub(r"\s+", "", l) + "\n" for l in luts))
print("WenTanGon", len(luts))

# Edel-Pott polynomials
MOD = 0x11d


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 0x100:
            a ^= MOD
    return r


def pw(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


ep = txt[txt.index("def all_quadratic_polynomials"):txt.index("def first_QAMs")]
polys = re.findall(r"poly_to_lut\((.*?)\),\n", ep)
rows = []
for p in polys:
    p = p.replace("+   +", "+")
    terms = []  # (coefficient, exponent)
    if "sum(" in p:  # X**a + sum(X**(b*2**i) for i in range(0, N))
        a_, b_ = re.findall(r"X\*\*\(?(\d+)", p)[:2]
        terms = [(1, int(a_))] + [(1, int(b_) * 2 ** i % 255 or 255) for i in range(8)]
    else:
        expr = p
        m = re.match(r"X\*\*(\d+) \+ X\*\*(\d+) \+ g\*\*(\d+)\*\(X\*\*(\d+) \+ X\*\*(\d+)\) \+ g\*\*(\d+)\*X\*\*(\d+)$", expr.strip())
        if m:
            e1, e2, c, e3, e4, c2, e5 = map(int, m.groups())
            terms = [(1, e1), (1, e2), (pw(2, c), e3), (pw(2, c), e4), (pw(2, c2), e5)]
        else:
            for t in expr.split("+"):
                t = t.strip()
                mc = re.match(r"(?:g(?:\*\*(\d+))?\*)?X\*\*(\d+)$", t)
                if not mc:
                    raise ValueError(f"cannot parse term {t!r} in {p!r}")
                cexp = mc.group(1)
                coef = 1 if t.startswith("X") else pw(2, int(cexp) if cexp else 1)
                terms.append((coef, int(mc.group(2))))
    lut = []
    for x in range(256):
        y = 0
        for c, e in terms:
            y ^= mul(c, pw(x, e)) if x else (0 if e else c)
        lut.append(y)
    rows.append(lut)
open(os.path.join(out, "sboxu_EdelPott.txt"), "w").write("".join(",".join(map(str, l)) + "\n" for l in rows))
print("EdelPott", len(rows))

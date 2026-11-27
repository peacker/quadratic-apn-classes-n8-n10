"""Decode the quadratic functions of the database of Baudrin, Galissant, Perrin (ePrint 2026/2799,
Zenodo 10.5281/zenodo.24006393). The blob is the sboxU quadratic_compact_representation: for each output
coordinate k, for each pair i < j (lexicographic), one bit for the monomial x_i x_j; read as a little-endian int.
Output TSV: id  ccz_id  source  aut_size  inv  lut
usage: decode_bgp26.py 8_bit_compact.db 8 "SQL condition" OUT.tsv
(the 3,833 new n = 8 classes: condition "source=7")"""
import os, sqlite3, sys, hashlib
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from apnsearch.invariants import ortho_invariant
def decode(blob,n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]; v=int.from_bytes(blob,'little'); t=0
    coef={p:0 for p in pairs}
    for k in range(n):
        for p in pairs:
            if (v>>t)&1: coef[p]|=1<<k
            t+=1
    F=[0]*(1<<n)
    for x in range(1<<n):
        r=0
        for (i,j),c in coef.items():
            if c and (x>>i)&1 and (x>>j)&1: r^=c
        F[x]=r
    return F
def isapn(F,n):
    N=1<<n
    return all(len({F[x]^F[x^a] for x in range(N)})==N//2 for a in range(1,N))
db,n,where,out=sys.argv[1],int(sys.argv[2]),sys.argv[3],sys.argv[4]
rows=sqlite3.connect(db).execute(f"select id,ccz_id,source,aut_size,representation from functions where type=0 and {where}").fetchall()
with open(out,'w') as f:
    for i,(rid,cz,src,aut,b) in enumerate(rows):
        F=decode(b,n)
        if i<3: assert isapn(F,n), rid
        h=hashlib.sha1(repr(ortho_invariant(F,n)).encode()).hexdigest()[:16]
        f.write(f"{rid}\t{cz}\t{src}\t{aut}\t{h}\t{','.join(map(str,F))}\n")
print(len(rows))

import sys, hashlib
for l in sys.stdin:
    i, st, r = l.rstrip("\n").split("\t")
    k = hashlib.sha1(r.encode()).hexdigest()[:16] if st == "Q" else "-"
    sys.stdout.write(f"{i}\t{st}\t{k}\n")

/*
 * wzs: Walsh zero spaces of an (n,n)-function with a given projection B, for functions whose Walsh supports
 * need not be affine (e.g. the Kasami function x^57 in dimension 10).
 *
 * Theorem 8 of our paper, parts (1)-(2), holds for every function: each n-dimensional subspace U of F_2^{2n}
 * with pi(U) = B (dim m) is U(B, P, phi) = {(a, b) : b in B, a in phi(b) + P^perp} for a unique m-dimensional
 * P (= U_0^perp) and a unique linear map phi : B -> F_2^n / P^perp. U is a Walsh zero space iff for every
 * b in B \ {0} the coset phi(b) + P^perp avoids the Walsh support S_b of the component b. (For b = 0 the fibre
 * P^perp consists of zeros automatically.) Only part (3) of the theorem needs affine supports; here it is
 * replaced by this direct check.
 *
 * The coset space F_2^n / P^perp is identified with F_2^m through c(a) = (y_1.a, ..., y_m.a), y_i a basis of P.
 *
 * Input: n, a lookup table file (2^n values), and a file of subspaces B (one per line: m basis vectors).
 * The m-dimensional subspaces P are enumerated through their reduced row echelon forms; --shard i --nshards k
 * splits that enumeration. Output: one line per Walsh zero space, "U <B index> <2n-bit basis vectors>",
 * where a vector (a, b) is written as a + (b << n).
 *
 * usage: wzs N LUTFILE BFILE [SHARD NSHARDS]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12
static int n, N;
static int *F;
static uint8_t *zero;          /* zero[b * N + a] = 1 iff W_F(a, b) = 0 */
static int **supp, *supplen;   /* Walsh support of each component b (a with W != 0) */

static void walsh(void) {
    zero = calloc((size_t)N * N, 1);
    supp = calloc(N, sizeof(int *));
    supplen = calloc(N, sizeof(int));
    int *w = malloc(sizeof(int) * N);
    for (int b = 1; b < N; b++) {
        for (int x = 0; x < N; x++) w[x] = (__builtin_popcount(b & F[x]) & 1) ? -1 : 1;
        for (int h = 1; h < N; h <<= 1)
            for (int i = 0; i < N; i += h << 1)
                for (int j = i; j < i + h; j++) { int u = w[j], v = w[j + h]; w[j] = u + v; w[j + h] = u - v; }
        int c = 0;
        for (int a = 0; a < N; a++) if (w[a] != 0) c++;
        supp[b] = malloc(sizeof(int) * (c ? c : 1));
        supplen[b] = 0;
        for (int a = 0; a < N; a++) {
            if (w[a] == 0) zero[(size_t)b * N + a] = 1;
            else supp[b][supplen[b]++] = a;
        }
    }
    free(w);
}

/* current B */
static int m, Bbasis[MAXN], Belem[1 << MAXN];  /* Belem[t] = XOR of basis vectors selected by t */
static int Bidx;
/* current P */
static int Py[MAXN];
static long long found;

static inline int coset(int a) {
    int c = 0;
    for (int i = 0; i < m; i++) c |= (__builtin_popcount(a & Py[i]) & 1) << i;
    return c;
}

static uint64_t allowed[1 << MAXN];   /* allowed[t] = cosets (bitmask over 2^m ids) for b = Belem[t] */
static int order[1 << MAXN];          /* test order of b's; move-to-front on rejection */
static int phi[1 << MAXN];

static void emit(void) {
    /* basis of U: (rep(phi(b_i)), b_i) for i < m, and (a, 0) for a basis of P^perp */
    int rep[64];
    for (int c = 0; c < (1 << m); c++) rep[c] = -1;
    int perp[MAXN], np = 0, piv[MAXN];
    for (int a = 0; a < N; a++) {
        int c = coset(a);
        if (rep[c] < 0) rep[c] = a;
        if (c == 0 && a && np < n - m) {
            int v = a;  /* insert into echelon basis of P^perp */
            for (int i = 0; i < np; i++) if ((v >> piv[i]) & 1) v ^= perp[i];
            if (v) { piv[np] = 31 - __builtin_clz(v); perp[np++] = v;
                     for (int i = 0; i < np - 1; i++) if ((perp[i] >> piv[np - 1]) & 1) perp[i] ^= v; }
        }
    }
    printf("U %d", Bidx);
    for (int i = 0; i < m; i++) printf(" %lld", (long long)rep[phi[1 << i]] + ((long long)Bbasis[i] << n));
    for (int i = 0; i < np; i++) printf(" %d", perp[i]);
    printf("\n");
    found++;
}

/* choose phi on the basis of B one vector at a time; check every b in the span so far */
static void bt(int A) {
    if (A == (1 << m) - 1) { emit(); return; }
    int j = __builtin_ctz(~A);
    int bj = 1 << j;
    for (int v = 0; v < (1 << m); v++) {
        int ok = 1;
        for (int beta = A;; beta = (beta - 1) & A) {
            int val = v ^ phi[beta];
            if (!((allowed[bj | beta] >> val) & 1)) { ok = 0; break; }
            if (!beta) break;
        }
        if (!ok) continue;
        for (int beta = A;; beta = (beta - 1) & A) { phi[bj | beta] = v ^ phi[beta]; if (!beta) break; }
        bt(A | bj);
    }
}

static long long nP;
static int shard, nshards;

static void test_P(void) {
    if (nP++ % nshards != shard) return;
    uint64_t full = (m == 6) ? ~0ULL : ((1ULL << (1 << m)) - 1);
    int T = 1 << m;
    for (int k = 0; k < T - 1; k++) {
        int t = order[k], b = Belem[t];
        uint64_t hit = 0;
        for (int i = 0; i < supplen[b] && hit != full; i++) hit |= 1ULL << coset(supp[b][i]);
        allowed[t] = full & ~hit;
        if (!allowed[t]) {  /* move to front: this b is likely to reject the next P too */
            memmove(order + 1, order, sizeof(int) * k);
            order[0] = t;
            return;
        }
    }
    allowed[0] = 1;  /* phi(0) = 0 */
    phi[0] = 0;
    bt(0);
}

/* enumerate m-dim subspaces of F_2^n by reduced row echelon form: pivots p_0 > p_1 > ... (leading bits) */
static int piv[MAXN];
static void rows_rec(int r) {
    if (r == m) { test_P(); return; }
    /* free positions of row r: bits below its pivot that are not pivots of other rows */
    int freemask = 0;
    for (int bit = 0; bit < piv[r]; bit++) {
        int isp = 0;
        for (int s = 0; s < m; s++) if (piv[s] == bit) isp = 1;
        if (!isp) freemask |= 1 << bit;
    }
    int nf = __builtin_popcount(freemask);
    for (int x = 0; x < (1 << nf); x++) {
        int v = 1 << piv[r], k = 0;
        for (int bit = 0; bit < piv[r]; bit++) if ((freemask >> bit) & 1) { if ((x >> k) & 1) v |= 1 << bit; k++; }
        Py[r] = v;
        rows_rec(r + 1);
    }
}
static void pivots_rec(int r, int maxbit) {
    if (r == m) { rows_rec(0); return; }
    for (int p = maxbit; p >= m - r - 1; p--) { piv[r] = p; pivots_rec(r + 1, p - 1); }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: wzs N LUTFILE BFILE [SHARD NSHARDS]\n"); return 1; }
    n = atoi(argv[1]); N = 1 << n;
    shard = argc > 5 ? atoi(argv[4]) : 0;
    nshards = argc > 5 ? atoi(argv[5]) : 1;
    F = malloc(sizeof(int) * N);
    FILE *f = fopen(argv[2], "r");
    for (int x = 0; x < N; x++) if (fscanf(f, " %d ,", &F[x]) != 1) { fprintf(stderr, "bad LUT\n"); return 1; }
    fclose(f);
    walsh();
    FILE *bf = fopen(argv[3], "r");
    char line[4096];
    Bidx = 0;
    while (fgets(line, sizeof line, bf)) {
        char *p = line;
        m = 0;
        while (*p) {
            while (*p == ' ' || *p == ',') p++;
            if (*p < '0' || *p > '9') break;
            Bbasis[m++] = (int)strtol(p, &p, 10);
        }
        if (!m) continue;
        for (int t = 0; t < (1 << m); t++) {
            int v = 0;
            for (int i = 0; i < m; i++) if ((t >> i) & 1) v ^= Bbasis[i];
            Belem[t] = v;
        }
        for (int k = 0; k < (1 << m) - 1; k++) order[k] = k + 1;
        nP = 0; found = 0;
        pivots_rec(0, n - 1);
        fprintf(stderr, "B %d (dim %d): %lld subspaces P enumerated (all shards), %lld Walsh zero spaces in this shard\n",
                Bidx, m, nP, found);
        fflush(stdout);
        Bidx++;
    }
    return 0;
}

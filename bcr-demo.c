#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct {
	uint64_t u, v;
} pair64_t;

#define rstype_t pair64_t
#define rskey(x) ((x).u)

#define RS_MIN_SIZE 64

typedef struct {
	rstype_t *b, *e;
} rsbucket_t;

void rs_sort(rstype_t *beg, rstype_t *end, int n_bits, int s) // integer radix sort
{
	rstype_t *i;
	int size = 1<<n_bits, m = size - 1;
	rsbucket_t *k, b[size], *be = b + size;

	for (k = b; k != be; ++k) k->b = k->e = beg;
	for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; // count radix
	for (k = b + 1; k != be; ++k) // set start and end of each bucket
		k->e += (k-1)->e - beg, k->b = (k-1)->e;
	for (k = b; k != be;) { // in-place classification based on radix
		if (k->b != k->e) { // the bucket is not full
			rsbucket_t *l;
			if ((l = b + (rskey(*k->b)>>s&m)) != k) { // destination different
				rstype_t tmp = *k->b, swap;
				do { // swap until we find an element in $k
					swap = tmp; tmp = *l->b; *l->b++ = swap;
					l = b + (rskey(tmp)>>s&m);
				} while (l != k);
				*k->b++ = tmp;
			} else ++k->b;
		} else ++k;
	}
	for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; // reset k->b
	if (s) { // if $s is non-zero, we need to sort buckets
		s = s > n_bits? s - n_bits : 0;
		for (k = b; k != be; ++k)
			if (k->e - k->b > RS_MIN_SIZE) rs_sort(k->b, k->e, n_bits, s);
			else if (k->e - k->b > 1) // then use an insertion sort
				for (i = k->b + 1; i < k->e; ++i)
					if (rskey(*i) < rskey(*(i - 1))) {
						rstype_t *j, tmp = *i;
						for (j = i; j > k->b && rskey(tmp) < rskey(*(j-1)); --j)
							*j = *(j - 1);
						*j = tmp;
					}
	}
}

/**
 * Append $T to existing BWT $B.
 *
 * @param Blen    length of the existing BWT
 * @param B       existing BWT; set to NULL if non-existing
 * @param Tlen    length of input string
 * @param T       input string; '\0' represents a sentinel
 *
 * @return  the new BWT string
 */
uint8_t *bcr_lite(long Blen, uint8_t *B, long Tlen, const uint8_t *T)
{
	long i, k, n, max, n0;
	uint8_t *p, *q, *B0;
	const uint8_t *end, **P;
	pair64_t *a;
	int c;
	// split $T into short strings at sentinels
	if (T == 0 || Tlen == 0) return B;
	for (p = q = (uint8_t*)T, end = T + Tlen, n = max = 0; p != end; ++p) {
		if (*p) continue;
		if (n == max) {
			max = max? max<<1 : 256;
			P = realloc(P, max * sizeof(void*));
		}
		P[n++] = q, q = p + 1;
	}
	P = realloc(P, (n + 1) * sizeof(void*));
	P[n] = q;
	// initialize
	for (p = B, end = B + Blen, i = 0; p < end; ++p) i += (*p == 0); // count # of sentinels
	a = malloc(sizeof(pair64_t) * n);
	for (k = 0; k < n; ++k) a[k].u = k + i, a[k].v = k<<8;
	B = realloc(B, Blen + Tlen);
	memmove(B + Tlen, B, Blen); // finished BWT is always placed at the end of $B
	B = B0 = B + Tlen;
	// core loop
	for (i = 0, n0 = n; n0; ++i) {
		long l, pre, ac[256], mc[256];
		for (c = 0; c != 256; ++c) mc[c] = 0;
		end = B0 + Blen; Blen += n0; B -= n0;
		for (n = k = 0, p = B0, q = B, pre = -1L; k < n0; ++k) {
			pair64_t *u = &a[k];
			c = P[(u->v>>8) + 1] - 2 - i >= P[u->v>>8]? *(P[(u->v>>8) + 1] - 2 - i) : 0; // symbol to insert
			u->v = (u->v&~0xffULL) | c;
			if (u->u - pre - 1) // then copy ($u->u - $pre - 1) symbols from B0 to B
				for (l = 0; l < u->u - pre - 1; ++l)
					++mc[*p], *q++ = *p++;
			*q++ = c;
			pre = u->u; u->u = mc[c]++;
			if (c) a[n++] = a[k];
		}
		while (p < end) ++mc[*p], *q++ = *p++; // copy the rest of $B0 to $B
		for (c = 1, ac[0] = 0; c != 256; ++c) ac[c] = ac[c-1] + mc[c-1]; // accumulative count
		for (k = 0; k < n; ++k) a[k].u += ac[a[k].v&0xff] + n; // compute the absolute positions
		for (k = Blen, c = 0; k; k >>= 1, ++c); // find how many bits needed for radix sort
		if (n) rs_sort(a, a + n, 8, c > 7? c - 7 : 0); // radix sort
		B0 = B; n0 = n;
	}
	free(P); free(a);
	return B;
}

#include <stdio.h>

int main(int argc, char *argv[])
{
	uint8_t *B, *s;
	long i, len;
	FILE *fp;
	if (argc == 1) {
		fprintf(stderr, "Usage: bcr-demo <in.txt>\n");
		return 1;
	}
	fp = fopen(argv[1], "rb"); // FIXME: check if $fp==0
	if (fseek(fp, 0, SEEK_END)) { // FIXME: can be adapted to work on a stream
		fprintf(stderr, "ERROR: fail to determine the file size; make sure the input is an ordinary file, not a stream.\n");
		return 2;
	}
	len = ftell(fp); // file size
	fseek(fp, 0, SEEK_SET);
	s = malloc(len);
	fread(s, 1, len, fp);
	fclose(fp);
	for (i = 0; i < len; ++i)
		if (s[i] == '\n') s[i] = 0;
	B = bcr_lite(0, 0, len, s);
	for (i = 0; i < len; ++i)
		putchar(B[i]? B[i] : '$');
	putchar('\n');
	free(s); free(B);
	return 0;
}

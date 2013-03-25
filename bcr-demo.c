#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct {
	uint64_t u, v;
} pair64_t;

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
	const uint8_t *end, **P = 0;
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
		long l, pre, ac[256], mc[256], mc2[256];
		pair64_t *b[256], *aa;
		for (c = 0; c != 256; ++c) mc[c] = mc2[c] = 0;
		end = B0 + Blen; Blen += n0; B -= n0;
		for (n = k = 0, p = B0, q = B, pre = 0; k < n0; ++k) {
			pair64_t *u = &a[k];
			c = P[(u->v>>8) + 1] - 2 - i >= P[u->v>>8]? *(P[(u->v>>8) + 1] - 2 - i) : 0; // symbol to insert
			u->v = (u->v&~0xffULL) | c;
			for (l = 0; l != u->u - pre; ++l) // copy ($u->u - $pre - 1) symbols from B0 to B
				++mc[*p], *q++ = *p++; // $mc: marginal counts of all processed symbols
			*q++ = c;
			pre = u->u + 1; u->u = mc[c]++;
			if (c) a[n++] = a[k], ++mc2[c]; // $mc2: marginal counts of the current column
		}
		while (p < end) ++mc[*p], *q++ = *p++; // copy the rest of $B0 to $B
		for (c = 1, ac[0] = 0; c != 256; ++c) ac[c] = ac[c-1] + mc[c-1]; // accumulative count
		for (k = 0; k < n; ++k) a[k].u += ac[a[k].v&0xff] + n; // compute positions for the next round
		// stable counting sort ($a[k].v&0xff); also possible with an in-place non-stable radix sort, which is slower
		aa = malloc(sizeof(pair64_t) * n);
		for (c = 1, b[0] = aa; c != 256; ++c) b[c] = b[c-1] + mc2[c-1];
		for (k = 0; k < n; ++k) *b[a[k].v&0xff]++ = a[k]; // this works because $a is already partially sorted
		free(a); a = aa; // $aa now becomes $a
		B0 = B; n0 = n;
	}
	free(P); free(a);
	return B;
}

#include "ksort.h"
#define pair64_lt(a, b) (((a).u<<8|((a).v&0xff)) < ((b).u<<8|((b).v&0xff)))
KSORT_INIT(tmp, pair64_t, pair64_lt)

uint8_t *bcr_rlo(long Tlen, const uint8_t *T)
{
	long i, k, n, max, n0, Blen = 0;
	uint8_t *p, *q, *B0, *B = 0;
	const uint8_t *end, **P = 0;
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
	a = malloc(sizeof(pair64_t) * n);
	for (k = 0; k < n; ++k) a[k].u = 0, a[k].v = k<<8;
	B = realloc(B, Blen + Tlen);
	memmove(B + Tlen, B, Blen); // finished BWT is always placed at the end of $B
	B = B0 = B + Tlen;
	// core loop
	for (i = 0, n0 = n; n0; ++i) {
		long l, pre, ac[256], mc[256], mc2[256];
		pair64_t *b[256], *aa;
		for (c = 0; c != 256; ++c) mc[c] = mc2[c] = 0;
		end = B0 + Blen; Blen += n0; B -= n0;
		for (k = 0; k < n0; ++k) { // set the base to insert
			pair64_t *u = &a[k];
			u->v = (u->v&~0xffULL) | (P[(u->v>>8) + 1] - 2 - i >= P[u->v>>8]? *(P[(u->v>>8) + 1] - 2 - i) : 0);
		}
		ks_introsort(tmp, n0, a); // we can use a partial counting sort, which is faster than this introsort
		for (k = 1, l = 0, pre = a[0].u; k < n0; ++k) { // differentiate identical ranks (*)
			if (a[k].u == pre) a[k].u += k - l;
			else l = k, pre = a[k].u;
		}

		for (n = k = 0, p = B0, q = B, pre = 0; k < n0; ++k) {
			pair64_t *u = &a[k];
			c = u->v & 0xff;
			for (l = 0; l != u->u - pre; ++l) // copy symbols from B0 to B
				++mc[*p], *q++ = *p++; // $mc: marginal counts of all processed symbols
			*q++ = c;
			pre = u->u + 1; u->u = mc[c]++;
			if (c) a[n++] = a[k], ++mc2[c]; // $mc2: marginal counts of the current column
		}
		while (p < end) ++mc[*p], *q++ = *p++; // copy the rest of $B0 to $B
		for (c = 1, ac[0] = 0; c != 256; ++c) ac[c] = ac[c-1] + mc[c-1]; // accumulative count
		for (k = 0; k < n; ++k) a[k].u += ac[a[k].v&0xff] + n; // compute positions for the next round
		// stable counting sort ($a[k].v&0xff); also possible with an in-place non-stable radix sort, which is slower
		aa = malloc(sizeof(pair64_t) * n);
		for (c = 1, b[0] = aa; c != 256; ++c) b[c] = b[c-1] + mc2[c-1];
		for (k = 0; k < n; ++k) *b[a[k].v&0xff]++ = a[k]; // this works because $a is already partially sorted
		free(a); a = aa; // $aa now becomes $a
		for (k = 1, pre = a[0].u<<8|(a[0].v&0xff), l = 0; k < n; ++k) { // largely the reverse of loop (*)
			if ((a[k].u<<8|(a[k].v&0xff)) == pre + ((k-l)<<8))
				a[k].u = a[l].u;
			else pre = a[k].u<<8|(a[k].v&0xff), l = k;
		}
		B0 = B; n0 = n;
	}
	free(P); free(a);
	return B;
}

#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
	uint8_t *B, *s;
	long i, len;
	FILE *fp;
	int c, is_rlo = 0;

	while ((c = getopt(argc, argv, "s")) >= 0)
		if (c == 's') is_rlo = 1; 

	if (argc == optind) {
		fprintf(stderr, "Usage: bcr-demo [-s] <in.txt>\n");
		return 1;
	}
	fp = fopen(argv[optind], "rb"); // FIXME: check if $fp==0
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
	B = is_rlo? bcr_rlo(len, s) : bcr_lite(0, 0, len, s);
	for (i = 0; i < len; ++i)
		putchar(B[i]? B[i] : '$');
	putchar('\n');
	free(s); free(B);
	return 0;
}
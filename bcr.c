#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "rld.h"
#include "ksort.h"

/*************************************************
 *** Data structure for long 2-bit encoded DNA ***
 *************************************************/

#define LD_SHIFT 20
#define LD_MASK  ((1U<<LD_SHIFT) - 1)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef struct {
	int max;
	uint64_t **a;
} longdna_t;

longdna_t *ld_init(void)
{
	return calloc(1, sizeof(longdna_t));
}

void ld_destroy(longdna_t *ld)
{
	int j;
	for (j = 0; j < ld->max; ++j) free(ld->a[j]);
	free(ld);
}

inline void ld_set(longdna_t *h, int64_t x, int c)
{
	int k = x >> LD_SHIFT, l = x & LD_MASK;
	if (k >= h->max) {
		int j, old_max = h->max;
		h->max = k + 1;
		kroundup32(h->max);
		h->a = realloc(h->a, sizeof(longdna_t) * h->max);
		for (j = old_max; j < h->max; ++j) h->a[j] = 0;
	}
	if (h->a[k] == 0) h->a[k] = calloc(1<<LD_SHIFT>>5, 8);
	h->a[k][l>>5] |= (uint64_t)(c&3)<<((l&31)<<1); // NB: we cannot set the same position multiple times
}

inline int ld_get(longdna_t *h, int64_t x)
{
	int l = x & LD_MASK;
	return h->a[x>>LD_SHIFT][l>>5]>>((l&31)<<1)&3;
}

/***********
 *** BCR ***
 ***********/

typedef struct {
	uint64_t u, v; // $u: position; $v: seq_id:61, base:3
} pair64_t;

#define bcr_lt(a, b) ((a).u < (b).u)
KSORT_INIT(bcr, pair64_t, bcr_lt)

typedef struct {
	rld_t *e;
	int64_t n, c[6];
	pair64_t *a;
} pbwt_t;

typedef struct {
	int max_len;
	uint64_t n_seqs, m_seqs, c[6], c0[6][6];
	uint8_t *len;
	longdna_t **seq;
	pair64_t *a;
	pbwt_t bwt[6];
} bcr_t;

bcr_t *bcr_init()
{
	bcr_t *b;
	int i;
	b = calloc(1, sizeof(bcr_t));
	for (i = 0; i < 6; ++i) b->bwt[i].e = rld_init(6, 6);
	return b;
}

void bcr_destroy(bcr_t *b)
{
	free(b->len); free(b->a); free(b->seq);
	free(b);
}

void bcr_append(bcr_t *b, int len, uint8_t *seq)
{
	int i;
	assert(len < 256 && len > 1);
	if (len > b->max_len) { // find a longer read
		b->seq = realloc(b->seq, len * sizeof(void*));
		for (i = b->max_len; i < len; ++i)
			b->seq[i] = ld_init();
		b->max_len = len;
	}
	if (b->n_seqs == b->m_seqs) {
		b->m_seqs = b->m_seqs? b->m_seqs<<1 : 256;
		b->len = realloc(b->len, b->m_seqs);
	}
	b->len[b->n_seqs] = len;
	for (i = 0; i < len; ++i)
		ld_set(b->seq[i], b->n_seqs, seq[len - 1 - i] - 1);
	++b->n_seqs;
}

static void print_bwt(rld_t *e, int endl)
{
	int64_t l, i;
	int c;
	rlditr_t itr;
	if (e->mcnt[0]) {
		rld_itr_init(e, &itr, 0);
		while ((l = rld_dec(e, &itr, &c, 0)) != -1)
			for (i = 0; i < l; ++i)
				fputc("$ACGTN"[c], stderr);
	}
	if (endl) fputc(endl, stderr);
}

static void set_bwt(bcr_t *bcr)
{
	int64_t k, c[6], i[6];
	int j, l;
	pair64_t *a;
	memset(c, 0, 48);
	for (k = 0; k < bcr->n_seqs; ++k) {
		pair64_t *u = &bcr->a[k];
		u->u += c[u->v&7];
		++c[u->v&7];
	}
	for (j = 0; j < 6; ++j) bcr->bwt[j].n = c[j];
	for (l = 0; l < 6; ++l) bcr->bwt[0].c[l] = 0;
	for (j = 1; j < 6; ++j)
		for (l = 0; l < 6; ++l)
			bcr->bwt[j].c[l] = bcr->bwt[j-1].e->mcnt[l+1];
	for (j = 1; j < 6; ++j)
		for (l = 0; l < 6; ++l)
			bcr->bwt[j].c[l] += bcr->bwt[j-1].c[l];
	memmove(c + 1, c, 40);
	for (k = 1, c[0] = 0; k < 6; ++k) c[k] += c[k - 1];
	a = malloc(bcr->n_seqs * 16);
	for (k = 0; k < 6; ++k) i[k] = c[k], bcr->c[k] += c[k], bcr->bwt[k].a = a + c[k];
	for (k = 0; k < bcr->n_seqs; ++k) {
		pair64_t *u = &bcr->a[k];
		int b = bcr->a[k].v&7;
		bcr->a[k].u += c[b];
		a[i[b]++] = bcr->a[k];
		fprintf(stderr, "[1] k=%lld, u=%lld, i=%lld, c=%c\n", k, u->u, u->v>>3, "$ACGTN"[u->v&7]);
	}
	free(bcr->a);
	bcr->a = a;
}

static void next_bwt(bcr_t *bcr, int class, int pos)
{
	int64_t c[6], k, l;
	rlditr_t ir, iw;
	pbwt_t *bwt = &bcr->bwt[class];
	rld_t *ew, *er = bwt->e;

	if (bwt->n == 0) return;
	if (class) ks_introsort(bcr, bwt->n, bwt->a);
	for (l = 0; l < 6; ++l) fprintf(stderr, "%lld, ", bwt->c[l]); fputc('\n', stderr);
	for (k = 0; k < bwt->n; ++k) {
		pair64_t *u = &bwt->a[k];
		u->u -= k + bcr->c[class];
		u->v = (u->v&~7ULL) | (pos == bcr->max_len? 0 : ld_get(bcr->seq[pos], u->v>>3) + 1);
		fprintf(stderr, "[2] class=%c, pos=%d, k=%lld, u=%lld, i=%lld, c=%c\n", "$ACGTN"[class], pos, k, u->u, u->v>>3, "$ACGTN"[u->v&7]);
	}
	ew = rld_init(6, 6);
	rld_itr_init(er, &ir, 0);
	rld_itr_init(ew, &iw, 0);
	memset(c, 0, 48);
	for (k = l = 0; k < bwt->n; ++k) {
		pair64_t *u = &bwt->a[k];
		if (u->u > l)
			l += rld_dec_enc_cnt(ew, &iw, er, &ir, u->u - l, c);
		rld_enc(ew, &iw, 1, u->v&7);
		u->u = c[u->v&7] + bcr->c[u->v&7] + bwt->c[u->v&7];
	}
	if (l < er->mcnt[0])
		l += rld_dec_enc_cnt(ew, &iw, er, &ir, er->mcnt[0] - l, c);
	rld_enc_finish_partial(ew, &iw);
//	rld_enc_finish(ew, &iw);
	rld_destroy(er);
	bwt->e = ew;
	for (k = 0; k < bwt->n; ++k) {
		pair64_t *u = &bwt->a[k];
		fprintf(stderr, "[3] class=%c, pos=%d, k=%lld, u=%lld, i=%lld, c=%c, bcr->c=%lld, bwt->c=%lld\n", "$ACGTN"[class], pos, k, u->u, u->v>>3, "$ACGTN"[u->v&7], bcr->c[u->v&7], bwt->c[u->v&7]);
	}
}

rld_t *bcr_build(bcr_t *b)
{
	int64_t k;
	int pos, c;
	rld_t *e;
	rlditr_t iw;

	b->m_seqs = b->n_seqs;
	b->len = realloc(b->len, b->n_seqs);
	b->a = malloc(b->n_seqs * 16);
	assert(b->a);
	for (k = 0; k < b->n_seqs; ++k) b->a[k].u = 0, b->a[k].v = k<<3;
	for (pos = 0; pos <= b->max_len; ++pos) {
		set_bwt(b);
		if (pos) {
			for (c = 1; c <= 4; ++c)
				next_bwt(b, c, pos);
		} else next_bwt(b, 0, pos);
		for (c = 0; c < 5; ++c) print_bwt(b->bwt[c].e, ','); fputc('\n', stderr);
	}

	e = rld_init(6, 3);
	rld_itr_init(e, &iw, 0);
	for (c = 0; c <= 4; ++c) {
		rld_t *er = b->bwt[c].e;
		rlditr_t ir;
		if (er->mcnt[0] == 0) continue;
		rld_itr_init(er, &ir, 0);
		rld_dec_enc(e, &iw, er, &ir, er->mcnt[0]);
		rld_destroy(er);
		b->bwt[c].e = 0;
	}
	rld_enc_finish(e, &iw);
	return e;
}

/*********************
 *** Main function ***
 *********************/

#include <zlib.h>
#include <unistd.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void seq_char2nt6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i) {
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
		if (s[i] >= 5) s[i] = (lrand48()&3) + 1;
	}
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	int c, for_only = 0, n_threads = 1;
	bcr_t *bcr;
	rld_t *e;

	while ((c = getopt(argc, argv, "ft:")) >= 0)
		if (c == 'f') for_only = 1;
		else if (c == 't') n_threads = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: bcr-mt [-f] [-t nThreads=1] <in.fq.gz>\n");
		return 1;
	}

	bcr = bcr_init();
	fp = gzopen(argv[optind], "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		uint8_t *s = (uint8_t*)ks->seq.s;
		seq_char2nt6(ks->seq.l, s);
		bcr_append(bcr, ks->seq.l, s);
		if (!for_only) {
			seq_revcomp6(ks->seq.l, s);
			bcr_append(bcr, ks->seq.l, s);
		}
	}
	kseq_destroy(ks);
	gzclose(fp);

	e = bcr_build(bcr);
	rld_dump(e, "-");
	bcr_destroy(bcr);
	return 0;
}

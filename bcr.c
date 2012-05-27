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

#define LD_SHIFT 25
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
	h->a[k][l>>5] |= (c&3)<<(l&31); // NB: we cannot set the same position multiple times
}

inline int ld_get(longdna_t *h, int64_t x)
{
	int k = x >> LD_SHIFT, l = x & LD_MASK;
	return h->a[k][l>>5]>>(l&31)&3;
}

/***********
 *** BCR ***
 ***********/

typedef struct {
	uint64_t u, v; // $u: position; $v: seq_id:61, base:3
} pair64_t;

#define bcr_lt(a, b) ((a).u < (b).u || ((a).u == (b).u && (a).v < (b).v))
KSORT_INIT(bcr, pair64_t, bcr_lt)

typedef struct {
	int max_len, n_threads;
	uint64_t n_seqs, m_seqs;
	uint8_t *len;
	longdna_t **seq;
	pair64_t *a;
	rlditr_t itr;
	rld_t *e;
} bcr_t;

bcr_t *bcr_init(int n_threads)
{
	bcr_t *b;
	b = calloc(1, sizeof(bcr_t));
	b->n_threads = n_threads > 0? n_threads : 1;
	b->e = rld_init(6, 3);
	rld_itr_init(b->e, &b->itr, 0);
	return b;
}

void bcr_destroy(bcr_t *b)
{
	free(b->len); free(b->a); free(b->seq);
	free(b);
}

void bcr_append(bcr_t *b, int len, uint8_t *seq)
{
	int i, c;
	assert(len < 256 && len > 1);
	if (len > b->max_len) { // find a longer read
		b->seq = realloc(b->seq, (len - 1) * sizeof(void*));
		for (i = b->max_len; i < len - 1; ++i)
			b->seq[i] = ld_init();
		b->max_len = len;
	}
	if (b->n_seqs == b->m_seqs) {
		b->m_seqs = b->m_seqs? b->m_seqs<<1 : 256;
		b->len = realloc(b->len, 1);
	}
	b->len[b->n_seqs] = len;
	c = (seq[len - 1] & 3) + 1;
	rld_enc(b->e, &b->itr, 1, c);
	for (i = 0; i < len - 2; ++i)
		ld_set(b->seq[i], b->n_seqs, seq[len - 2 - i]&3);
	++b->n_seqs;
}

void bcr_prepare(bcr_t *b)
{
	int64_t k;
	b->m_seqs = b->n_seqs;
	b->len = realloc(b->len, b->n_seqs);
	rld_enc_finish(b->e, &b->itr);
	b->a = malloc(b->n_seqs * 16);
	assert(b->a);
	for (k = 0; k < b->n_seqs; ++k)
		b->a[k].u = k, b->a[k].v = k<<3;
}

typedef struct {
	bcr_t *b;
	int start, step, i;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int64_t i, ok[6];
	for (i = w->start; i < w->b->n_seqs; i += w->step) {
		pair64_t *p = &w->b->a[i];
		int c = p->v&7;
		if (w->i < w->b->max_len - 1 && w->i > w->b->len[i] - 1) continue; // FIXME: check if this works for variable-length strings
		rld_rank1a(w->b->e, p->u, (uint64_t*)ok);
		p->u = w->b->e->cnt[c] + ok[c] - 1;
		p->v = p->v>>3<<3 | (w->i < w->b->max_len - 1? ld_get(w->b->seq[w->i], p->v>>3) : 0);
	}
	return data;
}

void bcr_build1(bcr_t *b, int which)
{
	rld_t *e0;
	rlditr_t itr0;
	int64_t i, last;
	int j;
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;

	// dispatch workers
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(b->n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(b->n_threads, sizeof(pthread_t));
	for (j = 0; j < b->n_threads; ++j)
		w[j].b = b, w[j].start = j, w[j].step = b->n_threads, w[j].i = which;
	for (j = 0; j < b->n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < b->n_threads; ++j) pthread_join(tid[j], 0);
	free(w); free(tid);
	if (which < b->max_len - 1) ld_destroy(b->seq[which]);

	// insert to the current BWT; similar to fm_merge_from_SA() in fermi
	e0 = b->e;
	b->e = rld_init(6, 3);
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(b->e, &b->itr, 0);
	ks_introsort(bcr, b->n_seqs, b->a);
	for (i = last = 0; i < b->n_seqs; ++i) {
		pair64_t *p = &b->a[i];
		if (p->u != last) {
			rld_dec_enc(b->e, &b->itr, e0, &itr0, p->u - last);
			last = p->u;
		}
		rld_enc(b->e, &b->itr, 1, p->v&7);
		p->u += i;
	}
	if (last != e0->mcnt[0] - 1)
		rld_dec_enc(b->e, &b->itr, e0, &itr0, e0->mcnt[0] - 1 - last);
	rld_destroy(e0);
	rld_enc_finish(b->e, &b->itr);
}

void bcr_build(bcr_t *b)
{
	int j;
	for (j = 0; j < b->max_len - 1; ++j) bcr_build1(b, j);
	bcr_build1(b, -1);
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
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
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

	while ((c = getopt(argc, argv, "ft:")) >= 0)
		if (c == 'f') for_only = 1;
		else if (c == 't') n_threads = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: bcr-mt [-f] [-t nThreads=1] <in.fq.gz>\n");
		return 1;
	}

	bcr = bcr_init(n_threads);
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

	bcr_build(bcr);
	rld_dump(bcr->e, "-");
	bcr_destroy(bcr);
	return 0;
}

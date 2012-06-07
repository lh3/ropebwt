#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <assert.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/**********************************************
 *** Lightweight run-length encoder/decoder ***
 **********************************************/

#define RLL_BLOCK_SIZE 0x100000

typedef struct {
	int c;
	int64_t l;
	uint8_t *q, **i;
} rllitr_t;

typedef struct {
	int n, m;
	uint8_t **z;
	int64_t l, mc[6];
} rll_t;

static rll_t *rll_init(void)
{
	rll_t *e;
	e = calloc(1, sizeof(rll_t));
	e->n = e->m = 1;
	e->z = malloc(sizeof(void*));
	e->z[0] = calloc(RLL_BLOCK_SIZE, 1);
	e->z[0][0] = 7;
	return e;
}

static void rll_destroy(rll_t *e)
{
	int i;
	if (e == 0) return;
	for (i = 0; i < e->n; ++i) free(e->z[i]);
	free(e->z); free(e);
}

static void rll_itr_init(const rll_t *e, rllitr_t *itr)
{
	itr->i = e->z; itr->q = *itr->i; itr->c = -1; itr->l = 0;
}

static inline void rll_enc0(rll_t *e, rllitr_t *itr, int l, uint8_t c)
{
	*itr->q++ = l<<3 | c;
	e->mc[c] += l;
	if (itr->q - *itr->i == RLL_BLOCK_SIZE) {
		if (e->n == e->m) {
			e->m <<= 1;
			e->z = realloc(e->z, e->m * sizeof(void*));
			memset(e->z + e->n, 0, (e->m - e->n) * sizeof(void*));
		}
		++e->n;
		itr->i = e->z + e->n - 1;
		itr->q = *itr->i = calloc(RLL_BLOCK_SIZE, 1);
	}
}

static inline void rll_enc(rll_t *e, rllitr_t *itr, int64_t l, uint8_t c)
{
	if (itr->c != c) {
		if (itr->l) {
			if (itr->l > 31) 
				for (; itr->l > 31; itr->l -= 31)
					rll_enc0(e, itr, 31, itr->c);
			rll_enc0(e, itr, itr->l, itr->c);
		}
		itr->l = l; itr->c = c;
	} else itr->l += l;
}

static void rll_enc_finalize(rll_t *e, rllitr_t *itr)
{
	int c;
	rll_enc(e, itr, 0, -1);
	*itr->q = 7; // end marker; there is always room for an extra symbol
	for (e->l = 0, c = 0; c < 6; ++c) e->l += e->mc[c];
}

static inline int64_t rll_dec(const rll_t *e, rllitr_t *itr, int *c, int is_free)
{
	int64_t l;
	if (*itr->q == 7) return -1;
	l = *itr->q>>3; *c = *itr->q&7;
	if (++itr->q - *itr->i == RLL_BLOCK_SIZE) {
		if (is_free) {
			free(*itr->i);
			*itr->i = 0;
		}
		itr->q = *++itr->i;
	}
	return l;
}

static inline void rll_copy(rll_t *e, rllitr_t *itr, const rll_t *e0, rllitr_t *itr0, int64_t k)
{
	if (itr0->l >= k) { // there are more pending symbols
		rll_enc(e, itr, k, itr0->c);
		itr0->l -= k; // l - k symbols remains
	} else { // use up all pending symbols
		int c = -1; // to please gcc
		int64_t l;
		rll_enc(e, itr, itr0->l, itr0->c); // write all pending symbols
		k -= itr0->l;
		for (; k > 0; k -= l) { // we always go into this loop because l0<k
			l = rll_dec(e0, itr0, &c, 1);
			rll_enc(e, itr, k < l? k : l, c);
		}
		itr0->l = -k; itr0->c = c;
	}
}

/*************************************************
 *** Data structure for long 2-bit encoded DNA ***
 *************************************************/

#define LD_SHIFT 20
#define LD_MASK  ((1U<<LD_SHIFT) - 1)

typedef struct {
	int max;
	uint64_t **a;
} longdna_t; // to allocate, simply call calloc()

void ld_destroy(longdna_t *ld)
{
	int j;
	for (j = 0; j < ld->max; ++j) free(ld->a[j]);
	free(ld->a); free(ld);
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
	return h->a[x>>LD_SHIFT][(x&LD_MASK)>>5]>>((x&31)<<1)&3;
}

/******************
 *** Radix sort ***
 ******************/

typedef struct {
	uint64_t u, v; // $u: position; $v: seq_id:61, base:3
} pair64_t;

#define rsint_t uint64_t
#define rstype_t pair64_t
#define rskey(x) ((x).u)

#define RS_MIN_SIZE 255

typedef struct {
	rstype_t *b, *e;
} rsbucket_t;

static inline void rs_insertsort(rstype_t *s, rstype_t *t)
{
	rstype_t *i;
	for (i = s + 1; i < t; ++i) {
		rstype_t *j, tmp = *i;
		for (j = i; j > s && rskey(tmp) < rskey(*(j-1)); --j)
			*j = *(j - 1);
		*j = tmp;
	}
}

void rs_combsort(size_t n, rstype_t a[])
{
	const double shrink_factor = 1. / 1.2473309501039786540366528676643;
	int do_swap;
	size_t gap = n;
	rstype_t tmp, *i, *j;
	do {
		if (gap > 2) {
			gap = (size_t)(gap * shrink_factor);
			if (gap == 9 || gap == 10) gap = 11;
		}
		do_swap = 0;
		for (i = a; i < a + n - gap; ++i) {
			j = i + gap;
			if (rskey(*j) < rskey(*i)) {
				tmp = *i; *i = *j; *j = tmp;
				do_swap = 1;
			}
		}
	} while (do_swap || gap > 2);
	if (gap != 1) rs_insertsort(a, a + n);
}

void rs_classify(rstype_t *beg, rstype_t *end, int n_bits, int s, rsbucket_t *b)
{
	rstype_t *i;
	int m = (1<<n_bits) - 1;
	rsbucket_t *k, *l, *be;

	be = b + (1<<n_bits);
	for (k = b; k != be; ++k) k->b = k->e = beg;
	for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e;
	if (b[0].e == end) return; // no need to sort
	for (k = b + 1; k != be; ++k)
		k->e += (k-1)->e - beg, k->b = (k-1)->e;
	for (k = b; k != be;) {
		rstype_t tmp, *p[2];
		int curr = 0;
		if (k->b == k->e) { ++k; continue; }
		l = b + (rskey(*k->b)>>s&m);
		if (k == l) { ++k->b; continue; }
		p[curr] = k->b; p[!curr] = &tmp;
		do {
			size_t tmp2 = l - b;
			while ((rskey(*l->b)>>s&m) == tmp2) ++l->b;
			*p[!curr] = *l->b;
			*l->b = *p[curr];
			curr = !curr;
			l = b + (rskey(*p[curr])>>s&m);
		} while (l != k);
		if (curr) *k->b = tmp;
		++k->b;
	}
	for (k = b + 1; k != be; ++k) k->b = (k-1)->e;
	b->b = beg;
}

void rs_sort(rstype_t *beg, rstype_t *end, int n_bits, int s)
{
	if (end - beg < 2) { // already sorted
		return;
	} else if (end - beg > RS_MIN_SIZE) {
		rsbucket_t *b;
		int i;
		b = malloc(sizeof(rsbucket_t) * (1<<n_bits));
		rs_classify(beg, end, n_bits, s, b);
		if (s) {
			s = s > n_bits? s - n_bits : 0;
			for (i = 0; i != 1<<n_bits; ++i)
				if (b[i].e > b[i].b + 1) rs_sort(b[i].b, b[i].e, n_bits, s);
		}
		free(b);
	} else rs_combsort(end - beg, beg);
}

/******************************
 *** Classify pair64_t::v&7 ***
 ******************************/

void rs_classify_alt(rstype_t *beg, rstype_t *end, int64_t *ac)
{
	rsbucket_t *b, *k, *l, *be;
	b = alloca(sizeof(rsbucket_t) * 8);
	be = b + 8;
	for (k = b; k != be; ++k) k->b = beg + ac[k-b];
	for (k = b; k != be - 1; ++k) k->e = (k+1)->b;
	k->e = end;
	for (k = b; k != be;) {
		rstype_t tmp, *p[2];
		int curr = 0;
		if (k->b == k->e) { ++k; continue; }
		l = b + ((*k->b).v&7);
		if (k == l) { ++k->b; continue; }
		p[curr] = k->b; p[!curr] = &tmp;
		do {
			size_t tmp2 = l - b;
			while (((*l->b).v&7) == tmp2) ++l->b;
			*p[!curr] = *l->b;
			*l->b = *p[curr];
			curr = !curr;
			l = b + (p[curr]->v&7);
		} while (l != k);
		if (curr) *k->b = tmp;
		++k->b;
	}
	for (k = b + 1; k != be; ++k) k->b = (k-1)->e;
	b->b = beg;
}

/***********
 *** BCR ***
 ***********/

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "bcr.h"

typedef struct {
	rll_t *e;
	int64_t n, c[6];
	pair64_t *a;
} bucket_t;

typedef struct {
	struct bcr_s *bcr;
	int class, pos, toproc;
} worker_t;

struct bcr_s {
	int max_len, n_threads;
	uint64_t n_seqs, m_seqs, c[6], tot;
	uint16_t *len;
	longdna_t **seq;
	bucket_t bwt[6];
#ifdef HAVE_PTHREAD
	pthread_t *tid;
	worker_t *w;
	pthread_mutex_t lock;
	pthread_cond_t cv;
	volatile int proc_cnt;
#endif
};

void bcr_append(bcr_t *b, int len, uint8_t *seq)
{
	int i;
	assert(len >= 1 && len < 65536);
	if (len > b->max_len) { // find a longer read
		b->seq = realloc(b->seq, len * sizeof(void*));
		for (i = b->max_len; i < len; ++i)
			b->seq[i] = calloc(1, sizeof(longdna_t));
		b->max_len = len;
	}
	if (b->n_seqs == b->m_seqs) {
		b->m_seqs = b->m_seqs? b->m_seqs<<1 : 256;
		b->len = realloc(b->len, b->m_seqs * 2);
	}
	b->len[b->n_seqs] = len;
	for (i = 0; i < len; ++i)
		ld_set(b->seq[i], b->n_seqs, seq[len - 1 - i] - 1);
	++b->n_seqs;
}

static pair64_t *set_bwt(bcr_t *bcr, pair64_t *a, int pos)
{
	int64_t k, c[8], m;
	int j, l;
	memset(c, 0, 64);
	if (pos == 0) {
		for (k = 0; k < bcr->n_seqs; ++k) {
			pair64_t *u = &a[k];
			u->u += c[u->v&7], ++c[u->v&7];
		}
	} else {
		for (k = m = 0; k < bcr->n_seqs; ++k) {
			pair64_t *u = &a[k];
			if ((u->v&7) == 0) continue;
			u->u += c[u->v&7], ++c[u->v&7];
			if (m == k) ++m;
			else a[m++] = a[k];
		}
		bcr->n_seqs = m;
	}
	bcr->tot += bcr->n_seqs;
	for (j = 0; j < 6; ++j) bcr->bwt[j].n = c[j];
	for (l = 0; l < 6; ++l) bcr->bwt[0].c[l] = 0;
	for (j = 1; j < 6; ++j)
		for (l = 0; l < 6; ++l)
			bcr->bwt[j].c[l] = bcr->bwt[j-1].e->mc[l];
	for (j = 1; j < 6; ++j)
		for (l = 0; l < 6; ++l)
			bcr->bwt[j].c[l] += bcr->bwt[j-1].c[l];
	memmove(c + 1, c, 40);
	for (k = 1, c[0] = 0; k < 8; ++k) c[k] += c[k - 1]; // NB: MUST BE "8"; otherwise rs_classify_alt() will fail
#if 0
	pair64_t *tmp = malloc(sizeof(pair64_t) * bcr->n_seqs);
	int64_t i[6];
	for (j = 0; j < 6; ++j) i[j] = c[j];
	for (k = 0; k < bcr->n_seqs; ++k)
		tmp[i[a[k].v&7]++] = a[k];
	free(a); a = tmp;
#else
	rs_classify_alt(a, a + bcr->n_seqs, c);
#endif
	for (j = 0; j < 6; ++j)
		bcr->c[j] += c[j], bcr->bwt[j].a = a + c[j];
	for (k = 0; k < bcr->n_seqs; ++k) a[k].u += c[a[k].v&7];
	return a;
}

static void next_bwt(bcr_t *bcr, int class, int pos)
{
	int64_t c[6], k, l;
	rllitr_t ir, iw;
	bucket_t *bwt = &bcr->bwt[class];
	rll_t *ew, *er = bwt->e;

	if (bwt->n == 0) return;
	for (k = bcr->tot, l = 0; k; k >>= 1, ++l);
	if (class) rs_sort(bwt->a, bwt->a + bwt->n, 8, l > 7? l - 7 : 0);
	for (k = 0; k < bwt->n; ++k) {
		pair64_t *u = &bwt->a[k];
		u->u -= k + bcr->c[class];
		u->v = (u->v&~7ULL) | (pos >= bcr->len[u->v>>3]? 0 : ld_get(bcr->seq[pos], u->v>>3) + 1);
	}
	ew = rll_init();
	rll_itr_init(er, &ir);
	rll_itr_init(ew, &iw);
	memset(c, 0, 48);
	for (k = l = 0; k < bwt->n; ++k) {
		pair64_t *u = &bwt->a[k];
		int a = u->v&7;
		if (u->u > l) rll_copy(ew, &iw, er, &ir, u->u - l);
		l = u->u;
		rll_enc(ew, &iw, 1, a);
		u->u = ((ew->mc[a] + iw.l - 1) - c[a]) + bcr->c[a] + bwt->c[a];
		++c[a];
	}
	if (l < er->l) rll_copy(ew, &iw, er, &ir, er->l - l);
	rll_enc_finalize(ew, &iw);
	rll_destroy(er);
	bwt->e = ew;
}
#ifdef HAVE_PTHREAD
static int worker_aux(worker_t *w)
{
	pthread_mutex_lock(&w->bcr->lock);
	while (!w->toproc)
		pthread_cond_wait(&w->bcr->cv, &w->bcr->lock);
	w->toproc = 0;
	pthread_mutex_unlock(&w->bcr->lock);
	next_bwt(w->bcr, w->class, w->pos);
	__sync_fetch_and_add(&w->bcr->proc_cnt, 1);
	return (w->bcr->max_len == w->pos);
}

static void *worker(void *data) { while (worker_aux(data) == 0); return 0; }
#endif
bcr_t *bcr_init(int is_thr)
{
	bcr_t *b;
	int i;
	b = calloc(1, sizeof(bcr_t));
	for (i = 0; i < 6; ++i) b->bwt[i].e = rll_init();
#ifdef HAVE_PTHREAD
	if (is_thr) {
		b->n_threads = 4;
		pthread_mutex_init(&b->lock, 0);
		pthread_cond_init(&b->cv, 0);
		b->tid = calloc(b->n_threads, sizeof(pthread_t)); // tid[0] is not used, as the worker 0 is launched by the master
		b->w = calloc(b->n_threads, sizeof(worker_t));
		for (i = 0; i < b->n_threads; ++i) b->w[i].class = i + 1, b->w[i].bcr = b;
		for (i = 1; i < b->n_threads; ++i) pthread_create(&b->tid[i], 0, worker, &b->w[i]);
	}
#endif
	return b;
}

void bcr_destroy(bcr_t *b)
{
	int i;
#ifdef HAVE_PTHREAD
	if (b->tid) for (i = 1; i < b->n_threads; ++i) pthread_join(b->tid[i], 0);
	pthread_mutex_destroy(&b->lock);
	pthread_cond_destroy(&b->cv);
	free(b->tid); free(b->w);
#endif
	for (i = 0; i < 6; ++i) rll_destroy(b->bwt[i].e);
	free(b->len); free(b->seq);
	free(b);
}

void bcr_build(bcr_t *b)
{
	int64_t k;
	int pos, c;
	pair64_t *a;

	b->m_seqs = b->n_seqs;
	b->len = realloc(b->len, b->n_seqs * 2);
	a = malloc(b->n_seqs * 16);
	for (k = 0; k < b->n_seqs; ++k) a[k].u = 0, a[k].v = k<<3;
	for (pos = 0; pos <= b->max_len; ++pos) {
		a = set_bwt(b, a, pos);
		if (pos) {
#if defined(HAVE_PTHREAD)
			if (b->n_threads > 1) {
				pthread_mutex_lock(&b->lock);
				for (c = 0; c < b->n_threads; ++c)
					b->w[c].toproc = 1, b->w[c].pos = pos;
				b->proc_cnt = 0;
				pthread_cond_broadcast(&b->cv);
				pthread_mutex_unlock(&b->lock);
				worker_aux(&b->w[0]);
				while (b->proc_cnt != b->n_threads);
			} else for (c = 1; c <= 4; ++c) next_bwt(b, c, pos);
#else
			for (c = 1; c <= 4; ++c) next_bwt(b, c, pos);
#endif
		} else next_bwt(b, 0, pos);
		if (pos != b->max_len) ld_destroy(b->seq[pos]);
	}
	free(a);
}

struct bcritr_s {
	const bcr_t *b;
	int c, i;
};

bcritr_t *bcr_itr_init(const bcr_t *b)
{
	bcritr_t *itr;
	itr = calloc(1, sizeof(bcritr_t));
	itr->b = b; itr->i = -1;
	return itr;
}

const uint8_t *bcr_itr_next(bcritr_t *itr, int *l)
{
	rll_t *e;
	const uint8_t *s;
	if (itr->c == 6) return 0;
	++itr->i;
	if (itr->i == itr->b->bwt[itr->c].e->n) {
		if (++itr->c == 6) return 0;
		itr->i = 0;
	}
	e = itr->b->bwt[itr->c].e;
	s = e->z[itr->i];
	if (itr->i == e->n - 1) {
		for (*l = 0; *l < RLL_BLOCK_SIZE; ++*l)
			if (s[*l] == 7) break;
	} else *l = RLL_BLOCK_SIZE;
	return s;
}

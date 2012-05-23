#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "bprope6.h"

#define MP_N_ELEMS 0x10000

typedef struct {
	int size, i;
	int64_t top, max;
	uint8_t **mem;
} mempool_t;

static mempool_t *mp_init(int size)
{
	mempool_t *mp;
	mp = calloc(1, sizeof(mempool_t));
	mp->size = size;
	mp->i = MP_N_ELEMS;
	mp->top = -1;
	return mp;
}

static void mp_destroy(mempool_t *mp)
{
	int64_t i;
	for (i = 0; i <= mp->top; ++i)
		free(mp->mem[i]);
	free(mp->mem); free(mp);
}

static inline void *mp_alloc(mempool_t *mp)
{
	if (mp->i == MP_N_ELEMS) {
		if (++mp->top == mp->max) {
			mp->max = mp->max? mp->max<<1 : 1;
			mp->mem = realloc(mp->mem, sizeof(void*) * mp->max);
		}
		mp->mem[mp->top] = calloc(MP_N_ELEMS, mp->size);
		mp->i = 0;
	}
	return mp->mem[mp->top] + (mp->i++) * mp->size;
}

static int insert_to_leaf(uint8_t *p, int a, int x)
{
#define MAX_RUNLEN 31
#define _insert_after(_n, _s, _i, _b) if ((_i) + 1 != (_n)) memmove(_s+(_i)+2, _s+(_i)+1, (_n)-(_i)-1); _s[(_i)+1] = (_b); ++(_n)

	int r[6], i, l, n = *(int32_t*)p;
	uint8_t *s = p + 4;
	if (n == 0) { // if $s is empty, that is easy
		s[n++] = 1<<3 | a;
		*(int32_t*)p = n;
		return 0;
	}
	memset(r, 0, 24);
	i = l = 0;
	do { // this loop is likely to be the bottleneck
		register int c = s[i++];
		l += c>>3;
		r[c&7] += c>>3;
	} while (l < x);
	assert(i <= n);
	r[s[--i]&7] -= l - x; // $i now points to the left-most run where $a can be inserted
	if (l == x && i != n - 1 && (s[i+1]&7) == a) ++i; // if insert to the end of $i, check if we'd better to the start of ($i+1)
	if ((s[i]&7) == a) { // insert to a long $a run
		if (s[i]>>3 == MAX_RUNLEN) { // the run is full
			for (++i; i != n && (s[i]&7) == a; ++i); // find the end of the long run
			--i;
			if (s[i]>>3 == MAX_RUNLEN) { // then we have to add one run
				_insert_after(n, s, i, 1<<3|a);
			} else s[i] += 1<<3;
		} else s[i] += 1<<3;
	} else if (l == x) { // insert to the end of run; in this case, neither this and the next run is $a
		_insert_after(n, s, i, 1<<3 | a);
	} else if (i != n - 1 && (s[i]&7) == (s[i+1]&7)) { // insert to a long non-$a run
		int rest = l - x, c = s[i]&7;
		s[i] -= rest<<3;
		_insert_after(n, s, i, 1<<3 | a);
		for (i += 2; i != n && (s[i]&7) == c; ++i); // find the end of the long run
		--i;
		if ((s[i]>>3) + rest > MAX_RUNLEN) { // we cannot put $rest to $s[$i]
			rest = (s[i]>>3) + rest - MAX_RUNLEN;
			s[i] = MAX_RUNLEN<<3 | (s[i]&7);
			_insert_after(n, s, i, rest<<3 | c);
		} else s[i] += rest<<3;
	} else { // insert to a short run
		memmove(s + i + 3, s + i + 1, n - i - 1);
		s[i]  -= (l-x)<<3;
		s[i+1] = 1<<3 | a;
		s[i+2] = (l-x)<<3 | (s[i]&7);
		n += 2;
	}
	*(int32_t*)p = n;
	return r[a];
}

typedef struct bpr_node_s {
	struct bpr_node_s *p;
	uint64_t l:54, n:9, is_bottom:1;
	uint64_t c[6];
} bpr_node_t;

struct bprope6_s {
	int max_nodes, max_runs; // both MUST BE even numbers
	uint64_t c[6];
	bpr_node_t *root;
	mempool_t *node, *leaf;
};

static void print_node(const bpr_node_t *p)
{
	if (p->is_bottom) {
		int i, j, k;
		putchar('(');
		for (i = 0; i < p->n; ++i) {
			uint8_t *q = (uint8_t*)p[i].p;
			int n = *(int32_t*)q;
			if (i) putchar(',');
			for (j = 0, q += 4; j < n; ++j)
				for (k = 0; k < q[j]>>3; ++k)
					putchar("$ACGTN"[q[j]&7]);
		}
		putchar(')');
	} else {
		int i;
		putchar('(');
		for (i = 0; i < p->n; ++i) {
			if (i) putchar(',');
			print_node(p[i].p);
		}
		putchar(')');
	}
}

void bpr_print(const bprope6_t *rope) { print_node(rope->root); putchar('\n'); }

static inline bpr_node_t *split_node(bprope6_t *rope, bpr_node_t *u, bpr_node_t *v)
{ // $u: the first node in the bucket; $v: whose child to be split; 
	int j, i = v - u;
	bpr_node_t *w; // this is the uncle
	if (u == 0) { // only happens at the root
		u = v = mp_alloc(rope->node);
		v->n = 1; v->p = rope->root;
		memcpy(v->c, rope->c, 48);
		for (j = 0; j < 6; ++j) v->l += v->c[j];
		rope->root = v;
	}
	if (i != u->n - 1)
		memmove(v + 2, v + 1, sizeof(bpr_node_t) * (u->n - i - 1));
	++u->n;
	w = v + 1;
	memset(w, 0, sizeof(bpr_node_t));
	w->p = mp_alloc(u->is_bottom? rope->leaf : rope->node);
	if (u->is_bottom) {
		uint8_t *p = (uint8_t*)v->p, *q = (uint8_t*)w->p;
		int32_t *np = (int32_t*)p, *nq = (int32_t*)q;
		*nq = *np - (rope->max_runs>>1); *np -= *nq;
		memcpy(q + 4, p + 4 + *np, *nq);
		for (i = 0, q += 4; i < *nq; ++i)
			w->c[q[i]&7] += q[i]>>3;
	} else {
		bpr_node_t *p = v->p, *q = w->p; // $p and $q are cousins
		p->n -= rope->max_nodes>>1;
		memcpy(q, p + p->n, sizeof(bpr_node_t) * (rope->max_nodes>>1));
		q->n = rope->max_nodes>>1;
		q->is_bottom = p->is_bottom;
		for (i = 0; i < q->n; ++i)
			for (j = 0; j < 6; ++j)
				w->c[j] += q[i].c[j];
	}
	for (j = 0; j < 6; ++j)
		w->l += w->c[j], v->c[j] -= w->c[j];
	v->l -= w->l;
	return v;
}

int64_t bpr_insert_symbol(bprope6_t *rope, int a, int64_t x)
{
	bpr_node_t *u = 0, *v = 0, *p = rope->root;
	int64_t y = 0, z;
	int i;
	for (i = 0, z = 0; i < a; ++i) z += rope->c[i];
	do {
		if (p->n == rope->max_nodes) { // node is full; split
			v = split_node(rope, u, v);
			if (y + v->l < x)
				y += v->l, z += v->c[a], ++v, p = v->p;
		}
		for (u = p; y + p->l < x; ++p) y += p->l, z += p->c[a];
		assert(p - u < u->n);
		if (v) ++v->c[a], ++v->l; // we should not change p->c[a] because when this will lead to problems when p's child is split
		v = p; p = p->p;
	} while (!u->is_bottom);
	++v->c[a]; ++v->l;
	++rope->c[a];
	z += insert_to_leaf((uint8_t*)p, a, x - y) + 1;
	if (*(uint32_t*)p + 2 > rope->max_runs)
		split_node(rope, u, v);
//	printf("%c,%lld\t", "$ACGTN"[a], x); bpr_print(rope); fflush(stdout);
	return z;
}

void bpr_insert_string(bprope6_t *rope, int l, uint8_t *str)
{
	uint64_t x = rope->c[0];
	for (--l; l >= 0; --l)
		x = bpr_insert_symbol(rope, str[l], x);
	bpr_insert_symbol(rope, 0, x);
}

bprope6_t *bpr_init(int max_nodes, int max_runs)
{
	bprope6_t *rope;
	rope = calloc(1, sizeof(bprope6_t));
	if (max_runs < 8) max_runs = 8;
	rope->max_nodes= (max_nodes+ 1)>>1<<1;
	rope->max_runs = ((max_runs + 1)>>1<<1) - 4;
	rope->node = mp_init(sizeof(bpr_node_t) * rope->max_nodes);
	rope->leaf = mp_init(rope->max_runs + 4);
	rope->root = mp_alloc(rope->node);
	rope->root->n = 1;
	rope->root->is_bottom = 1;
	rope->root->p = mp_alloc(rope->leaf);
	return rope;
}

void bpr_destroy(bprope6_t *rope)
{
	mp_destroy(rope->node);
	mp_destroy(rope->leaf);
	free(rope);
}

#define MAX_HEIGHT 80

struct bpriter_s {
	const bprope6_t *rope;
	const bpr_node_t *pa[MAX_HEIGHT];
	int k, ia[MAX_HEIGHT];
};

bpriter_t *bpr_iter_init(const bprope6_t *rope)
{
	bpriter_t *i;
	i = calloc(1, sizeof(bpriter_t));
	i->rope = rope;
	for (i->pa[i->k] = rope->root; !i->pa[i->k]->is_bottom;) // descend to the leftmost leaf
		++i->k, i->pa[i->k] = i->pa[i->k - 1]->p;
	return i;
}

const uint8_t *bpr_iter_next(bpriter_t *i, int *n)
{
	const uint8_t *ret;
	if (i->k < 0) return 0;
	*n = *(int32_t*)i->pa[i->k][i->ia[i->k]].p;
	ret = (uint8_t*)i->pa[i->k][i->ia[i->k]].p + 4;
	while (i->k >= 0 && ++i->ia[i->k] == i->pa[i->k]->n) i->ia[i->k--] = 0; // backtracking
	if (i->k >= 0)
		while (!i->pa[i->k]->is_bottom) // descend to the leftmost leaf
			++i->k, i->pa[i->k] = i->pa[i->k - 1][i->ia[i->k - 1]].p;
	return ret;
}

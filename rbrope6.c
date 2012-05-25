#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "rbrope6.h"

/***********************************
 *** Allocation-only memory pool ***
 ***********************************/

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

/******************************
 *** Red-black rope for DNA ***
 ******************************/

#define MAX_HEIGHT 80
#define MAX_RUNLEN 31

typedef struct rbrnode_s {
	union { // IMPORTANT: always make sure x[0].p is the first member; otherwise rb-tree rebalancing will fail
		struct rbrnode_s *p; // pointer to children; internal node only
		uint8_t *s; // string; leaf.x[1] only
		int n; // number of runs; leaf.x[0] only
	} x[2];
	uint64_t c[6];
} node_t;

#define is_red(_p) ((_p)->c[0]&1)
#define set_red(_p) ((_p)->c[0] |= 1)
#define set_black(_p) ((_p)->c[0] &= ~1U)

#define is_leaf(_p) ((_p)->c[1]&1)
#define set_leaf(_p) ((_p)->c[1] |= 1, (_p)->c[0] &= ~1U) // leaves are all black
#define set_internal(_p) ((_p)->c[1] &= ~1U)

#define rbr_strlen(_p) (((_p)->c[0]>>1) + ((_p)->c[1]>>1) + ((_p)->c[2]>>1) + ((_p)->c[3]>>1) + ((_p)->c[4]>>1) + ((_p)->c[5]>>1))

struct rbrope6_s {
	int max_runs;
	mempool_t *node, *str;
	node_t *root;
};

static node_t *rbr_leaf_init(rbrope6_t *rope)
{
	node_t *p;
	p = mp_alloc(rope->node); // $p has been filled with zeros
	set_leaf(p);
	p->x[1].s = mp_alloc(rope->str);
	return p;
}

rbrope6_t *rbr_init(int max_runs)
{
	rbrope6_t *rope;
	rope = calloc(1, sizeof(rbrope6_t));
	if (max_runs < 4) max_runs = 4;
	rope->max_runs = (max_runs + 1)>>1<<1; // make it an even number
	rope->node = mp_init(sizeof(node_t));
	rope->str  = mp_init(rope->max_runs);
	rope->root = rbr_leaf_init(rope);
	return rope;
}

void rbr_destroy(rbrope6_t *rope)
{
	mp_destroy(rope->node);
	mp_destroy(rope->str);
	free(rope);
}

static inline void split_leaf(rbrope6_t *rope, node_t *p)
{
	node_t *q[2];
	uint8_t *s;
	int i;
	q[0] = mp_alloc(rope->node);
	q[1] = rbr_leaf_init(rope);
	// compute q[1]
	s = p->x[1].s;
	memcpy(q[1]->x[1].s, s + (rope->max_runs>>1), rope->max_runs>>1); // copy the later half to q[1]
	memset(s + (rope->max_runs>>1), 0, rope->max_runs>>1); // clear the later half
	q[1]->x[0].n = p->x[0].n - (rope->max_runs>>1);
	for (i = 0, s = q[1]->x[1].s; i < q[1]->x[0].n; ++i) // compute q[1]->c[]
		q[1]->c[s[i]&7] += s[i]>>3<<1;
	// compute q[0]
	memcpy(q[0], p, sizeof(node_t)); // copy everything to q[0], including p->x[0].s and p->c[]
	q[0]->x[0].n = rope->max_runs>>1;
	for (i = 0; i < 6; ++i) q[0]->c[i] -= q[1]->c[i]&(~1ULL);
	// finalize p
	set_internal(p);
	p->x[0].p = q[0]; p->x[1].p = q[1];
}

static int probe_leaf(const node_t *p, int a, int x, int *_i, int *rest)
{
	int r[6], i, l = 0;
	const uint8_t *s = p->x[1].s;
	if (p->x[0].n == 0) {
		*_i = -1; *rest = 0;
		return 0;
	}
	for (i = 0; i < 6; ++i) r[i] = 0;
	do {
		l += *s>>3;
		r[*s&7] += *s>>3;
		++s;
	} while (l < x);
	r[*--s&7] -= l - x; // $i now points to the left-most run where $a can be inserted
	*_i = s - p->x[1].s;
	*rest = l - x;
	return r[a];
}

static void insert_at(node_t *p, int a, int i, int rest)
{
#define _insert_after(_n, _s, _i, _b) if ((_i) + 1 != (_n)) memmove(_s+(_i)+2, _s+(_i)+1, (_n)-(_i)-1); _s[(_i)+1] = (_b); ++(_n)

	uint8_t *s = p->x[1].s;
	if (i < 0) { // p is empty
		s[p->x[0].n++] = 1<<3 | a;
		return;
	}
	if (rest == 0 && i != p->x[0].n - 1 && (s[i+1]&7) == a) ++i; // if insert to the end of $i, check if we'd better to the start of ($i+1)
	if ((s[i]&7) == a) { // insert to a long $a run
		if (s[i]>>3 == MAX_RUNLEN) { // the run is full
			for (++i; i != p->x[0].n && (s[i]&7) == a; ++i); // find the end of the long run
			--i;
			if (s[i]>>3 == MAX_RUNLEN) { // then we have to add one run
				_insert_after(p->x[0].n, s, i, 1<<3|a);
			} else s[i] += 1<<3;
		} else s[i] += 1<<3;
	} else if (rest == 0) { // insert to the end of run; in this case, neither this and the next run is $a
		_insert_after(p->x[0].n, s, i, 1<<3 | a);
	} else if (i != p->x[0].n - 1 && (s[i]&7) == (s[i+1]&7)) { // insert to a long non-$a run
		int c = s[i]&7;
		s[i] -= rest<<3;
		_insert_after(p->x[0].n, s, i, 1<<3 | a);
		for (i += 2; i != p->x[0].n && (s[i]&7) == c; ++i); // find the end of the long run
		--i;
		if ((s[i]>>3) + rest > MAX_RUNLEN) { // we cannot put $rest to $s[$i]
			rest = (s[i]>>3) + rest - MAX_RUNLEN;
			s[i] = MAX_RUNLEN<<3 | (s[i]&7);
			_insert_after(p->x[0].n, s, i, rest<<3 | c);
		} else s[i] += rest<<3;
	} else { // insert to a short run
		memmove(s + i + 3, s + i + 1, p->x[0].n - i - 1);
		s[i]  -= rest<<3;
		s[i+1] = 1<<3 | a;
		s[i+2] = rest<<3 | (s[i]&7);
		p->x[0].n += 2;
	}
}

static inline void update_count(node_t *p) // recompute counts from the two children; p MUST BE internal
{
	p->c[0] = ((p->x[0].p->c[0]>>1) + (p->x[1].p->c[0]>>1))<<1 | (p->c[0]&1);
	p->c[1] = ((p->x[0].p->c[1]>>1) + (p->x[1].p->c[1]>>1))<<1;
	p->c[2] = ((p->x[0].p->c[2]>>1) + (p->x[1].p->c[2]>>1))<<1;
	p->c[3] = ((p->x[0].p->c[3]>>1) + (p->x[1].p->c[3]>>1))<<1;
	p->c[4] = ((p->x[0].p->c[4]>>1) + (p->x[1].p->c[4]>>1))<<1;
	p->c[5] = ((p->x[0].p->c[5]>>1) + (p->x[1].p->c[5]>>1))<<1;
}

static void rbr_print_node(const node_t *p)
{
	if (is_leaf(p)) {
		int i, j;
		const uint8_t *s = p->x[1].s;
		for (i = 0; i < p->x[0].n; ++i)
			for (j = 0; j < s[i]>>3; ++j)
				putchar("$ACGTN"[s[i]&7]);
	} else {
		putchar('(');
		rbr_print_node(p->x[0].p);
		putchar(',');
		rbr_print_node(p->x[1].p);
		putchar(')'); putchar("br"[is_red(p)]);
	}
}

void rbr_print(const rbrope6_t *rope) { rbr_print_node(rope->root); putchar('\n'); }

typedef struct {
	int i, rest, k, a;
	int64_t z;
	node_t *pa[MAX_HEIGHT]; // parents array (pa[0] is special cased)
	int8_t da[MAX_HEIGHT];
} probe1_t;

static int probe_rope(const rbrope6_t *rope, int a, int64_t x, probe1_t *t)
{
	const node_t *p;
	int dir, k, c;
	int64_t y;
	uint8_t *lock;
	for (c = 0, t->z = 0; c < a; ++c) t->z += rope->root->c[c]>>1;
	t->pa[0] = (node_t*)&rope->root; t->da[0] = 0; // this is a trick learnt from libavl
	for (p = rope->root, y = 0, k = 1; !is_leaf(p); p = p->x[dir].p) {
		int l = rbr_strlen(p->x[0].p);
		if (x > l + y) dir = 1, y += l, t->z += p->x[0].p->c[a]>>1;
		else dir = 0;
		t->pa[k] = (node_t*)p;
		t->da[k++] = dir;
	}
	lock = (uint8_t*)p->x[1].s + rope->max_runs - 1;
	if (!__sync_bool_compare_and_swap(lock, 0, 1)) return -1;
	t->pa[k] = (node_t*)p;
	t->z += probe_leaf(p, a, x - y, &t->i, &t->rest) + 1;
	t->k = k; t->a = a;
	return 0;
}

static void update_rope(rbrope6_t *rope, probe1_t *u)
{
	int i, k = u->k;
	for (i = 1; i <= k; ++i) u->pa[i]->c[u->a] += 2;
	insert_at(u->pa[k], u->a, u->i, u->rest);
	u->pa[k]->x[1].s[rope->max_runs - 1] = 0;
	if (u->pa[k]->x[0].n + 3 <= rope->max_runs) return;
	split_leaf(rope, u->pa[k]); set_red(u->pa[k]);
	while (k >= 3 && is_red(u->pa[k - 1])) { // rebalance the red-black tree
		int i = u->da[k - 2], j = !i; // $i: direction of the parent; $j: dir of uncle
		node_t *r = u->pa[k - 2]->x[j].p; // $r points to the uncle
		if (is_red(r)) { // if uncle is red, then grandparent must be black; switch colors and move upwards
			set_black(r);
			set_black(u->pa[k - 1]);
			set_red(u->pa[k - 2]); // grandparent to red
			k -= 2;
		} else {
			node_t *t;
			if (u->da[k - 1] != i) { // if the child and the parent are on different sides: 
				t = u->pa[k - 1]; // $t: parent node
				r = t->x[j].p; // $r: sibling node
				t->x[j].p = r->x[i].p; update_count(t); // rotate to the same side
				r->x[i].p = t; update_count(r);
				u->pa[k - 2]->x[i].p = r;
			} else r = u->pa[k - 1];
			t = u->pa[k - 2];
			set_red(t);
			set_black(r);
			t->x[i].p = r->x[j].p; update_count(t);
			r->x[j].p = t; update_count(r);
			u->pa[k - 3]->x[u->da[k - 3]].p = r; // when k==3, this line will automatically change the root
			break;
		}
	}
	set_black(rope->root); // $root is always black
}

// insert $a after $x characters in $rope and return "|{$rope[i]<$a}| + |{$rope[i]==$a:0<=i<$x}| + 1"
uint64_t rbr_insert_symbol(rbrope6_t *rope, int a, uint64_t x)
{
	int64_t z;
	probe1_t *u;
	u = alloca(sizeof(probe1_t));
	probe_rope(rope, a, x, u);
	z = u->z;
	update_rope(rope, u);
	return z;
}

void rbr_insert_string(rbrope6_t *rope, int l, uint8_t *str)
{
	uint64_t x = rope->root->c[0]>>1;
	for (--l; l >= 0; --l)
		x = rbr_insert_symbol(rope, str[l], x);
	rbr_insert_symbol(rope, 0, x);
}

struct rbriter_s {
	const rbrope6_t *rope;
	const node_t *pa[MAX_HEIGHT];
	int k, da[MAX_HEIGHT];
};

rbriter_t *rbr_iter_init(const rbrope6_t *rope)
{
	rbriter_t *iter;
	const node_t *p;
	iter = calloc(1, sizeof(rbriter_t));
	iter->rope = rope;
	for (p = rope->root; !is_leaf(p); p = p->x[0].p, ++iter->k) // descend to the left-most leaf
		iter->pa[iter->k] = p;
	iter->pa[iter->k] = p; // also add the leaf
	return iter;
}

const uint8_t *rbr_iter_next(rbriter_t *iter, int *n)
{
	const uint8_t *ret;
	if (iter->k < 0) return 0;
	*n = iter->pa[iter->k]->x[0].n;
	ret = iter->pa[iter->k]->x[1].s;
	// find the next leaf
	while (iter->k >= 1 && iter->da[iter->k - 1]) --iter->k;
	if (--iter->k >= 0) {
		const node_t *p = iter->pa[iter->k];
		iter->da[iter->k] = 1;
		p = iter->pa[iter->k++] = p->x[1].p;
		for (; !is_leaf(p); p = p->x[0].p, ++iter->k)
			iter->pa[iter->k] = p, iter->da[iter->k] = 0;
		iter->pa[iter->k] = p;
	}
	return ret;
}

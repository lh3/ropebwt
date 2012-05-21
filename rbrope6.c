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

typedef struct rbrnode_s {
	uint64_t c[6];
	union {
		struct rbrnode_s *p; // pointer to children; internal node only
		size_t n; // number of runs; leaf.x[0] only
		uint8_t *s; // string; leaf.x[1] only
	} x[2];
} rbrnode_t;

#define MAX_HEIGHT 80

#define is_red(_p) ((_p)->c[0]&1)
#define set_red(_p) ((_p)->c[0] |= 1)
#define set_black(_p) ((_p)->c[0] &= ~1U)

#define is_leaf(_p) ((_p)->c[1]&1)
#define set_leaf(_p) ((_p)->c[1] |= 1, (_p)->c[0] |= 1) // leaves are all black
#define set_internal(_p) ((_p)->c[1] &= ~1U)

#define rbr_strlen(_p) (((_p)->c[0]>>1) + ((_p)->c[1]>>1) + ((_p)->c[2]>>1) + ((_p)->c[3]>>1) + ((_p)->c[4]>>1) + ((_p)->c[5]>>1))

struct rbrope6_s {
	int max_runs;
	mempool_t *node, *str;
	rbrnode_t *root;
};

static rbrnode_t *rbr_leaf_init(rbrope6_t *rope)
{
	rbrnode_t *p;
	p = mp_alloc(rope->node); // $p has been filled with zeros
	set_leaf(p);
	p->x[1].s = mp_alloc(rope->str);
	return p;
}

rbrope6_t *rbr_init(int max_runs)
{
	rbrope6_t *rope;
	rope = calloc(1, sizeof(rbrope6_t));
	rope->max_runs = (max_runs + 1)>>1<<1; // make it an even number
	rope->node = mp_init(sizeof(rbrnode_t));
	rope->str  = mp_init(max_runs);
	rope->root = rbr_leaf_init(rope);
	return rope;
}

void rbr_destroy(rbrope6_t *rope)
{
	mp_destroy(rope->node);
	mp_destroy(rope->str);
	free(rope);
}

static void split_leaf(rbrope6_t *rope, rbrnode_t *p)
{
	rbrnode_t *q[2];
	uint8_t *s;
	int i;
	// compute q[1]
	q[0] = mp_alloc(rope->node);
	q[1] = rbr_leaf_init(rope);
	s = p->x[1].s;
	for (i = rope->max_runs>>1; i < p->x[0].n; ++i) // compute q[1]->c[]
		q[1]->c[s[i]&7] += s[i]>>3;
	memcpy(q[1]->x[1].s, s + (rope->max_runs>>1), rope->max_runs>>1); // copy the later half to q[1]
	memset(s + (rope->max_runs>>1), 0, rope->max_runs>>1); // clear the later half
	q[1]->x[0].n = p->x[0].n - (rope->max_runs>>1);
	// compute q[0]
	memcpy(q[0], p, sizeof(rbrnode_t)); // copy everything to q[0], including p->x[0].s and p->c[]
	q[0]->x[0].n = rope->max_runs>>1;
	for (i = 0; i < 6; ++i) q[0]->c[i] -= q[1]->c[i]&(~1ULL);
	// finalize p
	set_internal(p);
	p->x[0].p = q[0]; p->x[1].p = q[1];
}

static int insert_to_leaf(rbrnode_t *p, int a, int x)
{
#define _insert_after(_n, _s, _i, _b) if ((_i) + 1 != (_n)) memmove(_s+(_i)+2, _s+(_i)+1, (_n)-(_i)-1); _s[(_i)+1] = (_b); ++(_n)

	int r[6], i, l;
	uint8_t *s = p->x[1].s;
	if (p->x[0].n == 0) { // if $s is empty, that is easy
		s[p->x[0].n++] = 1<<3 | a;
		return 0;
	}
	memset(r, 0, 24);
	i = l = 0;
	do { // this loop is likely to be the bottleneck
		register int c = s[i++];
		l += c>>3;
		r[c&7] += c>>3;
	} while (l < x);
	assert(i <= p->x[0].n);
	r[s[--i]&7] -= l - x; // $i now points to the left-most run where $a can be inserted
	if (l == x && i != p->x[0].n - 1 && (s[i+1]&7) == a) ++i; // if insert to the end of $i, check if we'd better to the start of ($i+1)
	if ((s[i]&7) == a) { // insert to a long $a run
		if (s[i]>>3 == 31) { // the run is full
			for (++i; i != p->x[0].n && (s[i]&7) == a; ++i); // find the end of the long run
			--i;
			if (s[i]>>3 == 31) { // then we have to add one run
				_insert_after(p->x[0].n, s, i, 1<<3|a);
			} else s[i] += 1<<3;
		} else s[i] += 1<<3;
	} else if (i == p->x[0].n - 1) { // the end of block
		s[p->x[0].n++] = 1<<3 | a;
	} else if (l == x) { // insert to the end of run; in this case, neither this and the next run is $a
		_insert_after(p->x[0].n, s, i, 1<<3 | a);
	} else if ((s[i]&7) == (s[i+1]&7)) { // insert to a long non-$a run; note that $i<$n-1 always stands
		int i0 = i, rest = l - x, c = s[i]&7;
		s[i] -= rest<<3;
		for (++i; i != p->x[0].n && (s[i]&7) == a; ++i); // find the end of the long run
		--i;
		if ((s[i]>>3) + rest > 31) { // we cannot put $rest to $s[$i]
			rest = (s[i]>>3) + rest - 31;
			s[i] |= 31<<3;
			_insert_after(p->x[0].n, s, i, rest<<3 | c);
		} else s[i] += rest<<3;
		_insert_after(p->x[0].n, s, i0, 1<<3 | a);
	} else { // insert to a short run
		memmove(s + i + 3, s + i + 1, p->x[0].n - i - 1);
		s[i]  -= (l-x)<<3;
		s[i+1] = 1<<3 | a;
		s[i+2] = (l-x)<<3 | (s[i]&7);
		p->x[0].n += 2;
	}
	return r[a];
}

static inline void update_count(rbrnode_t *p) // recompute counts from the two children; p MUST BE internal
{
	p->c[0] = ((p->x[0].p->c[0]>>1) + (p->x[1].p->c[0]>>1))<<1;
	p->c[1] = ((p->x[0].p->c[1]>>1) + (p->x[1].p->c[1]>>1))<<1;
	p->c[2] = ((p->x[0].p->c[2]>>1) + (p->x[1].p->c[2]>>1))<<1;
	p->c[3] = ((p->x[0].p->c[3]>>1) + (p->x[1].p->c[3]>>1))<<1;
	p->c[4] = ((p->x[0].p->c[4]>>1) + (p->x[1].p->c[4]>>1))<<1;
	p->c[5] = ((p->x[0].p->c[5]>>1) + (p->x[1].p->c[5]>>1))<<1;
}

// insert $a after $x characters in $rope and return "|{$rope[i]<$a}| + |{$rope[i]==$a:0<=i<$x}| + 1"
uint64_t rbr_insert_symbol(rbrope6_t *rope, int a, uint64_t x)
{
	rbrnode_t *p, *pa[MAX_HEIGHT];
	uint64_t z, y, l;
	int da[MAX_HEIGHT], dir, k, c;

	fprintf(stderr, "%c,%lld\n", "$ACGTN"[a], x);
	for (c = 0, z = 0; c < a; ++c) z += rope->root->c[c]>>1; // $z equals the number of symbols smaller than $a
	// pinpoint the node where $a is inserted
	pa[0] = rope->root, da[0] = -1;
	for (p = rope->root, y = 0, k = 1; !is_leaf(p); p = p->x[dir].p) {
		l = rbr_strlen(p->x[0].p);
		if (x > l + y) dir = 1, y += l, z += p->c[a];
		else dir = 0;
		pa[k] = p;
		da[k++] = dir;
		p->c[a] += 2;
	}
	p->c[a] += 2; // the leaf count has not been updated
	z += insert_to_leaf(p, a, x - y) + 1; // NB: $p always has enough room for one insert; +1 to include $rope[$x], which equals $a
	if (p->x[0].n + 2 <= rope->max_runs) return z;
	// we need to split $p and rebalance the red-black rope
	split_leaf(rope, p); set_red(p);
	while (k >= 3 && is_red(pa[k - 1])) {
		int i = da[k - 2], j = !i; // $i: direction of the parent; $j: dir of uncle
		rbrnode_t *r = pa[k - 2]->x[j].p; // $r points to the uncle
		if (is_red(r)) { // if uncle is red, then grandparent must be black; switch colors and move upwards
			set_black(r);
			set_black(pa[k - 1]);
			set_red(pa[k - 2]); // grandparent to red
			k -= 2;
		} else {
			rbrnode_t *t;
			if (da[k - 1] != i) { // if the child and the parent are on different sides: 
				t = pa[k - 1]; // $t: parent node
				r = t->x[j].p; // $r: sibling node
				t->x[j].p = r->x[i].p; update_count(t); // rotate to the same side
				r->x[i].p = t; update_count(r);
				pa[k - 2]->x[i].p = r;
			} else r = pa[k - 1];
			t = pa[k - 2];
			set_red(t);
			set_black(r);
			t->x[i].p = r->x[j].p; update_count(t);
			r->x[j].p = t; update_count(r);
			pa[k - 3]->x[da[k - 3]].p = r;
			break;
		}
	}
	set_black(rope->root); // $root is always black
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
	const rbrnode_t *pa[MAX_HEIGHT];
	int k, da[MAX_HEIGHT];
};

rbriter_t *rbr_iter_init(const rbrope6_t *rope)
{
	rbriter_t *iter;
	const rbrnode_t *p;
	iter = calloc(1, sizeof(rbriter_t));
	iter->rope = rope;
	for (p = rope->root; !is_leaf(p); p = p->x[0].p, ++iter->k) // descend to the left-most leaf
		iter->pa[iter->k] = p;
	iter->pa[iter->k] = p; // also add the leaf
	return iter;
}

const uint8_t *rbr_iter_next(rbriter_t *iter, int *n, int *l)
{
	const uint8_t *ret;
	if (iter->k < 0) return 0;
	*l = rbr_strlen(iter->pa[iter->k]);
	*n = iter->pa[iter->k]->x[0].n;
	ret = iter->pa[iter->k]->x[1].s;
	// find the next leaf
	while (iter->k >= 1 && iter->da[iter->k - 1]) --iter->k;
	if (--iter->k >= 0) {
		const rbrnode_t *p = iter->pa[iter->k];
		iter->da[iter->k] = 1;
		p = iter->pa[iter->k++] = p->x[1].p;
		for (; !is_leaf(p); p = p->x[0].p, ++iter->k)
			iter->pa[iter->k] = p, iter->da[iter->k] = 0;
		iter->pa[iter->k] = p;
	}
	return ret;
}

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "ibwt.h"

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
			mp->mem = realloc(mp->mem, sizeof(void*));
		}
		mp->mem[mp->top] = calloc(MP_N_ELEMS, mp->size);
		mp->i = 0;
	}
	return mp->mem[mp->top] + (mp->i++) * mp->size;
}

/******************************
 *** Red-black rope for DNA ***
 ******************************/

typedef struct {
	uint64_t c[6];
	union {
		rbr_node_t *p; // pointer to children; internal node only
		size_t n; // length of RLE string; leaf.x[0] only
		uint8_t *s; // string; leaf.x[1] only
	} x[2];
} rbr_node_t;

#define MAX_HEIGHT 64
#define MAX_RUNS 512 // NB: this MUST BE an even number

#define is_red(_p) ((_p)->c[0]&1)
#define set_red(_p) ((_p)->c[0] |= 1)
#define set_black(_p) ((_p)->c[0] &= ~1U)

#define is_leaf(_p) ((_p)->c[1]&1)
#define set_leaf(_p) ((_p)->c[1] |= 1, (_p)->c[0] |= 1) // leaves are all black
#define set_internal(_p) ((_p)->c[1] &= ~1U)

#define leaf_len(_p) ((size_t)(_p)->p[1])

#define rbr_strlen(_p) (((_p)->c[0]>>1) + ((_p)->c[1]>>1) + ((_p)->c[2]>>1) + ((_p)->c[3]>>1) + ((_p)->c[4]>>1) + ((_p)->c[5]>>1))

typedef struct {
	mempool_t *node, *str;
	rbr_node_t *root;
} rbr_t;

static rbr_node_t *rbr_node_init(rbr_t *rope, int is_leaf)
{
	rbr_node_t *p;
	p = mp_alloc(rope->node); // $p has been initialized as zeros
	if (is_leaf) {
		set_leaf(p);
		p->x[1].s = mp_alloc(rope->str);
	}
	return p;
}

static rbr_t *rbr_init(void)
{
	rbr_t *rope;
	rope = calloc(1, sizeof(rbr_t));
	rope->node = mp_init(sizeof(rbr_node_t));
	rope->str = mp_init(MAX_LEAF_LEN);
	rope->root = rbr_node_init(rope, 1);
	return p;
}

static void rbr_destroy(rbr_t *rope)
{
	mp_destroy(mp->node);
	mp_destroy(mp->str);
	free(rope);
}

static void rbr_split(rbr_t *rope, rbr_node_t *p)
{
	int i;
	rbr_node_t *q[2];
	uint8_t *s;
	q[0] = mp_alloc(rope->node);
	memcpy(q[0], p, sizeof(rbr_node_t)); // copy everything to q[0], including p->p[0] and p->c[]
	q[1] = rbr_node_init(rope, 1);
	set_internal(p); set_red(p);
	p->x[0].p = q[0]; p->x[1].p = q[1];
	// set counts and length
	s = q[0]->x[1].s;
	for (i = MAX_LEAF_LEN>>1; i < MAX_LEAF_LEN; ++i) // compute q[1]->c[]
		q[1]->c[s[i]&7] += s[i]>>3;
	memcpy(q[1]->x[1].s, s + (MAX_LEAF_LEN>>1), MAX_LEAF_LEN>>1);
	q[0]->x[0].l = q[1]->x[0].l = MAX_LEAF_LEN>>1;
	for (i = 0; i < 6; ++i)
		q[0]->c[i] -= q[1]->c[i];
}

static int insert_to_block(int n, uint8_t *s, int a, int x)
{
#define _insert_after(_n, _s, _i, _b) if ((_i) + 1 != (_n)) memmove(_s+(_i)+2, _s+(_i)+1, (_n)-(_i)-1); _s[(_i)+1] = (_b); ++(_n)

	int r[6], i, l;
	if (n == 0) {
		s[n++] = 1<<3 | a;
		return 0;
	}
	memset(r, 0, 24);
	i = l = 0;
	do {
		register int c = s[i++];
		l += c>>3;
		r[c&7] += c>>3;
	} while (l < x);
	assert(i <= n);
	--i; // points to the first run where $a is possibly inserted to
	r[s[i]&7] -= l - x;
	if (l == x && i != n - 1 && (s[i+1]&7) == a) ++i; // if to the end of $i, then see if we'd better to the start of ($i+1)
	if ((s[i]&7) == a) { // insert to a long $a run
		if (s[i]>>3 == 31) {
			for (++i; i != n && (s[i]&7) == a; ++i); // find the end of the long run
			--i;
			if (s[i]>>3 == 31) { // then we have to add one run
				_insert_after(n, s, i, 1<<3|a);
			} else s[i] += 1<<3;
		} else s[i] += 1<<3;
	} else if (i == n - 1) { // the end of block
		s[i] = 1<<3 | a; ++n;
	} else if (l == x) { // insert to the end of run; in this case, neither this and the next run is $a
		_insert_after(n, s, i, 1<<3 | a);
	} else if ((s[i]&7) == (s[i+1]&7)) { // insert to a long non-$a run; note that $i<$n-1 always stands
		int rest = l - x, c = s[i]&7;
		s[i] -= rest<<3;
		_insert_after(n, s, i, 1<<3 | a);
		for (++i; i != n && (s[i]&7) == a; ++i); // find the end of the long run
		--i;
		if ((s[i]>>3) + rest > 31) {
			rest = (s[i]>>3) + rest - 31;
			s[i] |= 0xf8;
			_insert_after(n, s, i, rest<<3 | c);
		} else s[i] += rest<<3;
	} else { // insert to a short run
		memmove(s + i + 3, s + i + 1, n - i - 1);
		s[i]  -= (l-x)<<3;
		s[i+1] = 1<<3 | a;
		s[i+2] = (l-x)<<3 | (s[i]&7);
		n += 2;
	}
	return r[a];
}

static inline void update_count(rbr_node_t *p) // recompute counts from the two children; p MUST BE internal
{
	p->c[0] = ((p->p[0]->c[0]>>1) + (p->p[1]->c[0]>>1))<<1;
	p->c[1] = ((p->p[0]->c[1]>>1) + (p->p[1]->c[1]>>1))<<1;
	p->c[2] = ((p->p[0]->c[2]>>1) + (p->p[1]->c[2]>>1))<<1;
	p->c[3] = ((p->p[0]->c[3]>>1) + (p->p[1]->c[3]>>1))<<1;
	p->c[4] = ((p->p[0]->c[4]>>1) + (p->p[1]->c[4]>>1))<<1;
	p->c[5] = ((p->p[0]->c[5]>>1) + (p->p[1]->c[5]>>1))<<1;
}

// insert $a after $x characters in $rope and return the number of $a in $rope[0..$x-1]
uint64_t rbr_insert(rbr_t *rope, int a, uint64_t x)
{
	rbr_node_t *p, *q[2], *pa[MAX_HEIGHT];
	uint64_t z, y, l;
	int da[MAX_HEIGHT], dir, k;
	rbr_note_t *root = rope->root;
	// pinpoint the node where $a is inserted
	pa[0] = 0, da[0] = 0;
	for (p = root, y = 0, k = 1; !is_leaf(p); p = p->p[dir]) {
		l = rbr_strlen(p->x[0].p);
		if (x > l + y) dir = 1, y += l, z += p->c[a];
		else dir = 0;
		pa[k] = p;
		da[k++] = dir;
		++p->c[a];
	}
	++p->c[a]; // the leaf count has not been updated
	z += insert_to_block(p->x[0].n, p->x[1].s, a, x - l);
	if (p->x[0].n + 2 <= MAX_RUNS) return z;
	// now we need to split $p and rebalance the red-black rope
	rbr_split(rope, p);
	while (k >= 3 && is_red(pa[k - 1])) {
		int i = da[k - 2], j = !i; // $i: direction of the parent; $j: dir of uncle
		rbr_node_t *r = pa[k - 2]->p[j]; // $r points to the uncle
		if (is_red(r)) { // if uncle is red, then grandparent must be black; switch colors and move upwards
			set_black(r);
			set_black(pa[k - 1]);
			set_red(pa[k - 2]); // grandparent to red
			k -= 2;
		} else {
			rbr_node_t *t;
			if (da[k - 1] != i) { // if the child and the parent are on different sides: 
				t = pa[k - 1]; // $t: parent node
				r = t->p[j]; // $r: sibling node
				t->x[j].p = r->x[i].p; update_count(t); // rotate to the same side
				r->x[i].p = t; update_count(r);
				pa[k - 2]->x[i].p = r;
			} else r = pa[k - 1];
			t = pa[k - 2];
			set_red(t);
			set_black(r);
			t->x[i].p = r->x[j].p; update_count(t);
			r->x[j].p = t; update_count(r);
			pa[k - 3]->p[da[k - 3]] = r;
			break;
		}
	}
	set_black(root); // $root is always black
	return z;
}

/******************************************
 *** External APIs for constructing BWT ***
 ******************************************/

struct ibwt_s {
	uint64_t c[6];
	rbr_t *rope;
};

ibwt_t *ib_init()
{
	ibwt_t *ib;
	int c;
	ib = calloc(1, sizeof(ibwt_t));
	ib->rope = rbr_init();
	return ib;
}

void ib_destroy(ibwt_t *ib)
{
	int c;
	rbr_destroy(ib->rope[c]);
	free(ib);
}

uint64_t ib_push(ibwt_t *ib, int a, uint64_t x)
{
	int c;
	uint64_t y;
	y = rbr_insert(ib->rope, a, x);
	for (c = 0; c < a; ++c) y += ib->c[c];
	++ib->c[a];
	return y + 1;
}

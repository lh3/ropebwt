#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include "ksort.h"
#include "rbrope6-mt.h"

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
#define MIN_RECODE 10

typedef struct rbrnode_s {
	union {
		struct rbrnode_s *p; // pointer to children; internal node only
		uint8_t *s; // string; leaf.x[1] only
		int64_t n; // number of runs; leaf.x[0] only
	} x[2];
	struct rbrnode_s *parent;
	uint64_t c[6];
} node_t;

#define is_red(_p) ((_p)->c[0]&1)
#define set_red(_p) ((_p)->c[0] |= 1)
#define set_black(_p) ((_p)->c[0] &= ~1U)

#define is_leaf(_p) ((_p)->c[1]&1)
#define set_leaf(_p) ((_p)->c[1] |= 1, (_p)->c[0] &= ~1U) // leaves are all black
#define set_internal(_p) ((_p)->c[1] &= ~1U)

#define rbm_strlen(_p) (((_p)->c[0]>>1) + ((_p)->c[1]>>1) + ((_p)->c[2]>>1) + ((_p)->c[3]>>1) + ((_p)->c[4]>>1) + ((_p)->c[5]>>1))

typedef struct {
	node_t *p;
	uint32_t pos, i;
	uint64_t z:61, a:3;
} probe1_t;

#define probe_lt(a, b) ((a).z < (b).z)
KSORT_INIT(rbm, probe1_t, probe_lt)

typedef struct {
	int start, step, n, toproc;
	struct rbmope6_s *rope;
} worker_t;

struct rbmope6_s {
	int max_runs, n_threads, max_seqs, n_seqs;
	int *len;
	uint8_t **buf;
	probe1_t *u;
	mempool_t *node, *str;
	node_t *root;
	// for mt
	pthread_t *tid;
	worker_t *w;
	pthread_mutex_t lock;
	pthread_cond_t cv;
	volatile int proc_cnt;
	int done;
};

/**************
 *** Worker ***
 **************/

static int probe_leaf(const node_t *p, int a, int x, uint32_t *pos);
static int probe(const rbmope6_t *rope, probe1_t *u);
void rbm_update(rbmope6_t *rope);

static int worker_aux(worker_t *w)
{
	int i, stop = 0;
	pthread_mutex_lock(&w->rope->lock);
	while (!w->toproc && !w->rope->done)
		pthread_cond_wait(&w->rope->cv, &w->rope->lock);
	if (w->rope->done) stop = 1;
	w->toproc = 0;
	pthread_mutex_unlock(&w->rope->lock);
	if (stop) return 1; // to quit the thread
	for (i = w->start; i < w->n; i += w->step) {
		probe1_t *u = &w->rope->u[i];
		probe(w->rope, u);
		u->z += probe_leaf(u->p, u->a, u->pos, &u->pos) + w->rope->n_seqs;
	}
	__sync_fetch_and_add(&w->rope->proc_cnt, 1);
	return 0;
}

static void *worker(void *data)
{
	while (worker_aux(data) == 0);
	return 0;
}

void rbm_finalize(rbmope6_t *rope)
{
	int i;
	if (rope->done) return;
	if (rope->n_seqs) rbm_update(rope);
	if (rope->tid) {
		pthread_mutex_lock(&rope->lock);
		rope->done = 1; rope->proc_cnt = 0;
		pthread_cond_broadcast(&rope->cv);
		pthread_mutex_unlock(&rope->lock);
		for (i = 1; i < rope->n_threads; ++i) pthread_join(rope->tid[i], 0);
		pthread_cond_destroy(&rope->cv);
		pthread_mutex_destroy(&rope->lock);
		free(rope->tid); free(rope->w);
		rope->tid = 0; rope->w = 0;
	} else rope->done = 1;
}

/********************
 *** Key routines ***
 ********************/

static node_t *rbm_leaf_init(rbmope6_t *rope)
{
	node_t *p;
	p = mp_alloc(rope->node); // $p has been filled with zeros
	set_leaf(p);
	p->x[1].s = mp_alloc(rope->str);
	return p;
}

rbmope6_t *rbm_init(int n_threads, int max_seqs, int max_runs)
{
	rbmope6_t *rope;
	rope = calloc(1, sizeof(rbmope6_t));
	if (max_runs < 5) max_runs = 5; // cannot be smaller; otherwise fix() will not work
	rope->n_threads = n_threads;
	rope->max_seqs = max_seqs;
	rope->buf = malloc(sizeof(void*) * max_seqs);
	rope->len = malloc(sizeof(int) * max_seqs);
	rope->u = malloc(sizeof(probe1_t) * max_seqs);
	rope->max_runs = (max_runs + 1)>>1<<1; // make it an even number
	rope->node = mp_init(sizeof(node_t));
	rope->str  = mp_init(rope->max_runs);
	rope->root = rbm_leaf_init(rope);
	if (n_threads > 1) {
		int i;
		pthread_mutex_init(&rope->lock, 0);
		pthread_cond_init(&rope->cv, 0);
		rope->tid = calloc(n_threads, sizeof(pthread_t)); // tid[0] is not used, as the worker 0 is launched by the master
		rope->w = calloc(n_threads, sizeof(worker_t));
		for (i = 0; i < n_threads; ++i)
			rope->w[i].start = i, rope->w[i].step = n_threads, rope->w[i].rope = rope;
		for (i = 1; i < n_threads; ++i) pthread_create(&rope->tid[i], 0, worker, &rope->w[i]);
	}
	return rope;
}

void rbm_destroy(rbmope6_t *rope)
{
	int i;
	if (!rope->done) rbm_finalize(rope);
	for (i = 0; i < rope->n_seqs; ++i) free(rope->buf[i]);
	free(rope->buf); free(rope->len); free(rope->u);
	mp_destroy(rope->node);
	mp_destroy(rope->str);
	free(rope);
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

static inline node_t *insert_fix(node_t *q)
{
	while (q && q->parent) {
		node_t *p = q->parent, *u, *g, *gg;
		int i, j;
		if (!is_red(p)) break; // p: parent
		if ((g = p->parent) == 0) break; // g: grandparent
		i = (g->x[1].p == p); j = !i;
		u = g->x[j].p;
		if (is_red(u)) {
			set_black(u);
			set_black(p);
			set_red(g);
			q = g;
		} else {
			if (p->x[i].p != q) {
				node_t *s = p->x[j].p;
				p->x[j].p = s->x[i].p; s->x[i].p->parent = p; update_count(p);
				s->x[i].p = p; p->parent = s; update_count(s);
				g->x[i].p = s; s->parent = g;
				p = s;
			}
			gg = g->parent;
			set_red(g);
			set_black(p);
			g->x[i].p = p->x[j].p; p->x[j].p->parent = g; update_count(g);
			p->x[j].p = g; g->parent = p; update_count(p);
			if (gg == 0) {
				p->parent = 0;
				set_black(p);
				return p;
			} else gg->x[(gg->x[1].p == g)].p = p, p->parent = gg;
			break;
		}
	}
	return 0;
}
// split bottom node $p if it is too long and rebalance the tree if necessary; if $s0 is not NULL, p->x[1].s will be freed and then directed to s0
static inline void fix(rbmope6_t *rope, node_t *p, uint8_t *s0)
{
	while (p->x[0].n + 2 > rope->max_runs) {
		int i, n1 = p->x[0].n >> 1;
		node_t *q[2], *t;
		uint8_t *s;
		if (n1 + 2 > rope->max_runs) n1 = (rope->max_runs>>1) + (rope->max_runs>>2);
		q[0] = mp_alloc(rope->node);
		q[1] = rbm_leaf_init(rope);
		// get q[1]
		s = p->x[1].s;
		memcpy(q[1]->x[1].s, s + (p->x[0].n - n1), n1); // copy the later half to q[1]
		memset(s + (p->x[0].n - n1), 0, n1); // clear the later half
		q[1]->x[0].n = n1;
		for (i = 0, s = q[1]->x[1].s; i < n1; ++i) // compute q[1]->c[]
			q[1]->c[s[i]&7] += s[i]>>3<<1;
		// get q[0]
		*q[0] = *p; // copy everything to q[0], including p->x[0].s and p->c[]
		q[0]->x[0].n -= n1;
		for (i = 0; i < 6; ++i) q[0]->c[i] -= q[1]->c[i]&(~1ULL);
		// update p and rebalance
		set_internal(p); set_red(p);
		p->x[0].p = q[0]; q[0]->parent = p;
		p->x[1].p = q[1]; q[1]->parent = p;
		if ((t = insert_fix(p)) != 0) rope->root = t;
		p = q[0];
	}
	if (s0) {
		memcpy(s0, p->x[1].s, p->x[0].n);
		free(p->x[1].s);
		p->x[1].s = s0;
	}
}

static int probe_leaf(const node_t *p, int a, int x, uint32_t *pos)
{
	int r[6], i, l = 0, len;
	const uint8_t *s = p->x[1].s;
	if (p->x[0].n == 0) {
		*pos = 0;
		return 0;
	}
	len = rbm_strlen(p);
	if (x < len>>1) { // forward search
		for (i = 0; i < 6; ++i) r[i] = 0;
		do {
			l += *s>>3;
			r[*s&7] += *s>>3;
			++s;
		} while (l < x);
	} else { // backward search
		__builtin_prefetch(p->x[1].p, 0);
		for (i = 0; i < 6; ++i) r[i] = p->c[i]>>1;
		l = len, s += p->x[0].n;
		do {
			--s;
			l -= *s>>3;
			r[*s&7] -= *s>>3;
		} while (l >= x);
		l += *s>>3; r[*s&7] += *s>>3; ++s;
	}
	r[*--s&7] -= l - x; // $s now points to the left-most run where $a can be inserted
	*pos = (uint32_t)(s - p->x[1].s)<<16 | ((*s>>3) - (l - x));
	return r[a];
}

static void insert_at(node_t *p, int a, uint32_t pos)
{
#define _insert_after(_n, _s, _i, _b) if ((_i) + 1 != (_n)) memmove(_s+(_i)+2, _s+(_i)+1, (_n)-(_i)-1); _s[(_i)+1] = (_b); ++(_n)

	uint8_t *s = p->x[1].s;
	int i, rest;
	if (p->x[0].n == 0) { // p is empty
		assert(pos == 0);
		s[p->x[0].n++] = 1<<3 | a;
		return;
	}
	i = pos>>16; rest = (s[i]>>3) - (pos&0xffff);
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

static void rbm_print_node(const node_t *p)
{
	if (is_leaf(p)) {
		int i, j;
		const uint8_t *s = p->x[1].s;
		for (i = 0; i < p->x[0].n; ++i)
			for (j = 0; j < s[i]>>3; ++j)
				putchar("$ACGTN"[s[i]&7]);
	} else {
		putchar('(');
		rbm_print_node(p->x[0].p);
		putchar(',');
		rbm_print_node(p->x[1].p);
		putchar(')'); putchar("br"[is_red(p)]);
	}
}

void rbm_print(const rbmope6_t *rope) { rbm_print_node(rope->root); putchar('\n'); }

static int probe(const rbmope6_t *rope, probe1_t *u)
{
	const node_t *p;
	int dir, c;
	int64_t x = u->z, y;
	for (c = 0, u->z = 0; c < u->a; ++c) u->z += rope->root->c[c]>>1;
	for (p = rope->root, y = 0; !is_leaf(p); p = p->x[dir].p) {
		int l = rbm_strlen(p->x[0].p);
		if (x > l + y) dir = 1, y += l, u->z += p->x[0].p->c[u->a]>>1;
		else dir = 0;
	}
	u->p = (node_t*)p; u->pos = x - y;
	return 0;
}

typedef struct {
	int last_l, last_c;
	uint8_t *s;
} enc1_t;

static inline void enc1(enc1_t *e, int l, int c)
{
	if (c != e->last_c) {
		if (e->last_l) {
			while (e->last_l > 31) {
				*e->s++ = 31<<3 | e->last_c;
				e->last_l -= 31;
			}
			*e->s++ = e->last_l<<3 | e->last_c;
		}
		e->last_l = l; e->last_c = c;
	} else e->last_l += l;
}

static void modify_multi1(rbmope6_t *rope, int n, probe1_t *u)
{
	int i, j, c[6];
	node_t *p;
	memset(c, 0, 6 * sizeof(int));
	for (i = 0; i < n; ++i) c[u[i].a] += 2;
	for (p = u->p; p; p = p->parent)
		for (i = 0; i < 6; ++i) p->c[i] += c[i];
	p = u->p;
	if (n == 1) {
		insert_at(p, u->a, u->pos);
		fix(rope, p, 0);
	} else if (n < MIN_RECODE && (n<<1) + p->x[0].n <= rope->max_runs) {
		for (i = n - 1; i >= 0; --i)
			insert_at(p, u[i].a, u[i].pos);
		fix(rope, p, 0);
	} else { // the generic case
		enc1_t e;
		uint8_t *s, *t = p->x[1].s;
		e.s = s = malloc(n + rope->max_runs + (n > rope->max_runs? n : rope->max_runs) * 2 + 2); // s will be freed by fix()
		e.last_c = -1; e.last_l = 0;
		if (p->x[0].n) {
			for (i = j = 0; i < p->x[0].n; ++i) {
				if (j < n && i == u[j].pos>>16) {
					int l = 0;
					do {
						enc1(&e, (u[j].pos&0xffff) - l, t[i]&7);
						enc1(&e, 1, u[j].a);
						l = u[j++].pos&0xffff;
					} while (j < n && u[j].pos>>16 == i);
					enc1(&e, (p->x[1].s[i]>>3) - l, t[i]&7);
				} else enc1(&e, t[i]>>3, t[i]&7);
			}
		} else for (j = 0; j < n; ++j) enc1(&e, 1, u[j].a);
		enc1(&e, 0, -2);
		p->x[0].n = e.s - s;
		t = p->x[1].s; p->x[1].s = s;
		fix(rope, p, t);
	}
}

static void modify_multi(rbmope6_t *rope, int n, probe1_t *u)
{
	node_t *p = 0;
	int i, last = 0;
	for (i = 0; i < n; ++i) {
		if (u[i].p == p) continue;
		if (i) modify_multi1(rope, i - last, &u[last]);
		last = i; p = u[i].p;
	}
	modify_multi1(rope, i - last, &u[last]);
}

void rbm_update(rbmope6_t *rope)
{ // FIXME: *NOT* WORKING when input sequences are of different lengths!!!
	int i, l, m, n;
	for (i = 0; i < rope->n_seqs; ++i) {
		probe1_t *u = &rope->u[i];
		u->i = i;
		u->a = rope->buf[i][0];
		u->z = rope->root->c[0]>>1;
	}
	n = rope->n_seqs;
	for (l = 1; n; ++l) {
		int64_t c[6];
		// probe the insertion point and compute the pre-coordinate for the next insertion
		if (rope->n_threads > 1) {
			pthread_mutex_lock(&rope->lock);
			for (i = 0; i < rope->n_threads; ++i)
				rope->w[i].toproc = 1, rope->w[i].n = n;
			rope->proc_cnt = 0;
			pthread_cond_broadcast(&rope->cv);
			pthread_mutex_unlock(&rope->lock);
			worker_aux(&rope->w[0]);
			while (rope->proc_cnt < rope->n_threads);
		} else {
			for (i = 0; i < n; ++i) {
				probe1_t *u = &rope->u[i];
				probe(rope, u);
				u->z += probe_leaf(u->p, u->a, u->pos, &u->pos) + rope->n_seqs;
			}
		}
		// perform insertion
		modify_multi(rope, n, rope->u);
		// compute the insertion coordinate in the new BWT
		memset(c, 0, 48);
		for (i = 0; i < n; ++i) {
			probe1_t *u = &rope->u[i];
			u->z += c[u->a];
			++c[u->a];
		}
		memmove(c + 1, c, 40);
		for (i = 1, c[0] = 0; i < 6; ++i) c[i] += c[i - 1]; // c[] now keeps the accumulative counts
		for (i = 0; i < n; ++i) rope->u[i].z += c[rope->u[i].a];
		// compute the insertion coordinate in the current BWT
		ks_introsort(rbm, n, rope->u);
		for (i = 0; i < n; ++i) rope->u[i].z -= i;
		// get the base to insert
		for (i = 0; i < n; ++i) {
			probe1_t *u = &rope->u[i];
			u->a = l < rope->len[u->i]? rope->buf[u->i][l] : l == rope->len[u->i]? 0 : 7;
		}
		// exclude finished sequences
		for (i = m = 0; i < n; ++i)
			if (rope->u[i].a != 7) {
				if (m == i) ++m;
				else rope->u[m++] = rope->u[i];
			}
		assert(m == 0 || m == n); // assertion failure if input are of different lengths
		n = m;
	}
	for (i = 0; i < rope->n_seqs; ++i) free(rope->buf[i]);
	rope->n_seqs = 0;
}

void rbm_insert_string(rbmope6_t *rope, int l, uint8_t *str)
{
	int i;
	uint8_t *s;
	s = rope->buf[rope->n_seqs] = malloc(l);
	memcpy(s, str, l);
	for (i = 0; i < l>>1; ++i) { // reverse
		int c = s[i];
		s[i] = s[l - 1 - i];
		s[l - 1 - i] = c;
	}
	rope->len[rope->n_seqs++] = l;
	if (rope->n_seqs == rope->max_seqs) rbm_update(rope);
}

struct rbmiter_s {
	const rbmope6_t *rope;
	const node_t *pa[MAX_HEIGHT];
	int k, da[MAX_HEIGHT];
};

rbmiter_t *rbm_iter_init(const rbmope6_t *rope)
{
	rbmiter_t *iter;
	const node_t *p;
	iter = calloc(1, sizeof(rbmiter_t));
	iter->rope = rope;
	for (p = rope->root; !is_leaf(p); p = p->x[0].p, ++iter->k) // descend to the left-most leaf
		iter->pa[iter->k] = p;
	iter->pa[iter->k] = p; // also add the leaf
	return iter;
}

const uint8_t *rbm_iter_next(rbmiter_t *iter, int *n)
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

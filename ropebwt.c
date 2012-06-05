#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

/*********************
 *** Merge indices ***
 *********************/

#include "rld.h"

#define BLOCK_SIZE 0x40000
#define TIMER_INTV 64

typedef struct {
	int start, step;
	int64_t offset, shift;
	const rld_t *e0, *e1;
	uint64_t *bits;
	int64_t *buf;
} worker_t;

static inline void update_bits(int n, const int64_t *buf, uint64_t *bits)
{
	const int64_t *q, *end = buf + n;
	for (q = buf; q != end; ++q) {
		uint64_t *p = bits + (*q>>6);
		uint64_t x = 1ull<<(*q&0x3f);
		__sync_or_and_fetch(p, x); // SEE ALSO: http://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Atomic-Builtins.html
	}
}

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int n = 0;
	int64_t i, k, x, y;
	uint64_t *ok;
	ok = alloca(8 * w->e0->asize);
	k = x = w->start;
	y = w->offset + w->shift * (w->start + 1) - 1;
	i = w->shift? y : w->e0->mcnt[1] - 1;
	w->buf[n++] = i + k + 1;
	for (;;) {
		int c = rld_rank1a(w->e1, k, ok);
		if (c == 0) {
			x += w->step; y += w->shift * w->step;
			if (x >= w->e1->mcnt[1]) break;
			k = x;
			i = w->shift? y : w->e0->mcnt[1] - 1;
		} else {
			k = w->e1->cnt[c] + ok[c] - 1;
			rld_rank1a(w->e0, i, ok);
			i = w->e0->cnt[c] + ok[c] - 1;
		}
		if (n == BLOCK_SIZE) {
			update_bits(n, w->buf, w->bits);
			n = 0;
		}
		w->buf[n++] = i + k + 1;
	}
	if (n) update_bits(n, w->buf, w->bits);
	return 0;
}

static uint64_t *compute_bits(const rld_t *e0, const rld_t *e1, int n_threads, int64_t offset, int64_t shift)
{
	uint64_t *bits;
	pthread_t *tid;
	worker_t *w;
	int j;

	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	if ((bits = calloc((e0->mcnt[0] + e1->mcnt[0] + 63) / 64, 8)) == 0) return 0;
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e0 = e0; ww->e1 = e1;
		ww->step = n_threads;
		ww->start = j;
		ww->buf = malloc(BLOCK_SIZE * 8);
		ww->offset = offset;
		ww->shift = shift;
		ww->bits = bits;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], 0, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	for (j = 0; j < n_threads; ++j) free(w[j].buf);
	free(w); free(tid);
	return bits;
}

rld_t *rld_merge2(rld_t *e0, rld_t *e1, int n_threads, int64_t offset, int64_t shift)
{
	uint64_t i, k, n = e0->mcnt[0] + e1->mcnt[0], *bits;
	rlditr_t itr, itr0, itr1;
	rld_t *e;
	int last;

	// compute the gap array
	bits = compute_bits(e0, e1, n_threads, offset, shift);
	if (bits == 0) return 0;
	free(e0->frame); free(e1->frame); // deallocate the rank indexes of e0 and e1; they are not needed any more
	e0->frame = e1->frame = 0;
	// initialize the FM-index to be returned, and all the three iterators
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr, 0);
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(e1, &itr1, 0);
	// merge BWTs
	for (i = k = 1, last = bits[0]&1; i < n; ++i) {
		int c = bits[i>>6]>>(i&0x3f)&1;
		if (c != last) {
			if (last == 0) rld_dec_enc(e, &itr, e0, &itr0, k);
			else rld_dec_enc(e, &itr, e1, &itr1, k);
			last = c; k = 1;
		} else ++k;
	}
	if (k) {
		if (last == 0) rld_dec_enc(e, &itr, e0, &itr0, k);
		else rld_dec_enc(e, &itr, e1, &itr1, k);
	}
	// finalize the merge
	assert(itr0.l == 0 && itr1.l == 0); // both e0 and e1 stream should be finished
	free(bits);
	rld_destroy(e0); rld_destroy(e1);
	rld_enc_finish(e, &itr);
	return e;
}


/***********************
 *** Build sub-index ***
 ***********************/

#include <unistd.h>
#include <zlib.h>
#include "bprope6.h"
#include "ropebwt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static unsigned char nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

typedef struct {
	int i, n_threads;
	int64_t mem;
	FILE *fpr;
	const char *prefix;
} worker1_t;

static void dump_rope(const char *prefix, int id, bprope6_t *rope)
{
	FILE *fp;
	int n;
	const uint8_t *s;
	bpriter_t *itr;
	if (prefix) {
		char *fn;
		fn = calloc(strlen(prefix) + 6, 1);
		if (id >= 0) sprintf(fn, "%s.%.4d", prefix, id);
		else strcpy(fn, prefix);
		fp = fopen(fn, "wb");
		free(fn);
	} else fp = stdout;
	itr = bpr_iter_init(rope);
	fwrite("RL6\2", 1, 4, fp);
	while ((s = bpr_iter_next(itr, &n)) != 0) fwrite(s, 1, n, fp);
	fclose(fp);
	free(itr);
	bpr_destroy(rope);
}

static void *worker1(void *data)
{
	worker1_t *w = (worker1_t*)data;
	bprope6_t *rope;
	int len, max_len = 0, n_batches = 0;
	uint8_t *s = 0;
	rope = bpr_init(64, 512);
	while (fread(&len, sizeof(int), 1, w->fpr) == 1) {
		if (len == -1) {
			dump_rope(w->prefix, n_batches * w->n_threads + w->i, rope);
			rope = bpr_init(64, 512);
			++n_batches;
			w->mem = 0;
			continue;
		}
		if (len > max_len) {
			max_len = len;
			kroundup32(max_len);
			s = realloc(s, max_len);
		}
		fread(s, 1, len, w->fpr);
		bpr_insert_string(rope, len, s);
		w->mem = bpr_mem(rope);
	}
	free(s);
	fclose(w->fpr);
	dump_rope(w->prefix, n_batches * w->n_threads + w->i, rope);
	return 0;
}

rld_t *rld_build(const char *fn, int n_threads, int flag, long max_mem, const char *prefix)
{
	gzFile fp;
	FILE **fpw;
	kseq_t *ks;
	int i, j, which = 0, n_batches = 0;
	pthread_t *tid;
	worker1_t *w;
	rld_t *e;
	char *tmpfn;

	w = calloc(n_threads, sizeof(worker1_t));
	tid = calloc(n_threads, sizeof(pthread_t));
	fpw = calloc(n_threads, sizeof(void*));
	for (i = 0; i < n_threads; ++i) {
		int fd[2];
		pipe(fd);
		fpw[i] = fdopen(fd[1], "wb");
		w[i].fpr = fdopen(fd[0], "rb");
		w[i].prefix = prefix;
		w[i].i = i;
		w[i].n_threads = n_threads;
	}
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, worker1, &w[i]);

	fp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int l = ks->seq.l;
		uint8_t *s = (uint8_t*)ks->seq.s;
		int64_t mem;
		// compute encoded string
		for (i = 0; i < l; ++i)
			s[i] = s[i] < 128? nt6_table[s[i]] : 5;
		if ((flag & RB_F_ODD) && (l&1) == 0) { // then check reverse complement
			for (i = 0; i < l>>1; ++i) // is the reverse complement is identical to itself?
				if (s[i] + s[l-1-i] != 5) break;
			if (i == l>>1) --l; // if so, trim 1bp from the end
		}
		// insert
		if (flag & RB_F_FOR) {
			fwrite(&l, sizeof(int), 1, fpw[which]);
			fwrite(s, 1, l, fpw[which]);
			if (++which == n_threads) which = 0;
		}
		if (flag & RB_F_REV) {
			for (i = 0; i < l>>1; ++i) {
				int tmp = s[l-1-i];
				tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
				s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
				s[i] = tmp;
			}
			if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
			fwrite(&l, sizeof(int), 1, fpw[which]);
			fwrite(s, 1, l, fpw[which]);
			if (++which == n_threads) which = 0;
		}
		if (which == 0 && n_batches == 0) {
			for (i = 0, mem = 0; i < n_threads; ++i) mem += w[i].mem;
			if (mem >= max_mem) {
				l = -1;
				for (i = 0; i < n_threads; ++i) fwrite(&l, sizeof(int), 1, fpw[i]);
				++n_batches;
			}
		}
	}
	++n_batches;
	kseq_destroy(ks);
	gzclose(fp);

	for (i = 0; i < n_threads; ++i) fclose(fpw[i]);
	free(fpw);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	free(tid); free(w);

	tmpfn = calloc(strlen(prefix) + 6, 1);
	sprintf(tmpfn, "%s.%.4d", prefix, 0);
	if ((e = rld_restore(tmpfn)) != 0) unlink(tmpfn);
	else return 0;
	if (n_threads > 1) { // merge FM-index
		rld_t *e1, *e0;
		int64_t offset = 0;
		for (j = 0; j < n_batches; ++j) {
			if (j > 0) offset = e->mcnt[1];
			for (i = 0; i < n_threads; ++i) {
				if (j == 0 && i == 0) continue;
				sprintf(tmpfn, "%s.%.4d", prefix, j * n_threads + i);
				e0 = e;
				e1 = rld_restore(tmpfn);
				if (!(e = rld_merge2(e0, e1, n_threads, offset, i))) {
					sprintf(tmpfn, "%s.%.4d", prefix, j * n_threads + i - 1);
					rld_dump(e0, tmpfn);
					free(tmpfn);
					return 0;
				} else unlink(tmpfn);
			}
		}
	}
	free(tmpfn);
	return e;
}

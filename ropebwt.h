#ifndef ROPEBWT_H
#define ROPEBWT_H

#define RB_F_FOR    0x1
#define RB_F_REV    0x2
#define RB_F_ODD    0x4

struct __rld_t;

#define rld_merge(e0, e1, n_threads) rld_merge2((e0), (e1), (n_threads), 0, 0)

#ifdef __cplusplus
extern "C" {
#endif

	struct __rld_t *rld_merge2(struct __rld_t *e0, struct __rld_t *e1, int n_threads, int64_t offset, int64_t shift);
	struct __rld_t *rld_build(const char *fn, int n_threads, int flag, long max_mem, const char *prefix);

#ifdef __cplusplus
}
#endif

#endif

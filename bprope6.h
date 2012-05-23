#ifndef BPROPE6_H
#define BPROPE6_H

#include <stdint.h>

struct bprope6_s;
typedef struct bprope6_s bprope6_t;

struct bpriter_s;
typedef struct bpriter_s bpriter_t;

#ifdef __cplusplus
extern "C" {
#endif

	bprope6_t *bpr_init(int max_nodes, int max_runs);
	void bpr_destroy(bprope6_t *rope);
	int64_t bpr_insert_symbol(bprope6_t *rope, int a, int64_t x);
	void bpr_insert_string(bprope6_t *rope, int l, uint8_t *str);
	void bpr_print(const bprope6_t *rope);

	bpriter_t *bpr_iter_init(const bprope6_t *rope); // to free, simply call free()
	const uint8_t *bpr_iter_next(bpriter_t *iter, int *n, int *l);

#ifdef __cplusplus
}
#endif

#endif

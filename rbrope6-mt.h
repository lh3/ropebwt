#ifndef RBROPE6_MT_H
#define RBROPE6_MT_H

#include <stdint.h>

struct rbmope6_s;
typedef struct rbmope6_s rbmope6_t;

struct rbmiter_s;
typedef struct rbmiter_s rbmiter_t;

#ifdef __cplusplus
extern "C" {
#endif

	rbmope6_t *rbm_init(int max_runs);
	void rbm_destroy(rbmope6_t *rope);
	uint64_t rbm_insert_symbol(rbmope6_t *rope, int a, uint64_t x);
	void rbm_insert_string(rbmope6_t *rope, int l, uint8_t *str);
	void rbm_get_counts(rbmope6_t *rope, uint64_t c[6]);
	void rbm_print(const rbmope6_t *rope);

	rbmiter_t *rbm_iter_init(const rbmope6_t *rope); // to free, simply call free()
	const uint8_t *rbm_iter_next(rbmiter_t *iter, int *n);

#ifdef __cplusplus
}
#endif

#endif

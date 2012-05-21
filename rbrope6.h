#ifndef RBROPE6_H
#define RBROPE6_H

#include <stdint.h>

struct rbrope6_s;
typedef struct rbrope6_s rbrope6_t;

struct rbriter_s;
typedef struct rbriter_s rbriter_t;

#ifdef __cplusplus
extern "C" {
#endif

	rbrope6_t *rbr_init(void);
	void rbr_destroy(rbrope6_t *rope);
	uint64_t rbr_insert_symbol(rbrope6_t *rope, int a, uint64_t x);
	void rbr_insert_string(rbrope6_t *rope, int l, uint8_t *str);
	void rbr_get_counts(rbrope6_t *rope, uint64_t c[6]);

	rbriter_t *rbr_iter_init(const rbrope6_t *rope); // to free, simply call free()
	const uint8_t *rbr_iter_next(rbriter_t *iter, int *n, int *l);

#ifdef __cplusplus
}
#endif

#endif

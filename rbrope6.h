#ifndef RBROPE6_H
#define RBROPE6_H

#include <stdint.h>

struct rbrope6_s;
typedef struct rbrope_s rbrope6_t;

#ifdef __cplusplus
extern "C" {
#endif

	rbrope6_t *rbr_init(void);
	void rbr_destroy(rbrope6_t *rope);
	uint64_t rbr_insert(rbrope6_t *rope, int a, uint64_t x);

#ifdef __cplusplus
}
#endif

#endif

#ifndef BCR_H
#define BCR_H

#include <stdint.h>

struct bcr_s;
typedef struct bcr_s bcr_t;

struct bcritr_s;
typedef struct bcritr_s bcritr_t;

extern int bcr_verbose;

#ifdef __cplusplus
extern "C" {
#endif

	bcr_t *bcr_init(int is_thr, const char *tmpfn);
	void bcr_destroy(bcr_t *b);
	void bcr_append(bcr_t *b, int len, uint8_t *seq);
	void bcr_build(bcr_t *b);

	bcritr_t *bcr_itr_init(const bcr_t *b);
	const uint8_t *bcr_itr_next(bcritr_t *itr, int *l);

#ifdef __cplusplus
}
#endif

#endif

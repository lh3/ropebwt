#ifndef IBWT_H
#define IBWT_H

struct ibwt_s;
typedef struct ibwt_s ibwt_t;

#ifdef __cplusplus
extern "C" {
#endif

	ibwt_t *ib_init();
	void ib_destroy(ibwt_t *ib);
	uint64_t ib_push(ibwt_t *ib, int a, uint64_t x);

#ifdef __cplusplus
}
#endif

#endif

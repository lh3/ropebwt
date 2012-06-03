#ifndef AC_KSORT_H
#define AC_KSORT_H

#include <stdlib.h>
#include <string.h>

typedef struct {
	void *left, *right;
	int depth;
} ks_isort_stack_t;

#define KSORT_SWAP(type_t, a, b) { register type_t t=(a); (a)=(b); (b)=t; }

#define KSORT_INIT(name, type_t, __sort_lt)								\
	inline void __ks_insertsort_##name(type_t *s, type_t *t)			\
	{																	\
		type_t *i, *j, swap_tmp;										\
		for (i = s + 1; i < t; ++i)										\
			for (j = i; j > s && __sort_lt(*j, *(j-1)); --j) {			\
				swap_tmp = *j; *j = *(j-1); *(j-1) = swap_tmp;			\
			}															\
	}																	\
	void ks_combsort_##name(size_t n, type_t a[])						\
	{																	\
		const double shrink_factor = 1.2473309501039786540366528676643; \
		int do_swap;													\
		size_t gap = n;													\
		type_t tmp, *i, *j;												\
		do {															\
			if (gap > 2) {												\
				gap = (size_t)(gap / shrink_factor);					\
				if (gap == 9 || gap == 10) gap = 11;					\
			}															\
			do_swap = 0;												\
			for (i = a; i < a + n - gap; ++i) {							\
				j = i + gap;											\
				if (__sort_lt(*j, *i)) {								\
					tmp = *i; *i = *j; *j = tmp;						\
					do_swap = 1;										\
				}														\
			}															\
		} while (do_swap || gap > 2);									\
		if (gap != 1) __ks_insertsort_##name(a, a + n);					\
	}																	\
	void ks_introsort_##name(size_t n, type_t a[])						\
	{																	\
		int d;															\
		ks_isort_stack_t *top, *stack;									\
		type_t rp, swap_tmp;											\
		type_t *s, *t, *i, *j, *k;										\
																		\
		if (n < 1) return;												\
		else if (n == 2) {												\
			if (__sort_lt(a[1], a[0])) { swap_tmp = a[0]; a[0] = a[1]; a[1] = swap_tmp; } \
			return;														\
		}																\
		for (d = 2; 1ul<<d < n; ++d);									\
		stack = (ks_isort_stack_t*)malloc(sizeof(ks_isort_stack_t) * ((sizeof(size_t)*d)+2)); \
		top = stack; s = a; t = a + (n-1); d <<= 1;						\
		while (1) {														\
			if (s < t) {												\
				if (--d == 0) {											\
					ks_combsort_##name(t - s + 1, s);					\
					t = s;												\
					continue;											\
				}														\
				i = s; j = t; k = i + ((j-i)>>1) + 1;					\
				if (__sort_lt(*k, *i)) {								\
					if (__sort_lt(*k, *j)) k = j;						\
				} else k = __sort_lt(*j, *i)? i : j;					\
				rp = *k;												\
				if (k != t) { swap_tmp = *k; *k = *t; *t = swap_tmp; }	\
				for (;;) {												\
					do ++i; while (__sort_lt(*i, rp));					\
					do --j; while (i <= j && __sort_lt(rp, *j));		\
					if (j <= i) break;									\
					swap_tmp = *i; *i = *j; *j = swap_tmp;				\
				}														\
				swap_tmp = *i; *i = *t; *t = swap_tmp;					\
				if (i-s > t-i) {										\
					if (i-s > 16) { top->left = s; top->right = i-1; top->depth = d; ++top; } \
					s = t-i > 16? i+1 : t;								\
				} else {												\
					if (t-i > 16) { top->left = i+1; top->right = t; top->depth = d; ++top; } \
					t = i-s > 16? i-1 : s;								\
				}														\
			} else {													\
				if (top == stack) {										\
					free(stack);										\
					__ks_insertsort_##name(a, a+n);						\
					return;												\
				} else { --top; s = (type_t*)top->left; t = (type_t*)top->right; d = top->depth; } \
			}															\
		}																\
	}

#define ks_introsort(name, n, a) ks_introsort_##name(n, a)
#define ks_combsort(name, n, a) ks_combsort_##name(n, a)

#endif
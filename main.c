#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bprope6.h"
#include "rbrope6.h"
#include "rbrope6-mt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void seq_char2nt6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

int main(int argc, char *argv[])
{
	rbmope6_t *rbm = 0;
	rbrope6_t *rbr = 0;
	bprope6_t *bpr = 0;
	gzFile fp;
	kseq_t *ks;
	int c, i, j, for_only = 0, print_rope = 0, max_runs = 512, max_nodes = 64, n_threads = 1, batch_size = 0x10000;
	const uint8_t *s;

	while ((c = getopt(argc, argv, "Tfr:n:t:b:")) >= 0)
		if (c == 'f') for_only = 1;
		else if (c == 'T') print_rope = 1;
		else if (c == 'r') max_runs = atoi(optarg);
		else if (c == 'n') max_nodes= atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'b') batch_size = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: ropebwt [-fT] [-r maxRuns=%d] [-n maxNodes=%d] [-b batchSize=%d] <in.fq.gz>\n", max_runs, max_nodes, batch_size);
		return 1;
	}

	if (n_threads > 1) rbm = rbm_init(n_threads, batch_size, max_runs);
	else if (max_nodes <= 2) rbr = rbr_init(max_runs);
	else bpr = bpr_init(max_nodes, max_runs);
	fp = gzopen(argv[optind], "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		uint8_t *s = (uint8_t*)ks->seq.s;
		seq_char2nt6(ks->seq.l, s);
		if (rbm) rbm_insert_string(rbm, ks->seq.l, s);
		if (rbr) rbr_insert_string(rbr, ks->seq.l, s);
		if (bpr) bpr_insert_string(bpr, ks->seq.l, s);
		if (!for_only) {
			seq_revcomp6(ks->seq.l, s);
			if (rbm) rbm_insert_string(rbm, ks->seq.l, s);
			if (rbr) rbr_insert_string(rbr, ks->seq.l, s);
			if (bpr) bpr_insert_string(bpr, ks->seq.l, s);
		}
	}
	kseq_destroy(ks);
	gzclose(fp);

	if (rbm) {
		rbmiter_t *iter;
		rbm_update_bcr(rbm);
		iter = rbm_iter_init(rbm);
		while ((s = rbm_iter_next(iter, &c)) != 0)
			for (i = 0; i < c; ++i)
				for (j = 0; j < s[i]>>3; ++j)
					putchar("$ACGTN"[s[i]&7]);
		putchar('\n');
		if (print_rope) rbm_print(rbm);
		free(iter);
		rbm_destroy(rbm);
	}
	if (rbr) {
		rbriter_t *iter;
		iter = rbr_iter_init(rbr);
		while ((s = rbr_iter_next(iter, &c)) != 0)
			for (i = 0; i < c; ++i)
				for (j = 0; j < s[i]>>3; ++j)
					putchar("$ACGTN"[s[i]&7]);
		putchar('\n');
		if (print_rope) rbr_print(rbr);
		free(iter);
		rbr_destroy(rbr);
	}
	if (bpr) {
		bpriter_t *iter = bpr_iter_init(bpr);
		while ((s = bpr_iter_next(iter, &c)) != 0)
			for (i = 0; i < c; ++i)
				for (j = 0; j < s[i]>>3; ++j)
					putchar("$ACGTN"[s[i]&7]);
		putchar('\n');
		if (print_rope) bpr_print(bpr);
		free(iter);
		bpr_destroy(bpr);
	}
	return 0;
}

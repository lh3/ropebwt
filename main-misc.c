#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bprope6.h"
#include "rbrope6.h"
#include "bcr.h"
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

enum algo_e { BPR, RBR, BCR };

int main(int argc, char *argv[])
{
	rbrope6_t *rbr = 0;
	bprope6_t *bpr = 0;
	bcr_t *bcr = 0;
	gzFile fp;
	kseq_t *ks;
	enum algo_e algo = BPR;
	int c, i, j, for_only = 0, print_rope = 0, max_runs = 512, max_nodes = 64, n_threads = 1;
	int64_t cnt = 0;
	const uint8_t *s;

	while ((c = getopt(argc, argv, "Tfr:n:t:a:")) >= 0)
		if (c == 'a') {
			if (strcmp(optarg, "bpr") == 0) algo = BPR;
			else if (strcmp(optarg, "rbr") == 0) algo = RBR;
			else if (strcmp(optarg, "bcr") == 0) algo = BCR;
			else fprintf(stderr, "[W::%s] available algorithms: bpr, rbr or bcr; default to bpr\n", __func__);
		} else if (c == 'f') for_only = 1;
		else if (c == 'T') print_rope = 1;
		else if (c == 'r') max_runs = atoi(optarg);
		else if (c == 'n') max_nodes= atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: ropebwt [-fT] [-r maxRuns=%d] [-n maxNodes=%d] <in.fq.gz>\n", max_runs, max_nodes);
		return 1;
	}

	if (algo == BPR) bpr = bpr_init(max_nodes, max_runs);
	else if (algo == RBR) rbr = rbr_init(max_runs);
	else if (algo == BCR) bcr = bcr_init(n_threads > 1);
	fp = gzopen(argv[optind], "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		uint8_t *s = (uint8_t*)ks->seq.s;
		++cnt;
		seq_char2nt6(ks->seq.l, s);
		if (rbr) rbr_insert_string(rbr, ks->seq.l, s);
		if (bpr) bpr_insert_string(bpr, ks->seq.l, s);
		if (bcr) bcr_append(bcr, ks->seq.l, s);
		if (!for_only) {
			++cnt;
			seq_revcomp6(ks->seq.l, s);
			if (rbr) rbr_insert_string(rbr, ks->seq.l, s);
			if (bpr) bpr_insert_string(bpr, ks->seq.l, s);
			if (bcr) bcr_append(bcr, ks->seq.l, s);
		}
	}
	fprintf(stderr, "[M::%s] processed %ld sequences\n", __func__, (long)cnt);
	kseq_destroy(ks);
	gzclose(fp);

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
	if (bcr) {
		bcritr_t *itr;
		bcr_build(bcr);
		itr = bcr_itr_init(bcr);
		while ((s = bcr_itr_next(itr, &c)) != 0)
			for (i = 0; i < c; ++i)
				for (j = 0; j < s[i]>>3; ++j)
					putchar("$ACGTN"[s[i]&7]);
		putchar('\n');
		free(itr);
		bcr_destroy(bcr);
	}
	return 0;
}

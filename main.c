#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "rbrope6.h"
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
	rbrope6_t *rope;
	gzFile fp;
	kseq_t *ks;
	int n, l, c, for_only = 0, max_runs = 256;
	const uint8_t *s;
	rbriter_t *iter;

	while ((c = getopt(argc, argv, "fm:")) >= 0)
		if (c == 'f') for_only = 1;
		else if (c == 'm') max_runs = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: ropebwt [-f] [-m maxRuns=%d] <in.fq.gz>\n", max_runs);
		return 1;
	}

	rope = rbr_init(max_runs);
	fp = gzopen(argv[optind], "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		uint8_t *s = (uint8_t*)ks->seq.s;
		seq_char2nt6(ks->seq.l, s);
		rbr_insert_string(rope, ks->seq.l, s);
		if (!for_only) {
			seq_revcomp6(ks->seq.l, s);
			rbr_insert_string(rope, ks->seq.l, s);
		}
	}
	kseq_destroy(ks);
	gzclose(fp);

	iter = rbr_iter_init(rope);
	while ((s = rbr_iter_next(iter, &n, &l)) != 0) {
		int i, j;
//		printf("%d\t%d\t", n, l);
		for (i = 0; i < n; ++i)
			for (j = 0; j < s[i]>>3; ++j)
				putchar("$ACGTN"[s[i]&7]);
//		putchar('\n');
	}
	putchar('\n');
	rbr_destroy(rope);
	return 0;
}

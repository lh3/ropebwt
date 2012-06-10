#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bprope6.h"
#include "rbrope6.h"
#include "bcr.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

enum algo_e { BPR, RBR, BCR };

#define FLAG_FOR 0x1
#define FLAG_REV 0x2
#define FLAG_ODD 0x4
#define FLAG_BIN 0x8
#define FLAG_TREE 0x10
#define FLAG_THR 0x20
#define FLAG_NON 0x40

int main(int argc, char *argv[])
{
	rbrope6_t *rbr = 0;
	bprope6_t *bpr = 0;
	bcr_t *bcr = 0;
	gzFile fp;
	FILE *out = stdout;
	char *tmpfn = 0;
	kseq_t *ks;
	enum algo_e algo = BPR;
	int c, i, max_runs = 512, max_nodes = 64;
	int flag = FLAG_FOR | FLAG_REV | FLAG_ODD;

	while ((c = getopt(argc, argv, "TFRObNo:r:n:ta:f:v:")) >= 0)
		if (c == 'a') {
			if (strcmp(optarg, "bpr") == 0) algo = BPR;
			else if (strcmp(optarg, "rbr") == 0) algo = RBR;
			else if (strcmp(optarg, "bcr") == 0) algo = BCR;
			else fprintf(stderr, "[W::%s] available algorithms: bpr, rbr or bcr; default to bpr\n", __func__);
		} else if (c == 'o') out = fopen(optarg, "wb");
		else if (c == 'F') flag &= ~FLAG_FOR;
		else if (c == 'R') flag &= ~FLAG_REV;
		else if (c == 'O') flag &= ~FLAG_ODD;
		else if (c == 'T') flag |= FLAG_TREE;
		else if (c == 'b') flag |= FLAG_BIN;
		else if (c == 'N') flag |= FLAG_NON;
		else if (c == 't') flag |= FLAG_THR;
		else if (c == 'r') max_runs = atoi(optarg);
		else if (c == 'n') max_nodes= atoi(optarg);
		else if (c == 'f') tmpfn = optarg;
		else if (c == 'v') bcr_verbose = atoi(optarg);

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   ropebwt [options] <in.fq.gz>\n\n");
		fprintf(stderr, "Options: -a STR     algorithm: bpr, rbr or bcr [bpr]\n");
		fprintf(stderr, "         -r INT     max number of runs in leaves (bpr and rbr only) [%d]\n", max_runs);
		fprintf(stderr, "         -n INT     max number children per internal node (bpr only) [%d]\n", max_nodes);
		fprintf(stderr, "         -o FILE    output file [stdout]\n");
		fprintf(stderr, "         -f FILE    temporary sequence file name (bcr only) [null]\n");
		fprintf(stderr, "         -v INT     verbose level (bcr only) [%d]\n", bcr_verbose);
		fprintf(stderr, "         -b         binary output (5+3 runs starting after 4 bytes)\n");
		fprintf(stderr, "         -t         enable threading (bcr only)\n");
		fprintf(stderr, "         -F         skip forward strand\n");
		fprintf(stderr, "         -R         skip reverse strand\n");
		fprintf(stderr, "         -N         discard reads containing ambiguous bases\n");
		fprintf(stderr, "         -O         suppress end trimming when forward==reverse\n");
		fprintf(stderr, "         -T         print the tree stdout (bpr and rbr only)\n\n");
		return 1;
	}

	if (algo == BCR) {
		bcr = bcr_init(flag&FLAG_THR, tmpfn);
		if (!(flag&FLAG_NON)) fprintf(stderr, "Warning: With bcr, an ambiguous base will be converted to a random base\n");
	} else if (algo == BPR) bpr = bpr_init(max_nodes, max_runs);
	else if (algo == RBR) rbr = rbr_init(max_runs);
	fp = gzopen(argv[optind], "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int l = ks->seq.l;
		uint8_t *s = (uint8_t*)ks->seq.s;
		for (i = 0; i < l; ++i)
			s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
		if (flag & FLAG_NON) {
			for (i = 0; i < l; ++i) if (s[i] == 5) break;
			if (i != l) continue;
		} else if (algo == BCR)
			for (i = 0; i < l; ++i)
				if (s[i] == 5) s[i] = (lrand48()&3) + 1;
		if ((flag & FLAG_ODD) && (l&1) == 0) { // then check reverse complement
			for (i = 0; i < l>>1; ++i) // is the reverse complement is identical to itself?
				if (s[i] + s[l-1-i] != 5) break;
			if (i == l>>1) --l; // if so, trim 1bp from the end
		}
		if (flag & FLAG_FOR) {
			if (rbr) rbr_insert_string(rbr, ks->seq.l, s);
			if (bpr) bpr_insert_string(bpr, ks->seq.l, s);
			if (bcr) bcr_append(bcr, ks->seq.l, s);
		}
		if (flag & FLAG_REV) {
			for (i = 0; i < l>>1; ++i) {
				int tmp = s[l-1-i];
				tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
				s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
				s[i] = tmp;
			}
			if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
			if (rbr) rbr_insert_string(rbr, ks->seq.l, s);
			if (bpr) bpr_insert_string(bpr, ks->seq.l, s);
			if (bcr) bcr_append(bcr, ks->seq.l, s);
		}
	}
	kseq_destroy(ks);
	gzclose(fp);

#define print_bwt(itr_t, itr_set, itr_next_f, is_bin, fp) do { \
		itr_t *itr; \
		const uint8_t *s; \
		int i, j, l; \
		itr = (itr_set); \
		if (is_bin) { \
			fwrite("RLE\6", 4, 1, fp); \
			while ((s = itr_next_f(itr, &l)) != 0) \
				fwrite(s, 1, l, fp); \
		} else { \
			while ((s = itr_next_f(itr, &l)) != 0) \
				for (i = 0; i < l; ++i) \
					for (j = 0; j < s[i]>>3; ++j) \
						fputc("$ACGTN"[s[i]&7], fp); \
			fputc('\n', fp); \
		} \
		free(itr); \
	} while (0)

	if (rbr) {
		print_bwt(rbriter_t, rbr_iter_init(rbr), rbr_iter_next, flag&FLAG_BIN, out);
		if (flag&FLAG_TREE) rbr_print(rbr);
		rbr_destroy(rbr);
	}
	if (bpr) {
		print_bwt(bpriter_t, bpr_iter_init(bpr), bpr_iter_next, flag&FLAG_BIN, out);
		if (flag&FLAG_TREE) bpr_print(bpr);
		bpr_destroy(bpr);
	}
	if (bcr) {
		bcr_build(bcr);
		print_bwt(bcritr_t, bcr_itr_init(bcr), bcr_itr_next, flag&FLAG_BIN, out);
		bcr_destroy(bcr);
	}
	fclose(out);
	return 0;
}

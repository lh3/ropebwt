#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include "rld.h"
#include "bpr-mt.h"

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
static int64_t get_max_mem()
{
	int64_t mem;
	size_t len;
	int mib[4];
	mib[0] = CTL_HW; mib[1] = HW_MEMSIZE;
	len = 8;
	sysctl(mib, 2, &mem, &len, 0, 0);
	return mem;
}
#elif defined(__linux__)
#include <sys/time.h>
#include <sys/resource.h>
int64_t get_max_mem()
{
	FILE *fp;
	char buffer[64];
	int64_t mem = INT64_MAX;
	struct rlimit r;
	if ((fp = fopen("/proc/meminfo", "r")) != NULL) {
		while (fscanf(fp, "%s", buffer) > 0) {
			if (strstr(buffer, "MemTotal") == buffer) {
				long tmp;
				fscanf(fp, "%lu", &tmp);
				mem = 1024LL * tmp;
			}
		}
		fclose(fp);
	}
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
	if (mem > r.rlim_max) mem = r.rlim_max;
	return mem;
}
#endif

int main(int argc, char *argv[])
{
	int c, n_threads = 1, is_stdout = 0, flag = RB_F_FOR | RB_F_REV | RB_F_ODD;
	int64_t max_mem = INT64_MAX;
	rld_t *e;
#if defined(__linux) || defined(__APPLE__)
	max_mem = get_max_mem();
#endif
	while ((c = getopt(argc, argv, "FROot:m:")) >= 0)
		if (c == 'F') flag &= ~RB_F_FOR;
		else if (c == 'R') flag &= ~RB_F_REV;
		else if (c == 'O') flag &= ~RB_F_ODD;
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'm') {
			char *p = optarg;
			max_mem = strtol(p, &p, 10);
			if (*p == 'k' || *p == 'K') max_mem <<= 10;
			else if (*p == 'm' || *p == 'M') max_mem <<= 20;
			else if (*p == 'g' || *p == 'G') max_mem <<= 30;
		} else if (c == 'o') is_stdout = 1;
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: ropebwt [-FROo] [-t nThreads=1] [-m maxMem] <in.fq.gz> <out.fmd>\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Maximum memory is set as %ld\n", __func__, (long)max_mem);
	if ((e = rld_build(argv[optind], n_threads, flag, max_mem, argv[optind+1])) != 0) {
		if (!is_stdout) {
			char *fn;
			fn = calloc(strlen(argv[optind+1]) + 5, 1);
			sprintf(fn, "%s.fmd", argv[optind+1]);
			rld_dump(e, fn);
			free(fn);
		} else rld_dump(e, "-");
		rld_destroy(e);
	}
	return 0;
}

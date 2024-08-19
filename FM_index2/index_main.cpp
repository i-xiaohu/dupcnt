//
// Created by ixiaohu on 2024/8/19.
//

#include "bntseq.h"
#include "FMI_search.h"

int bwa_idx_build(const char *fa, const char *prefix)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	clock_t t;
	int64_t l_pac;

	{ // nucleotide indexing
		// BWA-MEM2 use the original bwa index code to manage packed nucleotides, annotations and ambiguous holes.
		gzFile fp = xzopen(fa, "r");
		t = clock();
		fprintf(stderr, "[bwa_index] Pack FASTA... ");
		bns_fasta2bntseq(fp, prefix, 1);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
		FMI_search *fmi = new FMI_search(prefix);
		fmi->build_index();
		delete fmi;
	}
	return 0;
}

int main(int argc, char *argv[]) // the "index" command
{
	int c;
	char *prefix = 0;
	while ((c = getopt(argc, argv, "p:")) >= 0) {
		if (c == 'p') prefix = optarg;
		else return 1;
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: fmidx2 [-p prefix] <in.fasta>\n");
		return 1;
	}
	if (prefix == 0) prefix = argv[optind];
	bwa_idx_build(argv[optind], prefix);
	return 0;
}


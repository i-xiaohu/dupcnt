//
// Created by ixiaohu on 2024/8/7.
//

#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include "dupcnt_core.h"

int usage() {
	fprintf(stderr, "Usage: dupcnt [options] <Sample1> [Sample2] ...\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -t <INT> Threads [16]\n");
	fprintf(stderr ,"  -i <STR> Index prefix\n");
	fprintf(stderr, "  -k <INT> Batch size\n");
	fprintf(stderr, "  -c <STR> Trie memory usage upperbound [100G]\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		return usage();
	}
	auto *opt = new Option();
	int c;
	while ((c = getopt(argc, argv, "t:i:k:c:")) >= 0) {
		if (c == 't') {
			opt->n_threads = atoi(optarg);
		} else if (c == 'i') {
			opt->index_prefix = optarg;
		} else if (c == 'k') {
			opt->batch_size = atoi(optarg);
		} else if (c == 'c') {
			opt->mem_cap = 0;
			for (int i = 0; optarg[i]; i++) {
				if (optarg[i] >= '0' and optarg[i] <= '9') {
					opt->mem_cap *= 10L;
					opt->mem_cap += optarg[i] - '0';
				} else if (optarg[i] == 'M' or optarg[i] == 'm') {
					opt->mem_cap *= 1024L * 1024L;
				} else if (optarg[i] == 'G' or optarg[i] == 'g') {
					opt->mem_cap *= 1024L * 1024L * 1024L;
				}
			}
		}
	}

	if (opt->index_prefix == nullptr) {
		fprintf(stderr, "Please provide index file\n");
		return 1;
	}

	fprintf(stderr, "%d samples to process\n", argc - optind);
	process(opt, argc - optind, argv + optind);

	return 0;
}
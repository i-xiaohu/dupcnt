//
// Created by ixiaohu on 2024/8/7.
//

#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include "dupcnt_core.h"

int usage() {
	fprintf(stderr, "Usage: dupcnt -t <threads> -i <index> <Sample1> [Sample2] ...\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		return usage();
	}
	int c, threads_n = 16;
	char *index_prefix = nullptr;
	while ((c = getopt(argc, argv, "t:i:")) >= 0) {
		if (c == 't') {
			threads_n = atoi(optarg);
		} else if (c == 'i') {
			index_prefix = optarg;
		}
	}

	if (index_prefix == nullptr) {
		fprintf(stderr, "Please provide index file\n");
		return 1;
	}

	fprintf(stderr, "%d samples to process\n", argc - optind);
	process(threads_n, index_prefix, argc - optind, argv + optind);

	return 0;
}
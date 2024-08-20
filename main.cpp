//
// Created by ixiaohu on 2024/8/7.
//

#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <algorithm>
#include "dupcnt_core.h"

int usage() {
	fprintf(stderr, "Usage: dupcnt [options] <Sample1> [Sample2] ...\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -t <INT> Threads [16]\n");
	fprintf(stderr ,"  -i <STR> Index prefix\n");
	fprintf(stderr, "  -k <INT> Batch size\n");
	fprintf(stderr, "  -a       Output all reads to stdout\n");
	fprintf(stderr, "  -s       Sort reads by frequency before output\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		return usage();
	}
	auto *opt = new Option();
	int c;
	while ((c = getopt(argc, argv, "t:i:k:a:s:h:")) >= 0) {
		if (c == 't') {
			opt->n_threads = atoi(optarg);
		} else if (c == 'i') {
			opt->index_prefix = optarg;
		} else if (c == 'k') {
			opt->batch_size = atoi(optarg);
		} else if (c == 'a') {
			opt->output = true;
		} else if (c == 's') {
			opt->sorted = true;
		} else {
			fprintf(stderr, "Unrecognized option `-%c`\n", c);
			delete opt;
			return usage();
		}
	}

	if (opt->index_prefix == nullptr) {
		fprintf(stderr, "Please provide index file\n");
		return 1;
	}

	fprintf(stderr, "%d sample(s) to process\n", argc - optind);
	process(opt, argc - optind, argv + optind);

	delete opt;

	return 0;
}
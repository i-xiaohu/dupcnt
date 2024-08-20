//
// Created by ixiaohu on 2024/8/7.
//

#ifndef DUPCNT_DEDUP_CORE_H
#define DUPCNT_DEDUP_CORE_H

#include <cstdint>
#include <utility>
#include <vector>
#include <string>

#define AVL_SHIFT  6
#define AVL_BUCKET 4096

struct Option {
	int n_threads = 16;
	int batch_size = 10 * 1000 * 1000; // 10M bases for each thread
	const char *index_prefix = nullptr;
	bool output = false;
	bool sorted = false;
};

void process(const Option *opt, int n_sample, char *files[]);

#endif //DUPCNT_DEDUP_CORE_H

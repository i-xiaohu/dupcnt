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

struct RepRead {
	int occ;
	std::string read;
	bool operator < (const RepRead &r) const {
		return occ > r.occ; // Put less frequent read at the top of heap
	}
};

struct Option {
	int n_threads = 16;
	int batch_size = 10 * 1000 * 1000; // 10M bases for each thread
	long most_rep = 0;
	const char *index_prefix = nullptr;
	size_t mem_cap = 100L * 1024L * 1024L * 1024L; // 100GB
};

void process(const Option *opt, int n_sample, char *files[]);

#endif //DUPCNT_DEDUP_CORE_H

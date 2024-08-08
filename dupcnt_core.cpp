//
// Created by ixiaohu on 2024/8/7.
//

#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <string>
#include <cassert>
#include <sys/time.h>
#include <algorithm>
#include "dupcnt_core.h"

#define MAX_LINE 10240
char inbuf[MAX_LINE];
int read_length_monitor = -1; // To avoid any read of different length

std::vector<std::string> input_reads(gzFile fp, int batch_size) {
	int bases_n = 0;
	std::vector<std::string> ret;
	while (gzgets(fp, inbuf, MAX_LINE)) {
		gzgets(fp, inbuf, MAX_LINE);
		int len = strlen(inbuf);
		inbuf[--len] = '\0';
		ret.emplace_back(std::string(inbuf));
		gzgets(fp, inbuf, MAX_LINE);
		gzgets(fp, inbuf, MAX_LINE);
		if (read_length_monitor == -1) read_length_monitor = len;
		else assert(read_length_monitor == len);
		bases_n += len;
		if (bases_n > batch_size) {
			break;
		}
	}
	return ret;
}

Trie::Trie() {
	nodes.clear();
	nodes.emplace_back(TrNode());
}

void Trie::add_read(int n, const char *s) {
	reads_n += 1;
	int32_t p = 0;
	for (int i = 0; i < n; i++) {
		uint8_t c;
		switch (s[i]) {
			case 'A': c = 0; break;
			case 'C': c = 1; break;
			case 'G': c = 2; break;
			case 'T': c = 3; break;
			default: c = 0; break;
		}
		if (nodes[p].x[c] == 0) { // No child
			nodes[p].x[c] = nodes.size();
			nodes.emplace_back(TrNode());
		}
		p = nodes[p].x[c]; // Move down to the next layer
	}
	// Remember x[0] is initialized as 0 thus can be used to track the occurrence number by +1
	nodes[p].x[0]++;
}

void Trie::stat() {
	fprintf(stderr, "Trie size: %ld\n", nodes.size());
	int32_t p = 0;
	// todo: switch back to vector
	std::vector<int> prev, curr;
	prev.push_back(p);
	for (int i = 0; i < read_length_monitor; i++) {
		curr.clear();
		for (int t : prev) {
			for (int k : nodes[t].x) {
				if (k) {
					curr.push_back(k);
				}
			}
		}
		std::swap(prev, curr);
	}
	unique_n = prev.size(); // Leaves are those unique reads
	int max_occ = 0, max_id = -1;
	for (int leaf : prev) {
		assert(nodes[leaf].x[0] > 0);
		if (nodes[leaf].x[0] > max_occ) {
			max_occ = nodes[leaf].x[0];
			max_id = leaf;
		}
	}
	// Linear scan to restore the most repetitive read
	std::string rep_read;
	int child = max_id;
	for (int i = max_id - 1; i >= 0; i--) {
		for (int j = 0; j < 4; j++) {
			if (nodes[i].x[j] == child) {
				rep_read += "ACGT"[j];
				child = i;
				break;
			}
		}
	}
	assert(rep_read.length() == read_length_monitor);
	std::reverse(rep_read.begin(), rep_read.end());
	fprintf(stderr, "The most repetitive read occurs %d times:\n", max_occ);
	fprintf(stderr, "%s\n", rep_read.c_str());
	fprintf(stderr, "%d unique reads in %d reads\n", unique_n, reads_n);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void process(int n_threads, int n_sample, char *files[]) {
	const int BATCH_SIZE = 10 * 1000 * 1000; // 10M bases for each thread
	int actual_batch_size = n_threads * BATCH_SIZE;
	for (int i = 0; i < n_sample; i++) {
		gzFile fp = gzopen(files[i], "r");
		if (fp == nullptr) {
			fprintf(stderr, "Open `%s` failed\n", files[i]);
			continue;
		}
		std::vector<std::string> all_seqs;
		std::vector<std::string> seqs;
		Trie counter;
		while (true) {
			seqs = input_reads(fp, actual_batch_size);
			if (seqs.empty()) break;
			for (auto &r : seqs) {
				for (auto &c : r) {
					if (c == 'N') c = 'A';
				}
				counter.add_read(r.length(), r.c_str());
				all_seqs.push_back(r);
			}
		}
		double t_start, t_end;
		t_start = realtime();
		counter.stat();
		t_end = realtime();
		fprintf(stderr, "Visiting trie in %.2f seconds\n", t_end - t_start);
	}
}

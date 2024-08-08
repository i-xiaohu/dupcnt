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
#include "bwalib/bwa.h"

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

std::vector<std::string> exact_matching(const bwaidx_t *idx, const std::vector<std::string> &seqs, int32_t *em_counter) {
	std::vector<std::string> ret;
	auto *b = (uint8_t*) malloc(read_length_monitor * sizeof(uint8_t));
	for (const auto &s : seqs) {
		for (int i = 0; i < read_length_monitor; i++) {
			b[i] = nst_nt4_table[(uint8_t)s[i]];
			assert(b[i] < 4);
		}
		bwtintv_t ik, ok[4];
		bwt_set_intv(idx->bwt, b[0], ik);
		for (int i = 1; i < read_length_monitor; i++) {
			bwt_extend(idx->bwt, &ik, ok, 0);
			ik = ok[3 - b[i]];
			if (ik.x[2] == 0) {
				break;
			}
		}
		if (ik.x[2] > 0) {
			// Use the first occurrence in suffix array
			int64_t x = bwt_sa(idx->bwt, ik.x[0]);
			em_counter[x]++;
		} else {
			ret.push_back(s);
		}
	}
	free(b);
	return ret;
}

void process(int n_threads, const char *index_prefix, int n_sample, char *files[]) {
	double t_start, t_end;
	/* FM-index from BWA-MEM; todo: do consider switch to BWA-MEM2 index */
	bwaidx_t *idx = bwa_idx_load_from_shm(index_prefix);
	if (idx == nullptr) {
		if ((idx = bwa_idx_load(index_prefix, BWA_IDX_ALL)) == nullptr) {
			fprintf(stderr, "Load index `%s` failed\n", index_prefix);
			return ;
		}
	}
	int64_t ref_len = idx->bns->l_pac * 2;
	fprintf(stderr, "Reference length: %ld\n", ref_len);
	auto *em_counter = (int32_t*) calloc(ref_len, sizeof(int32_t));

	const int BATCH_SIZE = 10 * 1000 * 1000; // 10M bases for each thread
	int actual_batch_size = n_threads * BATCH_SIZE;
	for (int i = 0; i < n_sample; i++) {
		gzFile fp = gzopen(files[i], "r");
		if (fp == nullptr) {
			fprintf(stderr, "Open `%s` failed\n", files[i]);
			continue;
		}
		double t_input = 0, t_match = 0, t_trie = 0;
		int64_t n_seqs = 0, n_match = 0;
		std::vector<std::string> seqs;
		Trie counter;
		while (true) {
			t_start = realtime();
			seqs = input_reads(fp, actual_batch_size);
			t_end = realtime();
			t_input += t_end - t_start;
			if (seqs.empty()) break;
			// Convert non-ACGT to A
			for (auto &r : seqs) {
				for (auto &c : r) {
					if (c == 'N') c = 'A';
				}
			}
			n_seqs += seqs.size();

			// 1. Identify exactly matched reads
			t_start = realtime();
			auto rest = exact_matching(idx, seqs, em_counter);
			t_end = realtime();
			t_match += t_end - t_start;
			fprintf(stderr, "%.2f %% matched reads\n", 100.0 * (double)(seqs.size() - rest.size()) / seqs.size());
			n_match += seqs.size() - rest.size();

			// 2. Process unmatched reads using trie
			// todo: trie is slower than sorting
			t_start = realtime();
			for (auto &r : rest) {
				counter.add_read(r.length(), r.c_str());
			}
			t_end = realtime();
			t_trie += t_end - t_start;
		}

		t_start = realtime();
		counter.stat();
		t_end = realtime();

		fprintf(stderr, "Input: %.2f s, Match: %.2f s, Trie: %.2f s, Post: %.2f s\n", t_input, t_match, t_trie, t_end - t_start);
		fprintf(stderr, "Exactly matched reads: %ld (%.2f %%)\n", n_match, 100.0 * n_match / n_seqs);

		int max_occ = 0;
		int64_t max_pos = -1;
		for (int64_t j = 0; j < ref_len; j++) {
			if (em_counter[j] > max_occ) {
				max_occ = em_counter[j];
				max_pos = j;
			}
		}
		fprintf(stderr, "Read matched at %ld occurs %d times\n", max_pos, max_occ);
	}

	free(em_counter);
}

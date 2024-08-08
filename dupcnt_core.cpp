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

#include "bwalib/kseq.h"
KSEQ_INIT(gzFile, gzread)

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
		uint8_t c = s[i];
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
	fprintf(stderr, "%ld unique reads in %ld reads\n", unique_n, reads_n);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

static inline void bseq1_destroy(const bseq1_t *b) {
	free(b->name);
	free(b->comment);
	free(b->seq);
	free(b->qual);
}

int exact_matching(const bwaidx_t *idx, int n, bseq1_t *seqs, int32_t *em_counter) {
	std::vector<bseq1_t> ret;
	int not_match_n = 0;
	for (int i = 0; i < n; i++) {
		bseq1_t *b = &seqs[i];
		auto *s = (uint8_t*)seqs[i].seq;
		bwtintv_t ik, ok[4];
		bwt_set_intv(idx->bwt, s[0], ik);
		for (int j = 1; j < b->l_seq; j++) {
			bwt_extend(idx->bwt, &ik, ok, 0);
			ik = ok[3 - s[j]];
			if (ik.x[2] == 0) {
				break;
			}
		}
		if (ik.x[2] > 0) {
			// Use the first occurrence in suffix array
			int64_t x = bwt_sa(idx->bwt, ik.x[0]);
			em_counter[x]++;
			bseq1_destroy(&seqs[i]);
		} else {
			seqs[not_match_n++] = seqs[i];
		}
	}
	return not_match_n;
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
		int total_seqs = 0, n_seqs = 0, n_match = 0;
		kseq_t *ks = kseq_init(fp);
		bseq1_t *seqs;
		Trie counter;
		while ((seqs = bseq_read(actual_batch_size, &n_seqs, ks, nullptr))) {
			total_seqs += n_seqs;
			for (int j = 0; j < n_seqs; j++) {
				bseq1_t *b = &seqs[j];
				if (read_length_monitor == -1) read_length_monitor = b->l_seq;
				if (read_length_monitor != b->l_seq) {
					fprintf(stderr, "Only fix-length reads are supported\n");
					std::abort();
				}
				for (int k = 0; k < b->l_seq; k++) {
					b->seq[k] = nst_nt4_table[(uint8_t)b->seq[k]];
					// Convert non-ACGT to A
					if (b->seq[k] > 3) b->seq[k] = 0;
				}
			}

			// 1. Identify exactly matched reads
			t_start = realtime();
			int n_rest = exact_matching(idx, n_seqs, seqs, em_counter);
			t_end = realtime();
			t_match += t_end - t_start;
			fprintf(stderr, "%.2f %% matched reads\n", 100.0 * (double)(n_seqs - n_rest) / n_seqs);
			n_match += n_seqs - n_rest;

			// 2. Process unmatched reads using trie
			// todo: trie is slower than sorting
			t_start = realtime();
			for (int j = 0; j < n_rest; j++) {
				bseq1_t *b = &seqs[j];
				counter.add_read(b->l_seq, b->seq);
				bseq1_destroy(b);
			}
			t_end = realtime();
			t_trie += t_end - t_start;
			free(seqs);
		}

		t_start = realtime();
		counter.stat();
		t_end = realtime();

		fprintf(stderr, "Input: %.2f s, Match: %.2f s, Trie: %.2f s, Post: %.2f s\n", t_input, t_match, t_trie, t_end - t_start);
		fprintf(stderr, "Exactly matched reads: %d (%.2f %%)\n", n_match, 100.0 * n_match / n_seqs);

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

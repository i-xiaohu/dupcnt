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
#include "cstl/kthread.h"

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

/** Parallelism over I/O and processing */
typedef struct {
	bwaidx_t *idx;
	kseq_t *ks;
	int actual_chunk_size;
	int n_threads;
	int64_t n_sample_seqs;
	int64_t n_matched_seqs;
	Trie *trie_counter;
	int32_t *em_counter;

	// Time profiler
	double t_input, t_match, t_trie;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static void *dual_pipeline(void *shared, int step, void *_data) {
	double t_start, t_end;
	auto *aux = (ktp_aux_t*)shared;
	auto *data = (ktp_data_t*)_data;
	if (step == 0) { // Input
		t_start = realtime();
		auto *ret = (ktp_data_t*) calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, nullptr);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		for (int i = 0; i < ret->n_seqs; i++) {
			bseq1_t *b = &ret->seqs[i];
			// Check read length
			if (read_length_monitor == -1) read_length_monitor = b->l_seq;
			if (read_length_monitor != b->l_seq) {
				fprintf(stderr, "Only fix-length reads are supported\n");
				std::abort();
			}
			// Convert non-ACGT to A
			for (int k = 0; k < b->l_seq; k++) {
				b->seq[k] = nst_nt4_table[(uint8_t)b->seq[k]];
				if (b->seq[k] > 3) b->seq[k] = 0;
			}
		}
		aux->t_input = realtime() - t_start;
		return ret;
	}
	else if (step == 1) { // Process
		const bwaidx_t *idx = aux->idx;
		aux->n_sample_seqs += data->n_seqs;
		// 1. Identify exactly matched reads
		t_start = realtime();
		int n_rest = exact_matching(idx, data->n_seqs, data->seqs, aux->em_counter);
		aux->n_matched_seqs += data->n_seqs - n_rest;
		t_end = realtime();
		aux->t_match += t_end - t_start;
		fprintf(stderr, "%.2f %% matched reads\n", 100.0 * (double)(data->n_seqs - n_rest) / data->n_seqs);

		// 2. Process unmatched reads using trie
		// todo: trie is slower than sorting
		t_start = realtime();
		for (int j = 0; j < n_rest; j++) {
			bseq1_t *b = &data->seqs[j];
			aux->trie_counter->add_read(b->l_seq, b->seq);
			bseq1_destroy(b);
		}
		t_end = realtime();
		aux->t_trie += t_end - t_start;
		free(data->seqs);
		return 0;
	}
	return 0;
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

	ktp_aux_t aux;
	const int BATCH_SIZE = 10 * 1000 * 1000; // 10M bases for each thread
	aux.idx = idx;
	aux.actual_chunk_size = n_threads * BATCH_SIZE;
	aux.n_threads = n_threads;
	aux.em_counter = (int32_t*) calloc(ref_len, sizeof(int32_t));
	aux.trie_counter = new Trie();

	for (int i = 0; i < n_sample; i++) {
		gzFile fp = gzopen(files[i], "r");
		if (fp == nullptr) {
			fprintf(stderr, "Open FASTA file `%s` failed\n", files[i]);
			continue;
		}
		aux.t_input = aux.t_match = aux.t_trie = 0;
		aux.ks = kseq_init(fp);
		aux.n_sample_seqs = 0;
		aux.n_matched_seqs = 0;

		kt_pipeline(2, dual_pipeline, &aux, 2);

		t_start = realtime();
		aux.trie_counter->stat();
		t_end = realtime();

		fprintf(stderr, "Input: %.2f s, Match: %.2f s, Trie: %.2f s, Post: %.2f s\n", aux.t_input, aux.t_match, aux.t_trie, t_end - t_start);
		fprintf(stderr, "Exactly matched reads: %ld (%.2f %%)\n", aux.n_matched_seqs, 100.0 * aux.n_matched_seqs / aux.n_sample_seqs);

		int max_occ = 0;
		int64_t max_pos = -1;
		for (int64_t j = 0; j < ref_len; j++) {
			if (aux.em_counter[j] > max_occ) {
				max_occ = aux.em_counter[j];
				max_pos = j;
			}
		}
		fprintf(stderr, "Read matched at %ld occurs %d times\n", max_pos, max_occ);
	}

	free(aux.em_counter);
	delete aux.trie_counter;
}

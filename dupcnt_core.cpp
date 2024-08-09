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
#include <sys/resource.h>
#include "dupcnt_core.h"
#include "bwalib/bwa.h"
#include "cstl/kthread.h"

int read_length_monitor = -1; // To avoid any read of different length

#include "bwalib/kseq.h"
KSEQ_INIT(gzFile, gzread)

Trie::Trie() {
	unique_n = 0;
	nodes.clear();
	nodes.emplace_back(TrNode());
}

void Trie::add_read(int n, const char *s) {
	int32_t p = 0;
	for (int i = TRIE_SHIFT; i < n; i++) {
		uint8_t c = s[i];
		if (nodes[p].x[c] == 0) { // No child
			nodes[p].x[c] = nodes.size();
			nodes.emplace_back(TrNode());
		}
		p = nodes[p].x[c]; // Move down to the next layer
	}
	// Remember x[0] is initialized as 0 thus can be used to track the occurrence number by +1
	if (nodes[p].x[0] == 0) unique_n++;
	nodes[p].x[0]++;
}

int Trie::get_max_occ() {
	int32_t p = 0, ret = 0;
	std::vector<int> prev, curr;
	prev.push_back(p);
	for (int i = TRIE_SHIFT; i < read_length_monitor; i++) {
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
	for (auto k : prev) {
		if (nodes[k].x[0] > ret) {
			ret = nodes[k].x[0];
		}
	}
	return ret;
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline void bseq1_destroy(const bseq1_t *b) {
	free(b->name);
	free(b->comment);
	free(b->seq);
	free(b->qual);
}

typedef struct {
	const bwt_t *bwt;
	const bseq1_t *seqs;
	int64_t *match_pos; // -1 for unmatched reads
	Trie **trie_counter;
	std::vector<bseq1_t> *unmatched_seqs;
} worker_t ;

void exact_matching(void *_data, long seq_id, int t_id) {
	auto *w = (worker_t*)_data;
	const bwt_t *bwt = w->bwt;
	const bseq1_t *b = &w->seqs[seq_id];

	auto *s = (uint8_t*)b->seq;
	bwtintv_t ik, ok[4];
	bwt_set_intv(bwt, s[0], ik);
	for (int j = 1; j < b->l_seq; j++) {
		bwt_extend(bwt, &ik, ok, 0);
		ik = ok[3 - s[j]];
		if (ik.x[2] == 0) {
			break;
		}
	}
	if (ik.x[2] > 0) {
		// Use the first occurrence in suffix array
		w->match_pos[seq_id] = bwt_sa(bwt, ik.x[0]);
		bseq1_destroy(b); // Matched reads are deallocated here.
	} else {
		w->match_pos[seq_id] = -1;
	}
}

void trie_insert(void *_data, long bucket_id, int t_id) {
	auto *w = (worker_t*)_data;
	Trie *trie_counter = w->trie_counter[bucket_id];
	std::vector<bseq1_t> &seqs = w->unmatched_seqs[bucket_id];
	trie_counter->unique_n = 0; // Reset to calculate unique reads added
	for (auto &seq : seqs) {
		const bseq1_t *b = &seq;
		trie_counter->add_read(b->l_seq, b->seq);
		bseq1_destroy(b); // Unmatched reads are deallocated here.
	}
	seqs.clear(); // Clear the bucket after each trie counting
}

/** Parallelism over I/O and processing */
typedef struct {
	bwaidx_t *idx;
	kseq_t *ks;
	int actual_chunk_size;
	int n_threads;
	int64_t n_sample_seqs;
	int64_t n_matched_seqs;
	int64_t n_unique; // Number of left reads in all samples after deduplication
	Trie **trie_counter;
	std::vector<bseq1_t> *unmatched_seqs;
	int32_t *em_counter;
	int32_t max_occ; // Occurrence number of most repetitive reads
	int sample_id;
	int batch_id;

	// Time profiler
	double t_start, t_input, t_match, t_trie;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static void *dual_pipeline(void *shared, int step, void *_data) {
	double t_start, t_end, c_start, c_end;
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
		aux->t_input += realtime() - t_start;
		return ret;
	}
	else if (step == 1) { // Process
		aux->n_sample_seqs += data->n_seqs;
		// 1. Identify exactly matched reads
		double cpu_ratio1, cpu_ratio2;
		t_start = realtime();
		c_start = cputime();
		worker_t w;
		w.bwt = aux->idx->bwt;
		w.seqs = data->seqs;
		w.trie_counter = aux->trie_counter;
		w.unmatched_seqs = aux->unmatched_seqs;
		w.match_pos = (int64_t*) malloc(data->n_seqs * sizeof(int64_t));
		kt_for(aux->n_threads, exact_matching, &w, data->n_seqs);
		int n_rest = 0;
		for (int i = 0; i < data->n_seqs; i++) {
			int64_t pos = w.match_pos[i];
			if (pos != -1) {
				if (aux->em_counter[pos] == 0) aux->n_unique++;
				aux->em_counter[pos]++;
				aux->max_occ = std::max(aux->max_occ, aux->em_counter[pos]);
			} else {
				uint32_t prefix = 0;
				for (int j = 0; j < TRIE_SHIFT; j++) {
					prefix <<= 2U;
					prefix += data->seqs[i].seq[j];
				}
				w.unmatched_seqs[prefix].push_back(data->seqs[i]);
				n_rest++;
			}
		}
		free(w.match_pos);
		aux->n_matched_seqs += data->n_seqs - n_rest;
		t_end = realtime();
		c_end = cputime();
		aux->t_match += t_end - t_start;
		cpu_ratio1 = (c_end - c_start) / (t_end - t_start);

		// 2. Process unmatched reads using trie
		t_start = realtime();
		c_start = cputime();
		kt_for(aux->n_threads, trie_insert, &w, TRIE_BUCKET_SIZE);
		for (int i = 0; i < TRIE_BUCKET_SIZE; i++) {
			aux->n_unique += aux->trie_counter[i]->unique_n;
//			aux->max_occ = std::max(aux->max_occ, aux->trie_counter[i]->get_max_occ());
		}
		t_end = realtime();
		c_end = cputime();
		aux->t_trie += t_end - t_start;
		cpu_ratio2 = (c_end - c_start) / (t_end - t_start);
		fprintf(stderr, "[Sample %d Batch %d] %ld reads processed; %.2f seconds elapsed (Input %.2f, Match %.2f, Trie %.2f); Match %.2f; Trie %.2f\n",
		        aux->sample_id, aux->batch_id++, aux->n_sample_seqs,
		        t_end - aux->t_start, aux->t_input, aux->t_match, aux->t_trie,
		        cpu_ratio1, cpu_ratio2);
		free(data->seqs);
		free(data);
		return 0;
	}
	return 0;
}

void process(int n_threads, const char *index_prefix, int n_sample, char *files[]) {
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
	aux.n_sample_seqs = 0;
	aux.n_matched_seqs = 0;
	aux.n_unique = 0;
	aux.max_occ = 0;
	aux.em_counter = (int32_t*) calloc(ref_len, sizeof(int32_t));
	aux.t_start = realtime();
	aux.trie_counter = new Trie*[TRIE_BUCKET_SIZE];
	for (int i = 0; i < TRIE_BUCKET_SIZE; i++) aux.trie_counter[i] = new Trie();
	aux.unmatched_seqs = new std::vector<bseq1_t>[TRIE_BUCKET_SIZE];

	for (int i = 0; i < n_sample; i++) {
		gzFile fp = gzopen(files[i], "r");
		if (fp == nullptr) {
			fprintf(stderr, "Open FASTA file `%s` failed\n", files[i]);
			continue;
		}
		aux.t_input = aux.t_match = aux.t_trie = 0;
		aux.ks = kseq_init(fp);
		aux.sample_id = i + 1;
		aux.batch_id = 1;

		kt_pipeline(2, dual_pipeline, &aux, 2);

		fprintf(stderr, "%d sample(s) processed\n", i + 1);
		fprintf(stderr, "  Time profile(s):       input %.2f; Match %.2f; Trie %.2f\n", aux.t_input, aux.t_match, aux.t_trie);
		fprintf(stderr, "  Number of reads:       %ld\n", aux.n_sample_seqs);
		fprintf(stderr, "  Exactly matched reads: %ld (%.2f %%)\n", aux.n_matched_seqs, 100.0 * aux.n_matched_seqs / aux.n_sample_seqs);
		fprintf(stderr, "  Unique reads:          %ld (%.2f %%)\n", aux.n_unique, 100.0 * aux.n_unique / aux.n_sample_seqs);
//		fprintf(stderr, "  Max occurrence:        %d\n", aux.max_occ);
		fprintf(stderr, "\n");

		kseq_destroy(aux.ks);
		gzclose(fp);
	}

	free(aux.em_counter);
	for (int i = 0; i < TRIE_BUCKET_SIZE; i++) delete aux.trie_counter[i];
	delete [] aux.trie_counter;
	delete [] aux.unmatched_seqs;
}

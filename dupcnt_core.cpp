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
#include <stack>
#include "dupcnt_core.h"
#include "cstl/kthread.h"
#include "FM_index2/FMI_search.h"
#include "input.h"
#include "bwalib/kseq.h"
KSEQ_INIT(gzFile, gzread)

int read_length_monitor = -1; // To avoid any read of different length

typedef struct {
	const Option *opt;
	kseq_t *ks; // Input stream
	const bseq1_t *seqs; // Input reads
	FMI_search *fmi; // FM-index2
	int64_t *match_pos; // -1 for unmatched reads
	int32_t *em_counter; // em_counter[i]: count of reads matched at the position i
	std::vector<bseq1_t> *unmatched_seqs; // Unmatched reads are distributed to different buckets
	int32_t sample_id; // Current sample
	int32_t batch_id; // Current batch
	int64_t n_sample_seqs; // Total number of input reads
	int64_t n_matched_seqs; // Number of matched reads in all samples
	int64_t n_unique; // Number of left reads in all samples after deduplication
	double t_input, t_match, t_avl; // Wall clock time of each stage
} ktp_aux_t;

static inline void bseq1_destroy(const bseq1_t *b) {
	free(b->name);
	free(b->comment);
	free(b->seq);
	free(b->qual);
}

void exact_matching(void *_data, long seq_id, int t_id) {
	auto *w = (ktp_aux_t*)_data;
	FMI_search *fmi = w->fmi;
	const bseq1_t *b = &w->seqs[seq_id];
	auto *s = (uint8_t*)b->seq;

	SMEM ik = fmi->setInterval(3 - s[0]);
	for (int j = 1; j < b->l_seq; j++) {
		ik = fmi->backwardExt(ik, 3 - s[j]);
		if (ik.s == 0) {
			break;
		}
	}
	if (ik.s > 0) {
		// Swap the forward and reverse match results
		std::swap(ik.k, ik.l);
		// Use the first occurrence in suffix array
		w->match_pos[seq_id] = fmi->get_sa_entry_compressed(ik.k);
		// Matched reads are deallocated here.
		bseq1_destroy(b);
	} else {
		w->match_pos[seq_id] = -1;
	}
}

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static void *dual_pipeline(void *shared, int step, void *_data) {
	double t_start, t_end, c_start, c_end;
	auto *aux = (ktp_aux_t*)shared;
	auto *opt = aux->opt;
	auto *data = (ktp_data_t*)_data;
	if (step == 0) { // Input
		t_start = realtime();
		auto *ret = (ktp_data_t*) calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(opt->n_threads * opt->batch_size, &ret->n_seqs, aux->ks, nullptr);
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
			char r1[read_length_monitor + 1], r2[read_length_monitor + 1];
			r1[read_length_monitor] = r2[read_length_monitor] = '\0';
			for (int k = 0; k < b->l_seq; k++) {
				b->seq[k] = nst_nt4_table[(uint8_t)b->seq[k]];
				if (b->seq[k] > 3) b->seq[k] = 0; // Convert non-ACGT to A
			}
			for (int k = 0; k < b->l_seq; k++) {
				r1[k] = "ACGT"[b->seq[k]];
				r2[b->l_seq - 1 - k] = "ACGT"[3 - b->seq[k]];
			}
			if (strcmp(r1, r2) > 0) {
				// todo: check if something funny happens
				// Choose the lexicographically smaller one
				for (int k = 0; k < b->l_seq; k++) {
					b->seq[k] = nst_nt4_table[(uint8_t)r2[k]];
				}
			}
		}
		aux->t_input += realtime() - t_start;
		return ret;
	}
	else if (step == 1) { // Process
		// 1. Identify exactly matched reads
		double cpu_ratio1, cpu_ratio2;
		t_start = realtime(); c_start = cputime();
		aux->n_sample_seqs += data->n_seqs;
		aux->seqs = data->seqs;
		aux->match_pos = (int64_t*) malloc(data->n_seqs * sizeof(int64_t));
		kt_for(opt->n_threads, exact_matching, aux, data->n_seqs);
		int n_rest = 0;
		for (int i = 0; i < data->n_seqs; i++) {
			int64_t pos = aux->match_pos[i];
			if (pos != -1) {
				if (pos >= aux->fmi->reference_seq_len) {
					fprintf(stderr, "%ld\n", pos);
				}
				if (aux->em_counter[pos] == 0) aux->n_unique++;
				aux->em_counter[pos]++;
			} else {
				uint32_t prefix = 0;
				for (int j = 0; j < AVL_SHIFT; j++) {
					prefix <<= 2U;
					prefix += data->seqs[i].seq[j];
				}
				aux->unmatched_seqs[prefix].push_back(data->seqs[i]);
				n_rest++;
			}
		}
		free(aux->match_pos);
		aux->n_matched_seqs += data->n_seqs - n_rest;
		t_end = realtime(); c_end = cputime();
		aux->t_match += t_end - t_start;
		cpu_ratio1 = (c_end - c_start) / (t_end - t_start);

		// 2. Process unmatched reads using trie
		free(data->seqs);
		free(data);
		return NULL;
	}
	return NULL;
}

/** Post-processing */
void process(const Option *opt, int n_sample, char *files[]) {
	/* FM-index from BWA-MEM2  */
	auto *fmi = new FMI_search(opt->index_prefix);
	fmi->load_index();
	int64_t ref_len = fmi->idx->bns->l_pac * 2;
	fprintf(stderr, "Reference length: %ld\n", ref_len);

	ktp_aux_t aux;
	aux.opt = opt;
	aux.fmi = fmi;
	aux.em_counter = (int32_t*) calloc(ref_len + 1, sizeof(int32_t));
	aux.unmatched_seqs = new std::vector<bseq1_t>[AVL_BUCKET];
	aux.n_sample_seqs = 0;
	aux.n_matched_seqs = 0;
	aux.n_unique = 0;

	for (int i = 0; i < n_sample; i++) {
		gzFile fp = gzopen(files[i], "r");
		if (fp == nullptr) {
			fprintf(stderr, "Open FASTA file `%s` failed\n", files[i]);
			continue;
		}
		aux.t_input = aux.t_match = aux.t_avl = 0;
		aux.ks = kseq_init(fp);
		aux.sample_id = i + 1;
		aux.batch_id = 1;

		kt_pipeline(2, dual_pipeline, &aux, 2);

		fprintf(stderr, "%s added, %d sample(s) processed\n", files[i], i + 1);
		fprintf(stderr, "  Time profile1(s):      Input %.2f; Match %.2f; Trie %.2f\n", aux.t_input, aux.t_match, aux.t_avl);
		fprintf(stderr, "  Number of reads:       %ld\n", aux.n_sample_seqs);
		fprintf(stderr, "  Exactly matched reads: %ld (%.2f %%)\n", aux.n_matched_seqs, 100.0 * aux.n_matched_seqs / aux.n_sample_seqs);
		fprintf(stderr, "  Unique reads:          %ld (%.2f %%)\n", aux.n_unique, 100.0 * aux.n_unique / aux.n_sample_seqs);
		fprintf(stderr, "\n");

		kseq_destroy(aux.ks);
		gzclose(fp);
	}

	delete fmi;
	free(aux.em_counter);
	delete [] aux.unmatched_seqs;
}

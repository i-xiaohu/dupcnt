//
// Created by ixiaohu on 2024/8/7.
//

#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <string>
#include <cassert>
#include <stack>
#include <queue>
#include "dupcnt_core.h"
#include "cstl/kthread.h"
#include "FM_index2/FMI_search.h"
#include "input.h"
#include "kavl.h"
#include "kalloc.h"
#include "bwalib/kseq.h"
KSEQ_INIT(gzFile, gzread)

int read_length_monitor = -1; // To avoid any read of different length
#define PACK_SIZE 4

typedef struct read_node_s {
	uint8_t *packed; // Sequence packed with 2-bit encoding
	int32_t count; // Occurrence number of read
	KAVL_HEAD(struct read_node_s) head;
} read_node_t;

inline int read_node_cmp(const read_node_t *a, const read_node_t *b) {
	int pack_len = (read_length_monitor - AVL_SHIFT + PACK_SIZE - 1) / PACK_SIZE;
	for (int i = 0; i < pack_len; i++) {
		if (a->packed[i] < b->packed[i]) return -1;
		else if (a->packed[i] > b->packed[i]) return 1;
	}
	return 0;
}
#define macro_node_cmp(a, b) (read_node_cmp(a, b))
KAVL_INIT(rn, read_node_t, head, macro_node_cmp)

typedef struct {
	const Option *opt;
	kseq_t *ks; // Input stream
	const bseq1_t *seqs; // Input reads
	FMI_search *fmi; // FM-index2
	int64_t *match_pos; // -1 for unmatched reads
	int32_t *em_counter; // em_counter[i]: count of reads matched at the position i
	std::vector<bseq1_t> *unmatched_seqs; // Unmatched reads are distributed to different buckets
	void **km; // Memory block for each thread
	read_node_t **roots; // AVL trees
	int32_t *tree_counter; // Unique reads found in binary trees
	int32_t sample_id; // Current sample
	int32_t batch_id; // Current batch
	int64_t n_sample_seqs; // Total number of input reads
	int64_t n_matched_seqs; // Number of matched reads in all samples
	int64_t n_unique; // Number of left reads in all samples after deduplication
	double t_input, t_match, t_avl; // Wall clock time of each stage
} ktp_aux_t;

/** Step1: Find all exactly matched reads */
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

/** Step2: Use AVL to process unmatched reads */
void avl_search(void *_data, long seq_id, int t_id) {
	auto *w = (ktp_aux_t*)_data;
	void *km = w->km[t_id];
	read_node_t *root = w->roots[seq_id];
	auto &reads = w->unmatched_seqs[seq_id];

	for (auto &b : reads) {
		read_node_t *r;
		KCALLOC(km, r, 1);
		int32_t pack_len = (b.l_seq - AVL_SHIFT + PACK_SIZE - 1) / PACK_SIZE;
		KMALLOC(km, r->packed, pack_len);
		for (int i = 0, k = AVL_SHIFT; i < pack_len; i++) {
			uint8_t x = 0;
			for (int j = 0; j < PACK_SIZE; j++) {
				x <<= 2U;
				x |= (k < b.l_seq ?b.seq[k] :0U);
				k++;
			}
			r->packed[i] = x;
		}
		read_node_t *t = kavl_size(head, root) > 0 ?kavl_find(rn, root, r, NULL) :NULL;
		if (t) {
			t->count++;
			kfree(km, r->packed);
			kfree(km, r);
		} else {
			r->count = 1;
			kavl_insert(rn, &root, r, NULL);
			w->tree_counter[t_id]++;
		}
		bseq1_destroy(&b); // Unmatched reads are deallocated here
	}
	w->roots[seq_id] = root; // Remember root may change
	reads.clear(); // Clear the vector at each call
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
				// Read appears at the position for the first time
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

		// 2. Process unmatched reads using AVL tree
		t_start = realtime(); c_start = cputime();
		for (int i = 0; i < opt->n_threads; i++) aux->tree_counter[i] = 0;
		kt_for(opt->n_threads, avl_search, aux, AVL_BUCKET);
		for (int i = 0; i < opt->n_threads; i++) aux->n_unique += aux->tree_counter[i];
		t_end = realtime(); c_end = cputime();
		aux->t_avl += t_end - t_start;
		cpu_ratio2 = (c_end - c_start) / (t_end - t_start);

		fprintf(stderr, "[Sample %d Batch %d] %ld reads processed; (Input %.2f, Match %.2f, AVL %.2f); Match %.2f; AVL %.2f\n",
		        aux->sample_id, aux->batch_id++, aux->n_sample_seqs,
		        aux->t_input, aux->t_match, aux->t_avl,
		        cpu_ratio1, cpu_ratio2);

		free(data->seqs);
		free(data);
		return NULL;
	}
	return NULL;
}


void output(const ktp_aux_t *aux) {
	const Option *opt = aux->opt;
	const FMI_search *fmi = aux->fmi;
	int64_t ref_len = fmi->reference_seq_len - 1; // Excluding sentinel
	auto *ref_seq = (uint8_t*) malloc(ref_len * sizeof(uint8_t));
	int32_t *em_counter = aux->em_counter;

	{
		// Load [prefix].0123
		FILE *fp = fopen((std::string(opt->index_prefix) + ".0123").c_str(), "r");
		assert(fp != nullptr);
		fseek(fp, 0, SEEK_END);
		int64_t filesize = ftell(fp);
		assert(filesize == ref_len);
		fseek(fp, 0, SEEK_SET);
		fread(ref_seq, sizeof(uint8_t), ref_len, fp);
		fclose(fp);
	}

	if (not opt->sorted) {
		// Output exacted matched reads
		char buf[read_length_monitor + 1];
		buf[read_length_monitor] = '\0';
		for (int i = 0; i < ref_len - read_length_monitor; i++) {
			if (em_counter[i] > 1) {
				memcpy(buf, ref_seq + i, read_length_monitor * sizeof(uint8_t));
				for (int j = 0; j < read_length_monitor; j++) {
					buf[j] = "ACGT"[buf[j]];
				}
				fprintf(stdout, "%s %d\n", buf, em_counter[i]);
			}
		}
		// Output reads in AVL trees
		for (int i = 0; i < AVL_BUCKET; i++) {
			uint32_t x = i;
			for (int j = AVL_SHIFT - 1; j >= 0; j--) {
				buf[j] = "ACGT"[x & 3U];
				x >>= 2U;
			}
			read_node_t *root = aux->roots[i];
			std::queue<read_node_t*> que;
			que.push(root);
			while (not que.empty()) {
				read_node_t *f = que.front();
				que.pop();
				if (f->count > 1) {
					for (int j = 0; j < read_length_monitor - AVL_SHIFT; j++) {
						uint8_t p = f->packed[j / PACK_SIZE];
						uint8_t q = (p >> (2 * (3 - (j % PACK_SIZE)))) & 3U;
						buf[j + AVL_SHIFT] = "ACGT"[q];
					}
					fprintf(stdout, "%s %d\n", buf, f->count);
				}
				if (f->head.p[0]) que.push(f->head.p[0]);
				if (f->head.p[1]) que.push(f->head.p[1]);
			}
		}
	}

	free(ref_seq);
}

void process(const Option *opt, int n_sample, char *files[]) {
	/* FM-index from BWA-MEM2 */
	auto *fmi = new FMI_search(opt->index_prefix);
	fmi->load_index();

	ktp_aux_t aux;
	aux.opt = opt;
	aux.fmi = fmi;
	aux.em_counter = (int32_t*) calloc(fmi->reference_seq_len, sizeof(int32_t));
	aux.unmatched_seqs = new std::vector<bseq1_t>[AVL_BUCKET];
	aux.km = (void**) malloc(opt->n_threads * sizeof(void*));
	for (int i = 0; i < opt->n_threads; i++) aux.km[i] = km_init();
	aux.roots = (read_node_t**) calloc(AVL_BUCKET, sizeof(read_node_t*));
	aux.tree_counter = (int32_t*) calloc(opt->n_threads, sizeof(int32_t));
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
		fprintf(stderr, "  Time profile1(s):      Input %.2f; Match %.2f; AVL %.2f\n", aux.t_input, aux.t_match, aux.t_avl);
		fprintf(stderr, "  Number of reads:       %ld\n", aux.n_sample_seqs);
		fprintf(stderr, "  Exactly matched reads: %ld (%.2f %%)\n", aux.n_matched_seqs, 100.0 * aux.n_matched_seqs / aux.n_sample_seqs);
		fprintf(stderr, "  Unique reads:          %ld (%.2f %%)\n", aux.n_unique, 100.0 * aux.n_unique / aux.n_sample_seqs);
		fprintf(stderr, "\n");

		kseq_destroy(aux.ks);
		gzclose(fp);
	}

	// Output repetitive reads with frequency
	output(&aux);

	delete fmi;
	free(aux.em_counter);
	delete [] aux.unmatched_seqs;
	for (int i = 0; i < opt->n_threads; i++) km_destroy(aux.km[i]);
	free(aux.km);
	free(aux.roots);
	free(aux.tree_counter);
}

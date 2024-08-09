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
#include "bwalib/bwa.h"
#include "cstl/kthread.h"

int read_length_monitor = -1; // To avoid any read of different length

#include "bwalib/kseq.h"
KSEQ_INIT(gzFile, gzread)

Trie::Trie() {
	unique_n = 0;
	overflow = true;
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

void Trie::auto_adjust_size() {
	if (overflow and nodes.size() > TRIE_SIZE_CAP) {
//		fprintf(stderr, "Trie size: %ld oversize (%ld bytes)\n", nodes.size(), nodes.capacity() * sizeof(TrNode));
		std::vector<int> parent; // For backtrace
		std::vector<bool> kept; // True if node is kept
		std::vector<int> shift; // For coordinate shift
		std::vector<TrNode> new_nodes;
		std::vector<int> prev, curr;
		int kept_n;
		parent.resize(nodes.size());
		kept.resize(nodes.size(), false);
		prev.push_back(0); // Root
		for (int i = TRIE_SHIFT; i < read_length_monitor; i++) {
			curr.clear();
			for (int p : prev) {
				for (int c : nodes[p].x) {
					if (c) {
						parent[c] = p;
						curr.push_back(c);
					}
				}
			}
			std::swap(prev, curr);
		}
		// Sort by frequency
		std::sort(prev.begin(), prev.end(),
			[&](int a, int b) -> bool { return nodes[a].x[0] > nodes[b].x[0]; }
		);
		kept[0] = true;
		kept_n = 1;
		for (int x: prev) {
			while (x != 0) {
				if (kept[x]) break; // All ancestors kept
				kept[x] = true;
				kept_n++;
				x = parent[x];
			}
			// Half the size
			if (kept_n > nodes.size() >> 1U) {
				break;
			}
		}
		parent.clear();
		prev.clear();
		curr.clear();

		// Shifting coordinate
		shift.resize(nodes.size());
		shift[0] = 0;
		kept_n = 1;
		for (int i = 1; i < nodes.size(); i++) {
			if (kept[i]) {
				shift[i] = kept_n++;
				new_nodes.push_back(nodes[i]);
			} else shift[i] = 0;
		}
		for (auto &p : new_nodes) {
			// Update child index
			// It is correct because c = 0, shift[c] = 0 and kept[c] = 0, shift[c] = 0;
			for (int &c : p.x) {
				c = shift[c];
			}
		}
		nodes = new_nodes;

		kept.clear();
		shift.clear();
	}
}

struct Tuple {
	int node_id;
	int depth;
	int label;
	Tuple(int n, int d, int l) : node_id(n), depth(d), label(l) {}
};

std::vector<RepRead> Trie::most_k_frequent(uint32_t bucket_id, int k) {
	std::vector<RepRead> heap;
	char seq[read_length_monitor + 1];
	seq[read_length_monitor] = '\0';
	std::stack<Tuple> st;
	for (int l = 3; l >= 0; l--) { // Preorder
		int c = nodes[0].x[l];
		if (c) st.push(Tuple(c, TRIE_SHIFT, l));
	}
	for (int i = TRIE_SHIFT-1; i >= 0; i--) {
		seq[i] = "ACGT"[bucket_id & 3U];
		bucket_id >> 2U;
	}

	while (not st.empty()) {
		auto t = st.top();
		st.pop();
		seq[t.depth] = "ACGT"[t.label];
		if (t.depth == read_length_monitor - 1) { // Leaf layer
			RepRead r;
			r.occ = nodes[t.node_id].x[0];
			r.read = std::string(seq);
			if (heap.size() < k) {
				heap.push_back(r);
				std::push_heap(heap.begin(), heap.end());
			} else if (r.occ > heap[0].occ) {
				heap[0] = r;
				std::make_heap(heap.begin(), heap.end());
			}
		} else {
			for (int l = 3; l >= 0; l--) {
				int c = nodes[t.node_id].x[l];
				if (c) st.push(Tuple(c, t.depth + 1, l));
			}
		}
	}

	return heap;
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

	trie_counter->auto_adjust_size();
}

/** Parallelism over I/O and processing */
typedef struct {
	bwaidx_t *idx;
	kseq_t *ks;
	const Option *opt;
	int64_t n_sample_seqs;
	int64_t n_matched_seqs;
	int64_t n_unique; // Number of left reads in all samples after deduplication
	Trie **trie_counter;
	std::vector<bseq1_t> *unmatched_seqs;
	int32_t *em_counter;
	int oversize_n;
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
		kt_for(opt->n_threads, exact_matching, &w, data->n_seqs);
		int n_rest = 0;
		for (int i = 0; i < data->n_seqs; i++) {
			int64_t pos = w.match_pos[i];
			if (pos != -1) {
				if (aux->em_counter[pos] == 0) aux->n_unique++;
				aux->em_counter[pos]++;
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
		size_t vram = 0;
		for (int i = 0; i < TRIE_BUCKET_SIZE; i++) vram += aux->trie_counter[i]->get_size();
		aux->oversize_n += vram > opt->mem_cap;
		for (int i = 0; i < TRIE_BUCKET_SIZE; i++) aux->trie_counter[i]->overflow = vram > opt->mem_cap;
		kt_for(opt->n_threads, trie_insert, &w, TRIE_BUCKET_SIZE);
		for (int i = 0; i < TRIE_BUCKET_SIZE; i++) aux->n_unique += aux->trie_counter[i]->unique_n;
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

void pickup_frequent(int k, Trie **trie_counter, const int *em_counter, const bwaidx_t *idx) {
	if (k <= 0) return ;
	std::vector<RepRead> heap;
	for (int i = 0; i < TRIE_BUCKET_SIZE; i++) {
		auto sub_heap = trie_counter[i]->most_k_frequent(i, k);
		for (auto &r : sub_heap) {
			if (heap.size() < k) {
				heap.push_back(r);
				std::push_heap(heap.begin(), heap.end());
			} else if (r.occ > heap[0].occ) {
				heap[0] = r;
				std::make_heap(heap.begin(), heap.end());
			}
		}
	}
	for (int64_t i = 0; i < idx->bns->l_pac << 1; i++) {
		RepRead r;
		r.occ = em_counter[i];
		if (heap.size() < k or r.occ > heap[0].occ) {
			int64_t len;
			uint8_t *seq = bns_get_seq(idx->bns->l_pac, idx->pac, i, i + read_length_monitor, &len);
			assert(len == read_length_monitor);
			r.read.resize(len);
			for (int j = 0; j < len; j++) r.read[j] = "ACGT"[seq[j]];
			free(seq);
		}
		if (heap.size() < k) {
			heap.push_back(r);
			std::push_heap(heap.begin(), heap.end());
		} else if (r.occ > heap[0].occ) {
			heap[0] = r;
			std::make_heap(heap.begin(), heap.end());
		}
	}
	std::sort(heap.begin(), heap.end());
	fprintf(stderr, "The %d most frequent read:\n", k);
	for (auto &r : heap) {
		fprintf(stderr, "  %s %d\n", r.read.c_str(), r.occ);
	}
	fprintf(stderr, "\n");
}

void process(const Option *opt, int n_sample, char *files[]) {
	/* FM-index from BWA-MEM; todo: do consider switch to BWA-MEM2 index */
	bwaidx_t *idx = bwa_idx_load_from_shm(opt->index_prefix);
	if (idx == nullptr) {
		if ((idx = bwa_idx_load(opt->index_prefix, BWA_IDX_ALL)) == nullptr) {
			fprintf(stderr, "Load index `%s` failed\n", opt->index_prefix);
			return ;
		}
	}
	int64_t ref_len = idx->bns->l_pac * 2;
	fprintf(stderr, "Reference length: %ld\n", ref_len);

	ktp_aux_t aux;
	aux.idx = idx;
	aux.opt = opt;
	aux.n_sample_seqs = 0;
	aux.n_matched_seqs = 0;
	aux.n_unique = 0;
	aux.em_counter = (int32_t*) calloc(ref_len, sizeof(int32_t));
	aux.t_start = realtime();
	aux.trie_counter = new Trie*[TRIE_BUCKET_SIZE];
	for (int i = 0; i < TRIE_BUCKET_SIZE; i++) aux.trie_counter[i] = new Trie();
	aux.unmatched_seqs = new std::vector<bseq1_t>[TRIE_BUCKET_SIZE];
	aux.oversize_n = 0;

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
		fprintf(stderr, "  Oversize times:        %d\n", aux.oversize_n);
		fprintf(stderr, "\n");

		// todo: parallel the function
		pickup_frequent(opt->most_rep, aux.trie_counter, aux.em_counter, aux.idx);

		kseq_destroy(aux.ks);
		gzclose(fp);
	}

	free(aux.em_counter);
	for (int i = 0; i < TRIE_BUCKET_SIZE; i++) delete aux.trie_counter[i];
	delete [] aux.trie_counter;
	delete [] aux.unmatched_seqs;
}

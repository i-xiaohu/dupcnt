//
// Created by ixiaohu on 2024/8/7.
//

#ifndef DUPCNT_DEDUP_CORE_H
#define DUPCNT_DEDUP_CORE_H

#include <cstdint>
#include <utility>
#include <vector>
#include <string>

/**
 * The first 6 bases of reads are used to distributed to 4^6 = 4096 tries.
 * Each trie has a maximum of 1,953,125 nodes.
 * If trie is oversized, some nodes and paths will be removed to restrict memory usage.
 **/

#define TRIE_SHIFT        6
#define TRIE_BUCKET_SIZE  4096
#define TRIE_SIZE_CAP     1953125

struct TrNode {
	// x[c] = 0 suggest no child down from branch c
	// It is correct since root node always occupies the index of 0
	int32_t x[4] = {0};
};

struct RepRead {
	int occ;
	std::string read;
	bool operator < (const RepRead &r) const {
		return occ > r.occ; // Put less frequent read at the top of heap
	}
};

class Trie {
private:
	std::vector<TrNode> nodes; /** Nodes in the trie */

public:
	int unique_n; // Number of reads that have no identical match in the trie
	bool overflow;

	/** Construction: init root node */
	Trie();

	/** Add a read into the trie that the first 6 bases of read are omitted and the occurrence
	 * numbers of leaves (reads) are recorded.
	 * @param n  read length
	 * @param s  2-bit encoded sequence
	 */
	void add_read(int n, const char *s);

	int get_max_occ();

	void auto_adjust_size();

	size_t get_size() { return nodes.capacity() * sizeof(TrNode); }

	std::vector<RepRead> most_k_frequent(uint32_t bucket_id, int k);
};

struct Option {
	int n_threads = 16;
	int batch_size = 10 * 1000 * 1000; // 10M bases for each thread
	int most_rep = 0;
	const char *index_prefix = nullptr;
	size_t mem_cap = 100L * 1024L * 1024L * 1024L; // 100GB
};

void process(const Option *opt, int n_sample, char *files[]);

#endif //DUPCNT_DEDUP_CORE_H

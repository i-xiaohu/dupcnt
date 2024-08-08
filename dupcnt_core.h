//
// Created by ixiaohu on 2024/8/7.
//

#ifndef DUPCNT_DEDUP_CORE_H
#define DUPCNT_DEDUP_CORE_H

#include <cstdint>
#include <vector>

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

class Trie {
private:
	std::vector<TrNode> nodes; /** Nodes in the trie */

public:
	int unique_n = 0;

	/** Construction: init root node */
	Trie();

	/** Add a read into the trie that the first 6 bases of read are omitted and the occurrence
	 * numbers of leaves (reads) are recorded.
	 * @param n  read length
	 * @param s  2-bit encoded sequence
	 */
	void add_read(int n, const char *s);

	/** Return unique reads in the trie */
	int count_unique();
};

void process(int n_threads, const char *index_prefix, int n_sample, char *files[]);

#endif //DUPCNT_DEDUP_CORE_H

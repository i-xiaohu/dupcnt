//
// Created by ixiaohu on 2024/8/19.
//

#ifndef DUPCNT_INPUT_H
#define DUPCNT_INPUT_H

// todo: Input function is confusing, code cleanup needed.

#include "FM_index2/bseq1.h"
#include <cstdlib>

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

inline void bseq1_destroy(const bseq1_t *b) {
	free(b->name);
	free(b->comment);
	free(b->seq);
	free(b->qual);
}

#endif //DUPCNT_INPUT_H

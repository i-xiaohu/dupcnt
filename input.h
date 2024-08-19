//
// Created by ixiaohu on 2024/8/19.
//

#ifndef DUPCNT_INPUT_H
#define DUPCNT_INPUT_H

// todo: Input function is confusing, code cleanup needed.

#include "bwalib/bwa.h"

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

#endif //DUPCNT_INPUT_H

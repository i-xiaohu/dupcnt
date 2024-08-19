//
// Created by ixiaohu on 2024/8/19.
//

#ifndef DUPCNT_BSEQ1_H
#define DUPCNT_BSEQ1_H

/** Input read structure from bwa.h */

typedef struct {
	int l_seq, id;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

#endif //DUPCNT_BSEQ1_H

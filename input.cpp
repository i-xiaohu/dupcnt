//
// Created by ixiaohu on 2024/8/19.
//

#include <cstdio>
#include <zlib.h>
#include "input.h"

/************************
 * Batch FASTA/Q reader *
 ************************/

#include "bwalib/kseq.h"
KSEQ_INIT(gzFile, gzread)

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

static inline char *dupkstring(const kstring_t *str, int dupempty)
{
	char *s = (str->l > 0 || dupempty)? (char*) malloc(str->l + 1) : NULL;
	if (!s) return NULL;

	memcpy(s, str->s, str->l);
	s[str->l] = '\0';
	return s;
}

static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
	s->name = dupkstring(&ks->name, 1);
	s->comment = dupkstring(&ks->comment, 0);
	s->seq = dupkstring(&ks->seq, 1);
	s->qual = dupkstring(&ks->qual, 0);
	s->l_seq = ks->seq.l;
}

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_)
{
	kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
			fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
			break;
		}
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t *) realloc(seqs, m * sizeof(bseq1_t));
		}
		trim_readno(&ks->name);
		kseq2bseq1(ks, &seqs[n]);
		seqs[n].id = n;
		size += seqs[n++].l_seq;
		if (ks2) {
			trim_readno(&ks2->name);
			kseq2bseq1(ks2, &seqs[n]);
			seqs[n].id = n;
			size += seqs[n++].l_seq;
		}
		if (size >= chunk_size && (n&1) == 0) break;
	}
	if (size == 0) { // test if the 2nd file is finished
		if (ks2 && kseq_read(ks2) >= 0)
			fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	}
	*n_ = n;
	return seqs;
}

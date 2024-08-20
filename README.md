# Introduction
**Count unique reads from multiple samples in an efficient way.**

*dupcnt* deduplicates reads in two steps.

1. Pick out repetitive reads that exactly match to HG19/HG38 reference utilizing FM-index2 (16GB). A 32-bit integer array is 
used to count the occurrence of reads at each reference position, requiring 24G for human genome.
2. A certain number of AVL trees are used to find duplications in the rest of unmatched reads. Specifically, reads in the tree 
share the same first 6 bases, facilitating parallelism over 4096 trees. Two-bit encoded reads as nodes are inserted into the binary
search trees, with each node occupying 40 plus (*read_length* / 4) bytes. For 6G reads length of 100bp (~200X), it needs < 400GB memory.
Considering not exactly matched reads accounts for <50% of total reads, a medium-end server of 256GB RAM could be fitted in.
It is a better solution to obtain absolutely accurate results than previously used trie that cannot guarantee correctness 
because of memory limitation.

# Usage
```
dupcnt [options] Sample_1 ... Sample_n
```

**Options**
- `-i` FM-index (BWA-MEM2 index, use `mem2idx` to create it if not available)
- `-t` Number of threads (16 by default)
- `-k` Input batch size (10M)

# Result
*dupcnt* is as fast as data input in overall performance.
After each sample processed, *dupcnt* outputs statistic information. Of note, `Number of reads`, `Exactly matched reads`
and `Unique reads` are accumulative values.
```
1 sample(s) processed
  Time profile(s):       input 3061.82; Match 1465.94; Trie 268.45
  Number of reads:       863941088
  Exactly matched reads: 693426503 (80.26 %)
  Unique reads:          779264875 (90.20 %)
```

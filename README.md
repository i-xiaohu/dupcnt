# Introduction
**How to count unique reads from multiple samples in an efficient way.**

*dupcnt* deduplicates reads in two steps.

1. Pick out repetitive reads that exactly match to HG19/HG38 reference utilizing suffix array. Such reads can be easily encoded with
64-bit integers (alignment position). Or use small integers to count the occurrence of reads at each reference position.
2. A certain number of tries are used to find duplications in the rest of reads. Here, reads in a trie share the same 6 bases, 
facilitating parallelism over 4096 tries.

To avoid excessive memory usage, the sizes of tries are capped. The less frequently visited nodes or paths will be deleted 
once the trie is oversized. Though it would make the duplication counting not accurate, it shouldn't have a significant 
impact on the overall result.

# Usage
**Options**
- `-i` FM-index (BWA index currently)
- `-t` Number of threads (16 by default)
- `-c` Trie maximum size to curb memory consumption (100GB by default, but needing ~150GB in practice)
- `-k` Input batch size (10M)
- `-j` Output the most p frequent reads (0 by default, -1 to print all reads appearing more than once)

# Result
*dupcnt* outputs time profile after each batch. Trie is super fast and well parallelized (64 threads used here).
Exact matching has a poor CPU usage, causing a major bottleneck. But *dupcnt* is as fast as data input in overall performance.
```
[Sample 1 Batch 1] 6336634 reads processed; 44.44 seconds elapsed (Input 42.17, Match 27.44, Trie 1.04); Match 15.46; Trie 58.16
[Sample 1 Batch 2] 12673268 reads processed; 64.16 seconds elapsed (Input 42.17, Match 46.04, Trie 2.15); Match 20.70; Trie 50.65
[Sample 1 Batch 3] 19009902 reads processed; 84.06 seconds elapsed (Input 63.88, Match 62.77, Trie 3.32); Match 23.22; Trie 52.81
...
[Sample 1 Batch 135] 855445590 reads processed; 3037.43 seconds elapsed (Input 3024.12, Match 1454.86, Trie 259.46); Match 32.78; Trie 43.24
[Sample 1 Batch 136] 861782224 reads processed; 3060.70 seconds elapsed (Input 3048.44, Match 1462.91, Trie 261.38); Match 36.30; Trie 45.84
[Sample 1 Batch 137] 863941088 reads processed; 3074.21 seconds elapsed (Input 3061.82, Match 1465.94, Trie 268.45); Match 32.07; Trie 4.22
```

After each sample processed, *dupcnt* outputs statistic information, Of note, `Number of reads`, `Exactly matched reads`
and `Unique reads` are accumulative values.
```
1 sample(s) processed
  Time profile(s):       input 3061.82; Match 1465.94; Trie 268.45
  Number of reads:       863941088
  Exactly matched reads: 693426503 (80.26 %)
  Unique reads:          779264875 (90.20 %)
```

After all samples added, *dupcnt* outputs the k most frequent reads with frequency.

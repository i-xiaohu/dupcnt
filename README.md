# Introduction
**How to count duplicated reads from multiple samples in an efficient way.**

*dupcnt* deduplicates reads in two steps.

1. Pick out repetitive reads that exactly match to HG38 reference utilizing suffix array. Such reads can be easily encoded with
64-bit integers (alignment position). Or use small integers to count the occurrence of reads at each reference base.
2. A certain number of tries are used to find the duplication in the rest of reads. Reads in a trie share a common characteristic, 
such as the same hash value or K prefix bases.

To avoid excessive memory usage, the sizes of tries are capped. The less frequently visited nodes or paths will be deleted 
once the trie is oversized. Though it would make the duplication counting not accurate, it shouldn't have a significant 
impact on the overall result.
=============================================
Pairwise Optimal Sequence Alignment
=============================================

This program uses the dynamic programming based sequence alignment algorithm, with the following three variations: global alignment, local alignment, and global alignment with affine gap penalties.

This takes O(n^2) time to construct the matrix and linear time to extract the optimal alignment from the matrix constructed. The alignment problems are only for DNA sequences.


Sample Input & Output

Input: ACCTG and CACTG

Output: GC_CTG
        _CACTG
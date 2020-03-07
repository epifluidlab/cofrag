# cf-HiC

Project goal: reconstruct 3D genome from cell-free DNA fragmentation patterns.

cell-free DNA fragmentation may carry much information revealing chromatin interactions. By analyzing the statistical distance between fragment distributions from two separated bins, we may infer how the two loci are interacting with each other. This is similar to the idea of Hi-C.

Seven steps are involved:

1. Calculate the distance matrix
2. Calculate the corrected correlation matrix, by using the Hi-C contact probability curve to correct the distance matrix
3. Calculate the similarity
### [WGCNA analysis](https://github.com/talimass/Cayman-translocation/tree/main/RNAseq_analysis/WGCNA)
VST-transformed counts were analyzed with WGCNA v.1.73 R package to identify modules of co-expressed genes and their association with experimental treatments. Network construction and module detection were performed with the `blockwiseModules` function. Identified modules were further grouped into clusters based on eigengene similarity. Differences between treatments within each cluster were tested using limma v.3.64.


# Bowtie2 Seminar
This is a repository for a seminar on Bowtie2 and its use of the Burrows-Wheeler Transform and FM-index in Short-read sequence mapping.

## Bowtie2

4 Steps:

1. Seed Extraction: 
    
    Substrings of read extracted with or without overlap. Seed length (parameter -L) can be set anywhere from 4-32 (20-25 with good rapport). When read length varies, it is better to use a seed length as a sublinear function of read length: 
    
    ```
    I(x) = max(1, floor(1 + 1.15 * √x)) 
    ```
    
    —> Bowtie2 end-to-end mode default (-i option to configure)
    
2. Seed alignment using FM-index:
    
    Finds ungapped alignments using the same reference pruning, policy pruning, and double indexing as Bowtie 1.
    
    Bowtie 2 also uses bidirectional BWT to switch between left-to-right and right-to-left alignment. Alignment can be done with up to 1 mismatch (-N being 0 or 1). For each seed, a set of 0 or more Burrows-Wheeler ranges, or seed-hit ranges, are returned (with multiple ranges for 1 seed if mismatches are allowed or if seed spans a repetitive sequence).
    
    When this step is done, Bowtie2 calculates the average number of seed hits per seed string. By default, if this average rises above 1000, reseeding occurs at successive rounds (for example seed from 3’ to 5’ starting at 0, next round at 5, next at 10...) whose maximum can be configured (-R option). This increases the chances that the best alignment is found with repetitive sequences.
    
3. Seed alignment prioritization:
    
    Each seed-hit range spans an amount or rows. The rows are scored according to the range of the hit: 
    
    ```
    1/r2
    ```
    
    where r is the total number of rows in the range. Then, rows are chosen in a random weighted fashion to be LF mapped to the reference by calculating the offset using the walk-left procedure. Each resolved offset as well as information on the seed that generated it are passed to step 4.
    
4. SIMD-accelerated dynamic programming:
    
    For the resolved seed-hit offsets, flanking sequences from the reference are considered for possibly gapped alignment using Needleman-Wunsch or Smith-Waterman dynamic programming algorithms. Although Smith-Waterman is optimal, it is quite slow (Farrar, 2007). This is why parallel computing is used. A basic example would be since only the values above, to the left, and diagonal upwards and to the left of the position at hand are needed, a matrix can be filled in an anti-diagonal manner using parallel programming. In the case of Bowtie2, a striped, vertical pattern is employed, as the one in striped swsse2 tool for protein alignments. Reads are processed into a query profile and diagonal scores are derived efficiently.
    
    Repetitive sequences can cause problems with dynamic programming. Bowtie2 allows adjustment of the number of trials that can “fail” (-D option) before reporting and moving on to the next read. Failure is when the attempt’s alignment score is smaller than the best or second best alignment score found so far.
    

Gap penalties for the read and reference can be separately adjusted. Penalties for gap continuation and initiation can also be configured separately.

## Bowtie version history

- Bowtie 2 is optimized to handle reads without a maximum read length, unlike Bowtie 1.
- Bowtie 1 does not permit gaps, whereas Bowtie 2 does.
- Bowtie 1 can only perform global alignment whereas Bowtie 2 can perform both global and local (by trimming either end of the read)
- Bowtie1 (Langmead et al., 2009):
    - assumptions mostly applicable to mammalian sequencing efforts with short reads
    - FM-index
    - employs quality-aware backtracking algorithm to allow few mismatches only within high-quality regions
    - quality-aware, greedy (outputs first alignment with good enough score unless —best option specified), randomized, depth-first search
    - Excessive backtracking can occur when reads are of low-quality

## Other Algorithms

- SOAP2 (Li et al., 2009):
    - Short oligonucleotide alignment program 2
    - BWT
    - reference indices hashed
    - allows for gaps and mismatches which can be configured
    - can support longer reads (up to 1024bp)
    - supports single and paired-end reads
- BWA (Li et al., 2010):
    - Burrows-Wheeler aligner
    - FM-index for both query and reference
    - prefix trie for reference and prefix Directed acyclic word graph (DAWG) for query
    - dynamic programming for alignment with heuristics: traverse the prefix DAWG first, and prune low scoring matches at each node.
    - SW extension on filtered seeds, except for highly repetitive sequences
    - support for repetitive sequences by keeping track of SA interval
    - assumes matches are high scoring (more correctly aligned)
- STAR (Dobin et al., 2012):
    - Spliced transcripts alignment to a reference
    - designed for RNAseq and tackling splice junctions as well as reads that span the length of a transcript (longer lengths)
    - Suffix array
    - employs Maximal Mappable Prefix (MMP) similar to MUM, this finds maximal length seed and repeats the process for unmapped portions of the read to map transcripts directly to the reference genome while accounting for splice sites. (This also accounts for mismatches and deletions.
    - support for high sequencing error rates
    - In the second step, seed matches are “stitched” together after clustering around a user-defined ‘anchor’ seed which defines maximal intron length. This allows support for paired-end reads by representing them as separate reads whose seed matches will fall within the same genomic window. Stitching is guided by user-defined scoring theme
    - Support for chimeric alignments
    - can make use for annotation of splice junctions as aid
- TopHat2 (Kim et al., 2013):
    - For transcriptomes and RNAseq data
    - can also make use annotation of splice junctions as aid
    - can handle datasets of reads with variable lengths (supports merging of datasets)
- MUMmer4 (Marçais et al., 2018):
    - Maximal Unique Match (MUM) and Maximal Exact Match (MEM; not suitable with large amounts of repeats) similar to MMP
    - Suffix array
    - clustering of matches
    - SW alignment
    - parallelized
- HISAT2 (Kim et al., 2019):
    - Hierarchical indexing scheme alignment tool 2
    - Long repeats projected to a single location, then later retrieved from FM index
    - support for HLA typing, DNA-fingerprinting analysis
    - builds linear graph of reference genome, and adds alternative paths for mutations, deletions, insertions. This graph is then turned into a prefix sorted graph
    - Hierarchical indexing at global and local levels helps spot repeats

## Comparison

Most algorithms can be configured through multiple parameters that control time and space complexity trade-off. Some methods also have multiple modes, each of which with a particular advantage conveyed through specific default parameters. Bowtie2, for example, has an end-to-end (global) mode as well as a local mode that is considerably slower but generates a higher % reads matched. Moreover, users can toggle between very sensitive, sensitive, fast, and very fast modes depending on application and preference (where the trade-off lies between performance and time complexity).

- Bowtie2 E2E: moderate runtime (~20-25 microsecs/read)
- Bowtie2 local: highest runtime (~47 microsecs/read)
- BWA: Similar % reads aligned to Bowtie2 local, moderate runtime (~20-25 microsecs/read)
- HISAT2: Similar % reads aligned to Bowtie2 E2E, lowest runtime/read (<10 microsecs/read)
- STAR: Similar % reads aligned to MUMmer4, less than Bowtie2 local and BWA but more than Bowtie2 E2E and HISAT2, moderate runtime (~20-25 microsecs/read)
- MUMmer4: Similar % reads aligned to STAR, less than Bowtie2 local and BWA but more than Bowtie2 E2E and HISAT2, second highest runtime (~38 microsecs/read)
- TopHat2: not used anymore, only used for comparison. lowest scores in terms of quality for time, space, and % reads aligned

Speed is usually not the primary concern, but if it is, such as in instantaneous alignment to monitor and treat ICU patients, HISAT2 represents the best option with very low runtime and performance similar to Bowtie2 E2E.

For performance, BWA seems to constitute a very good option, with runtime similar to Bowtie2 E2E and % aligned reads similar to Bowtie2 local (best of both Bowtie2 worlds).

For transcriptomic and RNAseq data where fragments can span splice junctions, tools like HISAT2 and STAR are designed for this purpose.

For further investigation of the alignments, such as with downstream analysis, SAM/BAM file outputs are needed. SAM (sequence alignment map) files are made for human readers. Packages like samtools are used to convert SAM to BAM (binary alignment map) files which are computer-readable. MUMmer4 is the only tool from the list that does not output a file in that category and requires some preprocessing before further analysis.

Aligners also differ in the way they handle multireads, which could occur with repetitive sequences or duplications.

For parallel processing, Bowtie2 local, BWA, HISAT2, and MUMmer4 achieved linear or almost linear parallel speedup. Bowtie2 E2E achieved slightly super-linear speedup. STAR and TopHat2 showed logarithmic speed up.

For transcript size, all except TopHat2 achieved >90% transcriptome coverage with alignments of at least 100 nt (BWA then Bowtie2 local and E2E). For longer cutoffs, such as >1000 nt, STAR and HISAT2 had better coverage (Musich et al., 2021).

- Bowtie
- Bowtie2: moderate peak memory
- SOAP2: highest peak memory
- BWA: lowest peak memory
- BWA-SW

In unpaired alignment, Bowtie2 showed greater % aligned reads and speed as well as intermediate memory usage.

In paired-end alignment, Bowtie is at a disadvantage because it does not perform gapped alignment. Bowtie2 also showed greater % aligned reads and speed as well as intermediate memory usage.

For longer reads (Roche 454 and Ion Torrent reads), Bowtie, BWA, and SOAP2 are also at a disadvantage as they are optimized for shorter read lengths. Bowtie2 showed greater % aligned reads and speed as well as slightly lower memory footprint than BWA-SW.

For shorter read length accuracy, BWA and Bowtie2 showed greater cumulative correct alignment than Bowtie and SOAP2 over a range of mapping quality cutoffs. This difference was much bigger for unpaired than paired-end reads. Nevertheless, in paired-end, BWA sometimes switched to local alignment whereas Bowtie2 did not.

For longer read length accuracy, Bowtie2 outperformed BWA-SW especially for average length of 250 nt. Moreover, Bowtie2 also trimmed less reads than BWA-SW (Langmead & L. Salzberg, 2012).

## Latest Algorithms

- Spark cluster parallelization of HISAT2 (Guo et al., 2022): efficiency is increased with more reads (tasks)
- Accel-Align (Yan et al., 2021): hash table, no support for local alignment, way faster than Bowtie2 and BWA, could be even used for error-prone long reads and parallelized in the seeding stage rather than with dynamic programming which has heavy data dependencies (future work prospects of authors)

## References

Burrows, M., & Wheeler, D. (1994). A block-sorting lossless data compression algorithm. In *Digital SRC Research Report.*

Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., ... & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, *29*(1), 15-21.

Farrar, M. (2007). Striped Smith–Waterman speeds database searches six times over other SIMD implementations. *Bioinformatics*, *23*(2), 156-161.

Guo, J., Gao, J., & Liu, Z. (2022). HISAT2 Parallelization Method Based on Spark Cluster. In *Journal of Physics: Conference Series* (Vol. 2179, No. 1, p. 012038). IOP Publishing.

Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature biotechnology*, *37*(8), 907-915.

Kim, D., Pertea, G., Trapnell, C., Pimentel, H., Kelley, R., & Salzberg, S. L. (2013). TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. *Genome biology*, *14*(4), 1-13.

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature methods*, *9*(4), 357-359. ( [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2) )

Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome biology*, *10*(3), 1-10.

Li, H., & Durbin, R. (2010). Fast and accurate long-read alignment with Burrows–Wheeler transform. *Bioinformatics*, *26*(5), 589-595.

Li, R., Yu, C., Li, Y., Lam, T. W., Yiu, S. M., Kristiansen, K., & Wang, J. (2009). SOAP2: an improved ultrafast tool for short read alignment. *Bioinformatics*, *25*(15), 1966-1967.

Marçais, G., Delcher, A. L., Phillippy, A. M., Coston, R., Salzberg, S. L., & Zimin, A. (2018). MUMmer4: A fast and versatile genome alignment system. *PLoS computational biology*, *14*(1), e1005944.

Musich, Ryan, Lance Cadle-Davidson, and Michael V. Osier. "Comparison of short-read sequence aligners indicates strengths and weaknesses for biologists to consider." *Frontiers in Plant Science* 12 (2021).

Yan, Y., Chaturvedi, N., & Appuswamy, R. (2021). Accel-align: a fast sequence mapper and aligner based on the seed–embed–extend method. *BMC bioinformatics*, *22*(1), 1-20.

# Bowtie2

bowtie2 can be used to:

    - Index reference FASTA nucleotide genomes/sequences
    - Align FASTQ sequencing reads to those genomes/sequences

## Differences between bowtie and bowtie2

    - bowtie2 has no upper limit on read length
    - bowtie2 can make gapped alignments
    - bowtie2 is more flexible for paired-end alignment
    - bowtie2 is faster and more memory efficient
    - bowtie is advantageous over bowtie2 for relatively short sequencing reads (50bp or less)


## Using Bowtie2 on the command line

The following DNA sequencing read data files were downloaded from the NCBI Sequence Read Archive. You can download it from [here](https://www.ebi.ac.uk/ena/browser/view/SRR030257?show=readsls).


### Indexing a reference genome/sequence using bowtie2-build

Before aligning reads to a reference genome with bowtie2, it must be indexed using bowtie2-build. This command will create six files with the extensions .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. These six files together are the index. Once an index has been created, the original reference genome/sequence is no longer needed to align reads. Here’s an example bowtie2-build command:

```
$ bowtie2-build reference_sequence.fasta index_name
```          

In this command, the reference_sequence.FASTA is the nucleotide FASTA sequence we want to index, and index_name is the name of the index. There will be six files beginning with the index_name in the output directory: index_name.1.bt2, index_name.2.bt2, index_name.3.bt2, index_name.4.bt2, index_name.rev.1.bt2, and index_name.rev.2.bt2. There’s no need to specify any of these files individually, just the index_name alone is enough to refer to the entire index.

### Aligning reads to an indexed genome/sequence using bowtie2

Now that the genome has been indexed, FASTQ sequencing reads can be aligned to it. This is done using the bowtie2 command. Here’s an example bowtie2 command:

```
$ bowtie2 --no-unal -p n -x index_name -1 reads_1.fastq -2 reads_2.fastq -S output.sam
```    

In this command…

`--no-unal` is an optional argument, meaning reads that do not align to the reference genome will not be written to sam output
`-p` is the number (n) of processors/threads used
`-x` is the genome index
`-1` is the file(s) containing mate 1 reads
`-2` is the file(s) containing mate 2 reads
`-S` is the output alignment in sam format

## Using Bowtie2 in R

### Installation 

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("Rbowtie2")
```

### Loading
```
library(Rbowtie2)
```

### Adapter Trimming
**Idetitify Adapter**

If you know the adapter sequence of reads files, you can skip this step. Besides,single end data is not support for this function yet so adapter sequence has to be known .

reads_1 and reads_2 are raw paired-end reads file with fastq format. adapters is two adapters character vector.
```
td <- tempdir()
reads_1 <- system.file(package="Rbowtie2", "extdata", "adrm", "reads_1.fq")
reads_2 <- system.file(package="Rbowtie2", "extdata", "adrm", "reads_2.fq")
(adapters <- 
    identify_adapters(file1=reads_1,file2=reads_2,
                      basename=file.path(td,"reads"),
                      "--threads 3",overwrite=TRUE))
```
**Remove Adapters**
```
(cmdout<-remove_adapters(file1=reads_1,file2=reads_2,adapter1 = adapters[1], 
                adapter2 = adapters[2],
output1=file.path(td,"reads_1.trimmed.fq"),
output2=file.path(td,"reads_2.trimmed.fq"),
basename=file.path(td,"reads.base"),overwrite=TRUE,"--threads 3"))
```
### Alignment

```
reads_1 <- system.file(package="Rbowtie2", "extdata", 
                       "bt2", "reads", "reads_1.fastq")
reads_2 <- system.file(package="Rbowtie2", "extdata", 
                       "bt2", "reads", "reads_2.fastq")
if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
    (cmdout<-bowtie2_samtools(bt2Index = file.path(td, "lambda_virus"),
        output = file.path(td, "result"),
        outputType = "sam",
        seq1=reads_1,
        seq2=reads_2,
        overwrite=TRUE,
        bamFile = NULL,
        "--threads 3"))
    head(readLines(file.path(td, "result.sam")))
}
```
# Differential-Expression-Analysis

We assume that there are fastq files (raw sequences) in our hand. If not, you can visit this https://www.ncbi.nlm.nih.gov/books/NBK242621/ website address to download fastq files from NCBI/SRA. For installation of SRA toolkit: https://ncbi.github.io/sra-tools/install_config.html.

# Quality control of fastq files (Step 1)

FASTQC tool is used for the quality control of fastq files. If base quality and adapter context are good, we can pass to further step, however, ıf not, bases with poor quality are removed from fastq files to reduce technical bias of RNA-seq analysis. To remove those bases, cutadapt is useful tool. We assume that our fastq files have a lot of adapter sequences (synthetic sequences) and so we want to remove them. For instance:

~$ cutadapt -a forward_adapter_sequence -A reverse_adapter_sequence -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq (paired-end trimming). In here, -a argument corresponds forward adapter sequences, -A argument corresponds reverse adapter sequences, -o argument corresponds output file (adapter trimmed fastq file), and -p argument corresponds paired-end mode. 
In the single-end trimming, ~$ cutadapt -a AACCGGTT -o output.fastq input.fastq or  cutadapt -g AACCGGTT -o output.fastq input.fastq. Cutadapt -a argument corresponds 3' adapter sequences, however -g corresponds 5' adapter sequences. Thus, overrepresent (contaminant) adapter sequences are removed by using cutadapt tool like examples above.

# Alignment Process (Step 2)

Each read in FASTQ files is aligned to reference genome to fit them to one or more genomic coordinates (namely, relevant genes). This step is achieved by several tools such as HISAT2, TopHat, so on. Refence genome is indexed by using aligners before aligment step. The purpose of reference genome index is sorting of genomic coordinates and this reduces alignment-process timing. We are utilizing HISAT2 aligner tool to make aligment process in here.

Reference genome indexing:

~$ hisat2-build reference_genome.fa index_name

(reference genome file should be fasta format (e.g: .fa or .fasta)

Aligment examples with HISAT2:

~$ hisat2 --dta -p --rna-strandness R -x reference_genome_indexes -U input.fastq -S output.sam 

The meaning of the arguments:

--dta argument: Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

-p argument: Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments.

--rna-strandness argument: Specify strand-specific information: the default is unstranded. For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF.

-x argument: The basename of the index for the reference genome.

-U argument: This argument is fitted for only single-end reads.

-S argument: This argument is used for generation of SAM files.

# The file conversion (Step 3)

In this step, SAM files are converted to BAM files by using Samtools. BAM files have less storage-load than SAM files and they can be thought as the compressed SAM files.

~$ samtools sort -@ 6 -o file.bam file.sam

The meaning of the arguments:

sort argument: SAM files are sorted by genomic coordinates or names. The sorting based on genomic coordinates or names is performed according to purpose of downstream analysis. 

-@ argument: The number of threads.

-o argument: This argument is used for generation of files with .bam extension.

# Gene-level quantification (Step 4)

All transcripts belonging to each gene are quantified across samples by using the featureCounts tool and this is made with bam files. After quantification, featureCounts generates an expression matrix in which each column represents individual sample, but each row represents individual gene.

~$ featureCounts patient1.bam patient2.bam control1.bam control2.bam -a annotation_file.gtf -o names_of_output_file -g gene_id -T 6 -s 2 -Q 50 --verbose

The meaning of arguments:

-a argument: This argument is required for using of gene annotation file (.GTF). According to annotation file, each transcript is assigned to its relevant genes.

-o argument: Output file is generated by this argument. Output is an expression matrix in which rows represent each gene, but columns represent each sample.

-g argument: This argument is used for gene-level quantification.

-T: This argument determines the number of threads.

-s argument: This argument assings strandness information. If it is used as the wrong then each transcript is quantified for wrong a gene and this effects downstream analysis results.

-Q argument: This argument filters transcripts with low mapping quality score. If mapping quality is low, it means that the same transcipt is aligned to at least two (or more) genomics coordinates (namely, genes) and this contributes to wrong gene-level quantification.

--verbose argument: Output verbose information for debugging, such as un-matched chromosome/contig names. During running analysis, information about process is simultaneously showed in the same window.

After gene-level quantification, the expression matrix is look like below:


|   | Control_1 | Control_2  | Patient_1  | Patient_2  |  
| ------------- | ------------- | ------------- | -------------  | ------------- 
| Content Cell  | Content Cell  | Content Cell  |  Content Cell  | Content Cell
| Content Cell  | Content Cell  | Content Cell  |  Content Cell  | Content Cell
















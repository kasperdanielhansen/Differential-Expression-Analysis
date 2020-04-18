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

-S argument: 





















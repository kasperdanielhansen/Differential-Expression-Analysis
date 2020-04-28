

# Differential-Expression-Analysis

We assume that we have FASTQ files (raw sequences). If not, you can visit this https://www.ncbi.nlm.nih.gov/books/NBK242621/ website address to download FASTQ files from NCBI/SRA. For installation of SRA toolkit: https://ncbi.github.io/sra-tools/install_config.html.

# Quality control of FASTQ files (Step 1)

FASTQC tool is used for the quality control of FASTQ files. If base quality and adapter context are good, we can move on next step, however, ıf not, bases with poor quality are removed to reduce technical bias of RNA-seq analysis. To remove those bases, cutadapt tool is used. We assume that our FASTQ files have a lot of adapter sequences (synthetic sequences) and so we want to remove them.

------------------------------------------------------------------------------------------------------

For paired-end FASTQ files:

~$ cutadapt -a forward_adapter_sequence -A reverse_adapter_sequence -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq 

The meanings of arguments:

-a argument: This argument corresponds to forward adapter sequences.
-A argument: This argument corresponds to reverse adapter sequences.
-o argument: This argument corresponds to output file (adapter trimmed fastq file)
-p argument: This argument corresponds to paired-end mode.

reads.1.fastq and reads.2.fastq file are input paired-end files, but out.1.fastq and out.2.fastq files are output files (trimmed fastq files).

Adapter sequences such as  ATCCGGAT and TTCGATCGA can be written for -a and -A arguments, respectively.

------------------------------------------------------------------------------------------------------

For single-end FASTQ files:

~$ cutadapt -a AACCGGTT -o output.fastq input.fastq or cutadapt -g AACCGGTT -o output.fastq input.fastq. 

The meanings of arguments:
-a argument corresponds to 3' adapter sequences.
-g argument: corresponds to 5' adapter sequences. 

Thus, overrepresent (contaminant) adapter sequences are removed by using cutadapt tool for both single and paired-end FASTQ files.

# Alignment process (Step 2)

In this step, trimmed (cleaned) FASTQ files are used. Each read in FASTQ files is aligned to reference genome to fit them to one or more genomic coordinates (namely, relevant genes). This step is achieved by several tools such as HISAT2, TopHat, so on. Reference genome is indexed by using aligners before alignment step. The purpose of reference genome index is sorting of genomic coordinates and this reduces alignment-process timing. We are utilizing HISAT2 aligner tool for index and alignment processes here.

Reference genome indexing:

~$ hisat2-build reference_genome.fa index_name

(reference genome file should be fasta format (e.g: .fa or .fasta)

Aligment examples with HISAT2:

~$ hisat2 --dta -p 4 --rna-strandness R -x reference_genome_indexes -U input.fastq -S output.sam 

------------------------------------------------------------------------------------------------------

The meaning of the arguments:

--dta argument: Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

-p argument: Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. We used four cores here.

--rna-strandness argument: Specify strand-specific information: the default is unstranded. For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF. We used 'R' argument for first-strand format FASTQ files, which means that our reads correspond to reverse complemented of transcripts. If our FASTQ files are second-strand format, then we should use 'F' argument instead of 'R' argument.

-x argument: The basename of the index for the reference genome.

-U argument: This argument is fitted for only single-end reads.

-S argument: This argument is used for generation of SAM files.

# File conversion (Step 3)

In this step, SAM files are converted to BAM files by using Samtools. BAM files have less storage-load than SAM files and they can be thought as the compressed SAM files.

~$ samtools sort -@ 6 -o file.bam file.sam

------------------------------------------------------------------------------------------------------

The meaning of the arguments:

sort argument: SAM files are sorted by genomic coordinates or names. The sorting based on genomic coordinates or names is performed according to purpose of downstream analysis. 

-@ argument: The number of threads. We used six cores here.

-o argument: This argument is used for generation of files with .bam extension.

# Gene-level quantification (Step 4)

All transcripts belonging to each gene are quantified across samples by using the featureCounts tool and this is made with bam files. After quantification, featureCounts generate an expression matrix in which each column represents individual sample, but each row represents individual gene.

~$ featureCounts control1.bam control2.bam patient1.bam patient2.bam -a annotation_file.gtf -o names_of_output_file -g gene_id -T 6 -s 2 -Q 50 --verbose

------------------------------------------------------------------------------------------------------

The meaning of arguments:

-a argument: This argument is required for using of gene annotation file (.GTF). According to annotation file, each transcript is assigned to its relevant genes.

-o argument: Output file is generated by this argument. Output is an expression matrix in which rows represent each gene, but columns represent each sample.

-g argument: This argument is used for gene-level quantification.

-T: This argument determines the number of threads. We used six core here again.

-s argument: This argument assings strandness information. If it is used as the wrong then transcripts can be assigned to wrong a gene and this effects downstream analysis results. We selected reversely stranded option (namely, -s 2) because we have first-stranded FASTQ files.

-Q argument: This argument filters transcripts with low mapping quality score. If mapping quality is low, it means that the same transcipt may be aligned to at least two (or more) genomics coordinates (namely, genes) and this contributes to wrong gene-level quantification.

--verbose argument: Output verbose information for debugging, such as un-matched chromosome/contig names. During running analysis, information about process is simultaneously shown in the same window in which analysis run.

------------------------------------------------------------------------------------------------------

After gene-level quantification, the expression matrix looks like below:


|   | Control_1 | Control_2  | Patient_1  | Patient_2  |  
| ------------- | ------------- | ------------- | -------------  | ------------- 
| Gene_1  | 3  | 5  |  7  | 9
| Gene_2  | 30  | 40  |  50  | 60
| Gene_3  | 100  | 120  |  130  | 95

# Differential expression analysis with DESeq2 (Step 5)

After finishing of gene-level quantification, expression matrix is read in R statistical computing environment to carry out differential expression analysis steps. Differential expression workflow is shown below:

featureCounts_expression_matrix=read.delim("featureCounts_expression_matrix", header = T, row.names = 1) # reading
expression matrix in R.

head(featureCounts_expression_matrix) # Looking at first six rows of expression matrix.

------------------------------------------------------------------------------------------------------

expression_matrix_assignment=as.matrix(featureCounts_expression_matrix) # DESeq2 recognizes matrix data instead of list or data.frame, so our expression data must be assigned as a matrix.

class(expression_matrix_assignment) # Controlling type of our expression data (Its matrix or not?)
typeof(expression_matrix_assignment) # Controlling type of our expression data (Its matrix or not?) If data is matrix then
we can move on next steps.

------------------------------------------------------------------------------------------------------

storage.mode(expression_matrix_assignment)="numeric" # DESeq2 recognizes matrices at numeric-type, so our matrix must be assigned as the numeric.

apply(expression_matrix_assignment, 2, typeof) # Controlling type of our expression data (Its numeric or not?)

If our expression data is both matrix format and numeric, we can move on other steps.

------------------------------------------------------------------------------------------------------

groups <- factor(c(rep("Controls",2),rep("Patients",2))) # Each sample group is assigned as the factor to make trait (phenotype) annotation. The first two columns of expression data are control groups, but last two columns are patient groups.

head(groups) # Controlling sample annotation. This step is very critical because If groups are assigned as the wrong, then differentially expressed gene results between groups could be wrong.

------------------------------------------------------------------------------------------------------

Filtration of genes with low expression:

min_read <- 1 
filtered_expression_data <- expression_matrix_assignment[apply(expression_matrix_assignment,1,function(x){max(x)}) > min_read,] # Row-wise filtration is performed in expression data to eliminate genes with low expression across samples. Low expressed genes may cause statistical noise and create a bias in results.

------------------------------------------------------------------------------------------------------

sampleInfo <- data.frame(groups,row.names=colnames(filtered_expression_data)) # Sample-trait annotation step.

------------------------------------------------------------------------------------------------------

Differential expression analysis:

dds <- DESeqDataSetFromMatrix(countData = filtered_expression_data, colData = sampleInfo, design = ~ groups)

dds$groups = relevel(dds$groups,"Controls")

dds <- DESeq(dds)

res <- results(dds,independentFiltering=F) # In this step, statistical analysis results are avaiable for all genes, but statistically meaningful and non-meaningful results are together.

------------------------------------------------------------------------------------------------------

To select only statistically meaningful results:

resSig <- res[(!is.na(res$padj) & (res$padj <= 0.05) & (abs(res$log2FoldChange)>= 1.5)), ] # Here, genes were selected according to adjusted p-value <=0.05 and log2FoldChange >=1.5. As the results, we obtained both statistically and biologically meaningful differentially expressed genes between patients and controls.

------------------------------------------------------------------------------------------------------

To write results to computer from R statistical computing environment:

write.csv(resSig, file="differential_expression_results_significant_genes.csv", quote=F)

------------------------------------------------------------------------------------------------------

Obtained file with .csv extension looks like below:


|   | baseMean | log2FoldChange  | IfcSE  | stat  | pvalue | padj 
| ------------- | ------------- | ------------- | -------------  | ------------- | ------------- | -------------
| Gene_1  | 85.8220655465971  | -1.874026295362  |  0.266611865022727  | -7.02904311930098  | 0.00000000000207954712246367  | 0.0000000000755367230374539 
| Gene_2  | 98.2824885614715  | 2.15753251459839  |  0.60573737089831  | 3.56182830753658  | 0.00036828120481343  | 0.00214459983313731 
| Gene_3  | 412.438051846698  | 2.18571599360682  |  0.231625110250678  | 6.86459918229613  | 0.00000000000666782010218906   | 0.000000000218449038969323 

If you desire more detailed workflow, can visit https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html website.

If you have any questions about this pipeline, you feel free to contact me via kaan.okay@msfr.ibg.edu.tr e-mail address.

Thanks for your interest!

Best wishes :)

-------------------------------------------------------- The End --------------------------------------------------------


# Differential-Expression-Analysis

We assume that there are fastq files (raw sequences) in our hand. If not, you can visit this https://www.ncbi.nlm.nih.gov/books/NBK242621/ address to download fastq files from NCBI/SRA. For installation of SRA toolkit: https://ncbi.github.io/sra-tools/install_config.html.

# Quality control of fastq files (Step 1)

FASTQC tool is used for the quality control of fastq files. If base quality and adapter context are good, we can pass to further step, however, ıf not, bases with poor quality are removed from fastq files to reduce technical bias of RNA-seq analysis. To remove those bases, cutadapt is useful tool. We assume that our fastq files have a lot of adapter sequences (synthetic sequences) and so we want to remove them. For instance:

~$ cutadapt -a forward_adapter_sequence -A reverse_adapter_sequence -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq (paired-end trimming). In here, -a argument corresponds forward adapter sequences, -A argument corresponds reverse adapter sequences, -o argument corresponds output file (adapter trimmed fastq file), and -p argument corresponds paired-end mode. 
In the single-end trimming, ~$ cutadapt -a AACCGGTT -o output.fastq input.fastq or  cutadapt -g AACCGGTT -o output.fastq input.fastq. Cutadapt -a argument corresponds 3' adapter sequences, however -g corresponds 5' adapter sequences. Thus, overrepresent (contaminant) adapter sequences are removed by using cutadapt tool like examples above.

# RNA-seq-analyzer
This repository contains a script (app.R) which has the code used to creat the shinny app "RNA-seq analyzer" (https://thaisharder.shinyapps.io/RNAseq/). The app was designed to analyze and organize RNA-seq raw data. The app also creates a MA plot and an interactive table, which helps to visualize your reuslts quickly and more efficiently.
To use the app it is necessary to have your owm RNAseq data. Once you have received the raw RNA-seq data as a FASTQ-format file containing the reads sequenced from the Next Generation Sequencing platform, you will align these reads to an annotated reference genome, and quantify the expression of genes. 

First step: building your index. Go to Bowtie's web site http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
Use the Command: bowtie2-build <reference_in> <bt2_base>

`<reference _in>` is a comma-separated list of FASTA files containing the reference sequences to be aligned to, and `<bt2_base>` is the basename of the index files to write.

Sencond step: aligning the data from the FASTAQ files with the genome which is the index that you built.

Use the command: bowtie2 -x ~/Directory/index -1 .ForwardSequence.fastq -2 ./ReverseSequence.fastq -S OutputName.sam -p 2

Third step: counting the genes that are being expressed using htseq-count and generate a table with your data.

Use the command: htseq-count -s no -t CDS -a 1 -m intersection-nonempty -i ID OutputName.sam HP_GCF_000210895.1_ASM21089v1_genomic.gff > OutputName.txt


You can also find these intructions on the app. Please see the tab named "Instructions" that explains how to use Bowtie which is a short read aligner and HtSeq to creat your table. 
After generating your table you can upload it on the "Start" tab and vizualize your result. 

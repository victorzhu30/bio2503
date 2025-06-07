# bio2503
Linux操作系统与Shell开发的理论与实践

## Script Overview

```bash
#!/bin/bash

# This script performs genome assembly using Illumina, PacBio, and Nanopore data.
# It includes steps for data download, quality control, k-mer analysis,
# assembly with different tools, and evaluation of assemblies.

# === 0. INITIAL SETUP AND DIRECTORY CREATION ===
echo "INFO: Starting script and setting up directories."
mkdir -p ~/bio2503/05.assembly/{data,log} # Use -p to create parent dirs if they don't exist and suppress errors if they do
cd ~/bio2503/05.assembly/data

# === 1. CONDA ENVIRONMENT SETUP AND TOOL INSTALLATION ===
echo "INFO: Setting up Conda environment 'assembly' and installing tools."
conda create -n assembly -y
conda activate assembly
conda install -c bioconda -c conda-forge sra-tools fastqc fastp jellyfish genomescope2 spades soapdenovo2 soapdenovo2-gapcloser seqkit seqtk quast flye

# === 2. DATA DOWNLOAD ===
prefetch SRR8482586 > ../log/SRR8482586.log 2>&1 
prefetch SRR8494912 > ../log/SRR8494912.log 2>&1 
prefetch SRR8494939 > ../log/SRR8494939.log 2>&1 
efetch -db nuccore -format fasta -id NC_009648.1 > MGH78578.fasta
efetch -db nuccore -format fasta -id NC_009649.1 > MGH78578_pKPN3.fasta
efetch -db nuccore -format fasta -id NC_009650.1 > MGH78578_pKPN4.fasta
efetch -db nuccore -format fasta -id NC_009651.1 > MGH78578_pKPN5.fasta
efetch -db nuccore -format fasta -id NC_009652.1 > MGH78578_pKPN6.fasta
efetch -db nuccore -format fasta -id NC_009653.1 > MGH78578_pKPN7.fasta

echo "INFO: Organizing downloaded SRA files."
mv SRR8482586/SRR8482586.sra ./
mv SRR8494912/SRR8494912.sra ./
mv SRR8494939/SRR8494939.sra ./
rm -r SRR8482586
rm -r SRR8494912
rm -r SRR8494939
mv SRR8482586.sra illumina.sra
mv SRR8494912.sra pacbio.sra
mv SRR8494939.sra nanopore.sra

# === 3.CLUSTER JOB SUBMISSION AND ENVIRONMENT REACTIVATION ===
qsub -I -q cpuq -l nodes=2:ppn=8
conda activate assembly
cd ~/bio2503/05.assembly/data

# === 4. SRA TO FASTQ CONVERSION AND COMPRESSION ===
echo "INFO: Converting SRA files to FASTQ format using fasterq-dump."
fasterq-dump illumina.sra
fasterq-dump pacbio.sra
fasterq-dump nanopore.sra
gzip *.fastq
mkdir qc

# === 5. ILLUMINA READ QUALITY CONTROL AND FILTERING ===
echo "INFO: Performing QC and filtering on Illumina reads."
fastqc -f fastq -o qc illumina_1.fastq.gz illumina_2.fastq.gz
fastp -i illumina_1.fastq.gz -I illumina_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -D -z 4 -q 20 -u 30 -n 10 -f 20 -t 10 -F 20 -T 10 -h clean.html

# === 6. K-MER ANALYSIS FOR GENOME CHARACTERISTICS ESTIMATION (ILLUMINA) ===
echo "INFO: Performing k-mer analysis using Jellyfish and GenomeScope2."
mkdir ../34.kmer
cd ../34.kmer/
echo "INFO: Counting k-mers (k=15) from clean.1.fq.gz."
jellyfish count -m 15 -s 2G -o kmer15.count -C <(zcat ../data/clean.1.fq.gz) # 不支持压缩格式
jellyfish stats kmer15.count -o kmer15.stats
jellyfish histo kmer15.count > kmer15.histo
genomescope2 -i kmer15.histo -o gscope15 -p 1 -k 15
# seqkit stats ../data/MGH78578.fasta
echo "INFO: Counting k-mers (k=17) from clean.2.fq.gz."
jellyfish count -m 17 -s 2G -o kmer17.count -C <(zcat ../data/clean.2.fq.gz)
jellyfish stats kmer17.count -o kmer17.stats
jellyfish histo kmer17.count > kmer17.histo
genomescope2 -i kmer17.histo -o gscope -p 1 -k 17
echo "INFO: Counting k-mers (k=21) from clean.1.fq.gz."
jellyfish count -m 21 -s 2G -o kmer21.count -C <(zcat ../data/clean.1.fq.gz)
jellyfish stats kmer21.count -o kmer21.stats
jellyfish histo kmer21.count > kmer21.histo
genomescope2 -i kmer21.histo -o gscope -p 1 -k 21

# === 7. ILLUMINA-ONLY ASSEMBLY (SOAPDENOVO2) ===
mkdir ../36.illumina
cd ../36.illumina
echo "INFO: Uncompressing cleaned Illumina reads for SOAPdenovo2."
gzip -d -c ../data/clean.1.fq.gz > ../data/clean.1.fastq
gzip -d -c ../data/clean.2.fq.gz > ../data/clean.2.fastq
echo -e "max_rd_len=150\n\n[LIB]\navg_ins=439\nreverse_seq=0\nasm_flags=3\nrank=1\npair_num_cutoff=3\nq1=/full/path/clean.1.fastq\nq2=/full/path/clean.2.fastq" > soapdenovo2.config
mkdir kmer45 kmer85
echo "INFO: Running SOAPdenovo2 with K=45."
SOAPdenovo-63mer all -s soapdenovo2.config -K 45 -o kmer45/kmer45 -D 1 -d 1 -u 2 -p 16 > kmer45/kmer45.log 2>&1
echo "INFO: Running SOAPdenovo2 with K=85."
SOAPdenovo-127mer all -s soapdenovo2.config -K 85 -o kmer85/kmer85 -D 1 -d 1 -u 2 -p 16 > kmer85/kmer85.log 2>&1
echo "INFO: Running GapCloser for K=45 assembly."
GapCloser -a kmer45/kmer85.scafSeq -b soapdenovo2.config -o kmer45/kmer4
5.fill.fa -t 12
echo "INFO: Running GapCloser for K=85 assembly."
GapCloser -a kmer85/kmer85.scafSeq -b soapdenovo2.config -o kmer85/kmer8
5.fill.fa -t 12

# === 8. EVALUATION OF ILLUMINA ASSEMBLIES ===
echo "INFO: Evaluating Illumina assemblies."
mkdir ../37.evaluate
cd ../37.evaluate
cp ../36.illumina/kmer45/kmer45.scafSeq .
# seqkit stats kmer45.scafSeq
seqkit fx2tab kmer45.scafSeq -o kmer45.scafSeq.table
awk '{print $1 "\t" length($3)}' kmer45.scafSeq.table | sort -nr -k2 | head
seqtk comp kmer45.scafSeq

cp ../data/MGH78578.fasta .
cp ../36.illumina/kmer85/kmer85.scafSeq .
mv kmer45.scafSeq soapdenovo2_kmer45.fasta
mv kmer85.scafSeq soapdenovo2_kmer85.fasta
mkdir quast
echo "INFO: Running QUAST with reference for Illumina assemblies."
quast -o quast -r MGH78578.fasta -t 12 soapdenovo2_kmer45.fasta soapdenovo2_kmer85.fasta > quast/quast.log 2>&1

# === 9. LONG-READ ASSEMBLY (FLYE - PACBIO AND RAW NANOPORE) ===
mkdir ../40.nanopore/
cd ../40.nanopore/
mkdir flye
echo "INFO: Assembling PacBio reads with Flye."
flye --pacbio-raw ../data/pacbio.fastq.gz -g 5.5m -t 12 -o flye_pacbio > flye_pacbio.log 2>&1
echo "INFO: Assembling raw Nanopore reads with Flye."
flye --nano-raw ../data/nanopore.fastq.gz -g 5.5m -t 12 -o flye_nanopore > flye_nanopore.log 2>&1

# === 10. NANOPORE READ QC AND FILTERING ===
echo "INFO: Performing QC and filtering on Nanopore reads."
conda deactivate
conda create -n nanopore -c bioconda -c conda-forge nanoplot filtlong -y
conda activate nanopore
echo "INFO: Generating NanoPlot QC for raw Nanopore reads."
NanoPlot --fastq ../data/nanopore.fastq.gz -o nanoplot -t 12
echo "INFO: Filtering Nanopore reads with Filtlong."
filtlong --min_length 1000 --min_mean_q 80 ../data/nanopore.fastq.gz > clean.filter.fastq
echo "INFO: Generating NanoPlot QC for filtered Nanopore reads."
NanoPlot --fastq clean.filter.fastq -o clean -t 12

# === 11. RE-ASSEMBLY OF FILTERED NANOPORE READS ===
echo "INFO: Re-assembling filtered Nanopore reads with Flye."
conda deactivate
conda activate assembly
flye --nano-raw clean.filter.fastq -g 5.5m -t 12 -o flye_nanopore_clean > flye_nanopore_clean.log 2>&1
# === 12. EVALUATION OF LONG-READ ASSEMBLIES ===
echo "INFO: Evaluating long-read assemblies with QUAST."
mkdir quast
cd quast/
cp ../flye_pacbio/assembly.fasta pacbio.fasta
cp ../flye_nanopore/assembly.fasta nanopore.fasta
cp ../flye_nanopore_clean/assembly.fasta nanopore_clean.fasta
cp ../../data/MGH78578.fasta .
quast -o quast -r MGH78578.fasta -t 12 nanopore.fasta nanopore_clean.fasta pacbio.fasta > quast.log 2>&1
echo "INFO: Script finished!!!"
```
## Usage

1. Save the script to a file, e.g., `assembly_pipeline.sh`.
2. Make the script executable:
   ```bash
   chmod +x assembly_pipeline.sh
   ```
3. Run the script:
   ```bash
    ./assembly_pipeline.sh
    ```

## Advice

- Although the script is written, it is still recommended to run it line by line, observing the output and log files at each step to ensure no errors occur.
- Ensure you have sufficient disk space and memory available, as genome assembly can be resource-intensive.
- Adjust the number of threads (`-t`) based on your system's capabilities.
- Check the paths in the script to ensure they match your directory structure.

## References
[1] https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000294

[2] https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/

[3] Guillaume Marcais and Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics (2011) 27(6): 764-770

[4] [aquaskyline/SOAPdenovo2: Next generation sequencing reads de novo assembler.](https://github.com/aquaskyline/SOAPdenovo2)

[5] https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/de_novo_assembly_tools/soapdenovo2/

[6] Ranallo-Benavidez, T.R., Jaron, K.S. & Schatz, M.C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* **11**, 1432 (2020)

[7] https://vcru.wisc.edu/simonlab/bioinformatics/programs/soap/GapCloser_Manual.pdf

[8] https://www.biostars.org/p/95803/#9496314

[9] Alla Mikheenko, Vladislav Saveliev, Pascal Hirsch, Alexey Gurevich, WebQUAST: online evaluation of genome assemblies, *Nucleic Acids Research* (2023) 51 (W1): W601–W606. doi: [10.1093/nar/gkad406](http://dx.doi.org/10.1093/nar/gkad406)

[10] Vurture, GW, Sedlazeck, FJ, Nattestad, M, Underwood, CJ, Fang, H, Gurtowski, J, Schatz, MC (2017) *Bioinformatics* doi: https://doi.org/10.1093/bioinformatics/btx153

[11] [[算法学习 1\] 基因组组装算法 De Bruijn Graph - 知乎](https://zhuanlan.zhihu.com/p/57177938)

[12] [kmer估计基因组大小](https://www.dxy.cn/bbs/newweb/pc/post/46489294)

[13] [MSELab](http://47.57.89.98/detail/25/)

[14] [Meta组装小课题2:GapCloser - 简书](https://www.jianshu.com/p/8cf8b961afbf)

[15] Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. *iMeta* e191. [doi:10.1002/imt2.191](https://doi.org/10.1002/imt2.191).

[16] [lh3/seqtk: Toolkit for processing sequences in FASTA/Q formats](https://github.com/lh3/seqtk)

[17] [rrwick/Filtlong: quality filtering tool for long reads](https://github.com/rrwick/Filtlong)

[18] Ruan, J. and Li, H. (2019) Fast and accurate long-read assembly with wtdbg2. *Nat Methods* doi:10.1038/s41592-019-0669-3

[19] Ruan, J. and Li, H. (2019) Fast and accurate long-read assembly with wtdbg2. *bioRxiv*. doi:10.1101/530972

[20] https://space.bilibili.com/404172354
# Create and install dependencies as conda env
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n rnaseq_demo python==3.10
conda activate rnaseq_demo
conda install -c bioconda samtools subread samtools

# create directories
mkdir mapping_data
mkdir ref_data

# run QC once for the chr1 subsetted reads
./tools/FastQC/fastqc raw_data/test1_R*.fq

java -jar /workspaces/bulkRNAseq_demo/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -trimlog ./raw_data/sam_new_NEXTERA_trimlog.txt \
    ./raw_data/test1_R1.fq  ./raw_data/test1_R1.fq \
   ./clean_data/sam_new_trimmed_R1.fastq.gz ./clean_data/sam_new_trimmed_R2.fastq.gz \
    LEADING:20 \
    ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:7:2:True \
    MINLEN:36

# download the chr1 references from UCSC and uncompress it
wget -P ref_data https://hgdownload.gi.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
wget -P ref_data https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
cat hg38.refGene.gtf | grep chr1 >> chr1_hg38.gtf # gtf file for chr1 alone

## build index and map the reads, add sample id and read groups for easy access
hisat2-build -p 4 chr1.fa chr1_hisat2
hisat2 -p 4 \
  -x ./ref_data/chr1_hisat2 \
  -1 ./clean_data/sam_new_trimmed_R1.fastq.gz \
  -2 ./clean_data/sam_new_trimmed_R2.fastq.gz \
  -S ./mapping_data/sam_new.sam \
  --summary-file ./mapping_data/hisat2_summary.txt \
  --rg-id sam1 \
  --rg SM:sam1 \
  --rg PL:ILLUMINA \
  --rg LB:lib1 \
  --rg PU:unit1

# sort the chr1_full sam file
# sort by read names with max mem at 2GB and 2 threads
samtools addreplacerg \
  -r ID:sam1 \
  -r SM:sam1 \
  -r PL:ILLUMINA \
  -r LB:lib1 \
  -r PU:unit1 \
  -o ./mapping_data/chr1_demo_rg.sam \
  ./mapping_data/chr1_demo.sam
samtools sort -n --threads 2 -m 2G -o ./mapping_data/chr1_sorted_demo.sam ./mapping_data/chr1_demo_rg.sam

# determine strandedness
# create bed file from gtf file
 awk 'BEGIN{OFS="\t"} 
$3=="exon" {
    gene=$10
    gsub(/"/,"",gene)
    gsub(/;/,"",gene)
    print $1,$4-1,$5,gene,".",$7
}' ./ref_data/chr1_hg38.gtf > ./ref_data/chr1_exons.bed

# remove chr tag in bed file
sed 's/^chr//' ./ref_data/chr1_exons.bed > ./ref_data/chr1_exons_nochr.bed
infer_experiment.py -r ./ref_data/chr1_exons_nochr.bed -i ./mapping_data/chr1_sorted_demo.sam

# output is reverse stranded, use 2 in featureCounts -s
featureCounts -F GTF -t exon -g gene_id -s 2 \
-Q 30 -p -P -B --ignoreDup -T 4 -R SAM --verbose \
-a ./ref_data/chr1_hg38.gtf \
-o sam1_new_counts \
./mapping_data/chr1_sorted_demo.sam

# sam1_new_counts has the counts matrix for sample
awk 'BEGIN{OFS="\t"} 
NR==1{print; next} 
NR==2{$7="sam1"} 
{print}' sam1_new_counts > sam1_counts_clean.txt
awk -F '\t' '{print $1,$7}' sam1_counts_clean.txt | head
### Sample SRA: ERR12917628	
### Paired end Illumina NovaSeq X

# download basic tools
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

# get raw data
../tools/sratoolkit.3.3.0-ubuntu64/bin/fastq-dump --split-files --gzip -X 100000 ERR12917628

# fastqc pre

# trimmomatic
java -jar trimmomatic-0.39.jar PE -threads 4 -trimlog ../../raw_data/sam1_NEXTERA_trimlog.txt \
    ../../raw_data/ERR12917628_1.fastq.gz  ../../raw_data/ERR12917628_2.fastq.gz \
    -baseout ../../clean_data/sam1_trimmed_nextera \
    LEADING:20 \
    ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:7:2:True \
    MINLEN:36

# fastqc post
./fastqc ../../clean_data/sam1_trimmed_nextera_*
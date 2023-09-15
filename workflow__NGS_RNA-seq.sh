# Fastqc the 6 files
fastqc -t 4 *fastq -o ../fastqc_output

# Trimmomatic --> Do it for every fastq file
 trimmomatic PE -threads 4 
    ../data_chr19/SRR10045016_1.fastq 
    ../data_chr19/SRR10045016_2.fastq 
    ../fastqc_trimmed_output/SRR10045016_1P.fastq 
    ../fastqc_trimmed_output/SRR10045016_1U.fastq 
    ../fastqc_trimmed_output/SRR10045016_2P.fastq 
    ../fastqc_trimmed_output/SRR10045016_2U.fastq 
    ILLUMINACLIP:trimmomaticAdapters/TruSeq3-PE-2.fa:2:30:10 TRAILING:30


# Create the genome index
 STAR --runThreadN 2 --runMode genomeGenerate --genomeDir STAR_index_chr19/ --genomeFastaFiles genome/chr19.fa --sjdbGTFfile genome/chr19_Homo_sapiens.GRCh38.95.gtf
        STAR --runThreadN 2 --runMode genomeGenerate --genomeDir STAR_index_chr19/ --genomeFastaFiles genome/chr19.fa --sjdbGTFfile genome/chr19_Homo_sapiens.GRCh38.95.gtf


# Create subfolders for the STAR's mapping outputs
for i in {16..21}; do
    mkdir "SRR100450$i"
done                      


###  NOTE: STAR's mapping may take ONE DAY to finish. It depents from the  ###
###  size of the files and also the genome. It is advisable to create a    ###
###  tmux-session first and then create a .sh file, running                ###
###  a for loop.                                                           ###
###  If not using a cluster then don't shut down your laptop               ###

# Mapping with STAR(just the first sequence)
STAR --runThreadN 8 --genomeDir STAR_index_chr19/ --readFilesIn trim_output/SRR10045016_1P.fastq trim_output/SRR10045016_2P.fastq --sjdbGTFfile genome/chr19_Homo_sapiens.GRCh38.95.gtf --outFileNamePrefix star_output/SRR10045016/


# Now create a for loop for all the rest of the sequences
for i in {17..21}; do
        STAR --runThreadN 5 --genomeDir STAR_index_chr19/ --readFilesIn trim_output/"SRR100450$i"_1P.fastq trim_output/"SRR100450$i"_2P.fastq --sjdbGTFfile genome/chr19_Homo_sapiens.GRCh38.95.gtf --outFileNamePrefix STAR_output/"SRR100450$i"/
done


# Given a file with aligned sequencing reads 
# and a list of genomic features, 
# a common task is to count how many reads map to each feature.
htseq-count STAR_output/SRR10045016/Aligned.out.sam STAR_output/SRR10045017/Aligned.out.sam STAR_output/SRR10045018/Aligned.out.sam STAR_output/SRR10045019/Aligned.out.sam STAR_output/SRR10045020/Aligned.out.sam STAR_output/SRR10045021/Aligned.out.sam genome/chr19_Homo_sapiens.GRCh38.95.gtf > htseq_output/SRR10045016-17-18-19-20-21_counts.csv


#!/bin/bash
#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=48:00:00
#PBS -N your_title
#PBS -d /path/to/working/directory
#PBS -m abe
#PBS -M your.email@domain.com
#PBS -q your_queue

# Define paths (update these placeholders before use)
fqpath="/path/to/raw_fastq_files"
linkpath="/path/to/working/directory/data"
adapter="/path/to/illumina_adapter.fa"
hs2index="/path/to/hisat2_index/chr2"
gtf="/path/to/annotations/chr2.gtf"

mkdir -p "${linkpath}"

# Create symbolic links for raw fastq files and generate test fastq files with 25K reads
for sample in s1.c2 s2.c2 s3.c2 s4.c2 s5.c2 s6.c2 s7.c2 s8.c2 s9.c2 s10.c2 s11.c2 s12.c2
do
    raw_file="${fqpath}/${sample}.fq"
    link_name="${linkpath}/${sample}.fq"
    ln -s "${raw_file}" "${link_name}"
    test_file="${sample}.25K.fq"  # Fastq with 25K reads for testing
    head -n 100000 "${link_name}" > "${test_file}"
done

# Create output directories
hisat_dir="./hisat"
stringtie_dir="./stringtie"
mkdir -p "${hisat_dir}"
mkdir -p "${stringtie_dir}"

gtflist="list.gtf.txt"
rm -f "${gtflist}"

# Loop over samples to trim, align, and assemble transcripts
for sample in s1.c2 s2.c2 s3.c2 s4.c2 s5.c2 s6.c2 s7.c2 s8.c2 s9.c2 s10.c2 s11.c2 s12.c2
do
    fastq="${linkpath}/${sample}.fq"
    trim1="${linkpath}/${sample}.t1.fq"
    trim2="${linkpath}/${sample}.t2.fq"
    sam="${hisat_dir}/${sample}.sam"
    bam="${hisat_dir}/${sample}.bam"
    sorted_bam="${hisat_dir}/${sample}.sort.bam"

    # Adapter trimming with Scythe
    scythe -q sanger -a "${adapter}" -o "${trim1}" "${fastq}"

    # Quality trimming with Sickle
    sickle se -f "${trim1}" -t sanger -o "${trim2}" -q 10 -l 50

    # Align reads with HISAT2
    hisat2 -p 4 --rna-strandness RF --phred33 -x "${hs2index}" -U "${trim2}" -S "${sam}"

    # Convert SAM to BAM, sort, and clean intermediate files
    samtools view -b -o "${bam}" "${sam}"
    samtools sort -o "${sorted_bam}" "${bam}"
    rm -f "${sam}" "${bam}" "${trim1}" "${trim2}"

    # Run StringTie for transcript assembly
    str_sample_dir="${stringtie_dir}/${sample}"
    mkdir -p "${str_sample_dir}"
    sample_gtf="${str_sample_dir}/${sample}_transcripts.gtf"
    stringtie -p 4 -t -e -B -G "${gtf}" -o "${sample_gtf}" "${sorted_bam}"

    echo "${sample} ${sample_gtf}" >> "${gtflist}"
done

# Generate count matrices from assembled transcripts
prepDE -i "${gtflist}"

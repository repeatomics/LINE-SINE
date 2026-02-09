#!/usr/bin/env bash
set -euo pipefail

mkdir -p fastqc_reports
fastqc *.fastq.gz -o fastqc_reports --quiet

mkdir -p trimmed_fastq

for i in {19..24}; do
  trimmomatic PE -threads 2 \
    ${i}_S${i}_L001_R1_001.fastq.gz ${i}_S${i}_L001_R2_001.fastq.gz \
    trimmed_fastq/${i}_S${i}_L001_R1_paired.fastq.gz \
    trimmed_fastq/${i}_S${i}_L001_R1_unpaired.fastq.gz \
    trimmed_fastq/${i}_S${i}_L001_R2_paired.fastq.gz \
    trimmed_fastq/${i}_S${i}_L001_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}/TruSeq3-PE.fa:2:30:10:8:true \
    SLIDINGWINDOW:4:20 \
    MINLEN:30
done

mkdir -p fastqc_trimmed

fastqc trimmed_fastq/*_paired.fastq.gz \
       -o fastqc_trimmed \
       -t 4 \
       --quiet

bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

for i in {19..24}; do
  bwa mem -t 8 -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    trimmed_fastq/${i}_S${i}_L001_R1_paired.fastq.gz \
    trimmed_fastq/${i}_S${i}_L001_R2_paired.fastq.gz \
  | samtools view -b - \
  | samtools sort -@ 4 -o ${i}.sorted.bam
done

for i in {19..24}; do
	samtools index ${i}.sorted.bam
done	

#per-base coverage
for i in {19..24}; do
	samtools depth ${i}.sorted.bam \
	| awk '{print $1"\t"$2-1"\t"$2"\t"$3}' \
    > ${i}.covered.base.bed
done

for i in {19..24}; do
	bedtools merge \
    -i ${i}.covered.base.bed \
    > ${i}.covered.regions.bed
done

#Сколько регионов на каждую хромосому
awk '{count[$1]++} END {for (c in count) print c, count[c]}' \
19.covered.regions.bed | sort -k1,1

#Файл с посчитанным расстоянием между ампликонами по кадой хромосоме отдельно
sort -k1,1 -k2,2n 19.covered.regions.bed \
| awk '
$1==prev_chr {
  dist = $2 - prev_end;
  print $1, dist
}
{
  prev_chr=$1;
  prev_end=$3
}' > 19.inter_amplicon_distances.txt

#Создание таблицы для оценки плотности покрытия хромосом.
awk '
{
  chr=$1; d=$2;
  if (d > 0) {
    count[chr]++;
    sum[chr]+=d;
    if (!(chr in min) || d < min[chr]) min[chr]=d;
    if (d > max[chr]) max[chr]=d;
  }
}
END {
  printf "Chromosome\tN_intervals\tMean_distance(bp)\tMin(bp)\tMax(bp)\n";
  for (c in count)
    printf "%s\t%d\t%.1f\t%d\t%d\n", c, count[c], sum[c]/count[c], min[c], max[c];
}' 19.inter_amplicon_distances.txt | sort -k1,1 > 19.inter_amplicon_distances.tsv

#переведем в tsv, там две колонки название хромосомы и расстояние между ампликонами

mv 19.inter_amplicon_distances.txt 19.inter_amplicon_just_distances.tsv
#будет косяк с разделителем колонок, решается так, но уже написано решение в юпитере
#awk '{print $1 "\t" $2}' 19.inter_amplicon_just_distances.tsv \
#> 19.inter_amplicon_just_distances.fixed.tsv
#Делаем понятные имена двум колонкам

sed -i '1i Chromosome\tDistance_bp' 19.inter_amplicon_just_distances.tsv

#Разбираемся с глубинйо покрытия
#Создаем таблицу с с хромосомой, позицией и глубиной покрытия каждой позиции
samtools depth 19.sorted.bam \
| awk '$3 > 0 {print $1 "\t" $2 "\t" $3}' \
> 19.covered_positions_depth.tsv

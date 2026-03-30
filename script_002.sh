#!/usr/bin/env bash
set -euo pipefail

#режем адпатеры, праймеры не трогаем
fastp -i *R1*.fastq.gz -I *R2*.fastq.gz -o 19_R1.trim.fastq.gz -O 19_R2.trim.fastq.gz --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim_poly_g --cut_tail --qualified_quality_phred 20 --length_required 30

#мапим на чистый референс, индексируем
bwa mem reference.primary.fna 19_R1.trim.fastq.gz 19_R2.trim.fastq.gz | samtools sort -o aligned.sorted.bam  
samtools index aligned.sorted.bam

#tagged.sorted.bam создается через питон, 003_script.ipynb
#индексируем
samtools index tagged.sorted.bam

#tagged.sorted.bam создается через питон, где к ридам цепляем тэги прямой обратный праймер
#собираем только покрытые участки, оставляем колонку с глубиной покрытия
bedtools genomecov -ibam tagged.sorted.bam -bga | awk '$4 > 0' > continuous_coverage.bed

#rmsk.bed -Файл RepeatMasker, фильтруем его чтобы оставить только LINE и SINE
awk '$4 ~ /#LINE/ || $4 ~ /#SINE/' rmsk.bed > line_sine.bed

#участки, пересекающиеся с LINE/SINE
bedtools intersect -a continuous_coverage.bed -b line_sine.bed > covered_in_repeats.bed

#участки вне повторов
bedtools intersect -a continuous_coverage.bed -b line_sine.bed -v > de_novo_candidates.bed
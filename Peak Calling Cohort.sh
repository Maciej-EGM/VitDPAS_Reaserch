#!/bin/bash
set -euo pipefail
export LC_ALL=C

MIN=24   # peak must be present in at least MIN of the 47 participants (change if needed)

#Peak calling (output name = full sample name, as in your previous run)#
for bam in *.sorted.bam; do
    n=${bam%.sorted.bam}
    macs3 callpeak -t "$bam" -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n "$n" -B
done

wc -l *_peaks.narrowPeak

##Chromosome Lengh reference##
fetchChromSizes hg38 > hg38.chrom.sizes

for d in 0 1; do
    #Mergin Bedgraph files in condition order#
    bedtools unionbedg -i *d${d}_treat_pileup.bdg > Day${d}.bdg

    ##Mean of all samples (columns 4..NF) -> 4-column bedGraph, sorted##
    awk 'BEGIN{OFS="\t"} {s=0; for(i=4;i<=NF;i++) s+=$i; print $1,$2,$3,s/(NF-3)}' Day${d}.bdg \
        | sort -k1,1 -k2,2n > Day${d}.bedGraph

    ##Converting to BigWig##
    bedGraphToBigWig Day${d}.bedGraph hg38.chrom.sizes Day${d}.bw

    #Overlap Beetween Peaks: regions with a peak in >= MIN samples#
    for f in *d${d}_peaks.narrowPeak; do
        awk -v s="$f" 'BEGIN{OFS="\t"} {print $1,$2,$3,s}' "$f"
    done | sort -k1,1 -k2,2n \
         | bedtools merge -i - -c 4 -o count_distinct \
         | awk -v m="$MIN" '$4>=m' | cut -f1-3 > Day${d}_Overlap.bed
done

wc -l Day0_Overlap.bed Day1_Overlap.bed

##Create gene matrix##
bedops -m Day0_Overlap.bed Day1_Overlap.bed > ALL_Peaks.bed

##Compute matrix for plotHeatmap##
computeMatrix scale-regions -S Day0.bw \
                               Day1.bw \
                              -R ALL_Peaks.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros \
                              -o matrix.mat.gz

## PlotHeatmap ##
plotHeatmap -m matrix.mat.gz \
      -out Heatmap_D0_D1.png

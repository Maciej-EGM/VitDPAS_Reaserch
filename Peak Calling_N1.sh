#Peak Calling#


macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R1D0.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R1D0 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R1D2.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R1D1 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R1D2.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R1D2 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R2D0.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R2D0 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R2D1.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R2D1 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R2D2.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R2D2 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R3D0.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R3D0 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R3D1.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R3D1 -B
macs3 callpeak -t HCNFWBGXN_ATAC1_09012023_23s000035-1-1_Carlberg_lane1R3D2.sorted.bam -f BAMPE --verbose=3 --shift 0 -g hs -q 0.01 -n R3D2 -B


wc -l bed
wc -l Peak




#Mergin Bedgraph files in condition order#
bedtools unionbedg -i R1D0_treat_pileup.bdg R2D0_treat_pileup.bdg R3D0_treat_pileup.bdg > Day0.bdg
bedtools unionbedg -i R1D1_treat_pileup.bdg R2D1_treat_pileup.bdg R3D1_treat_pileup.bdg > Day1.bdg
bedtools unionbedg -i R1D2_treat_pileup.bdg R2D2_treat_pileup.bdg R3D2_treat_pileup.bdg > Day2.bdg
##Remove Header##
awk 'NR!=1' Day0.bdg > Day0.bedGraph
awk 'NR!=1' Day1.bdg > Day1.bedGraph
awk 'NR!=1' Day2.bdg > Day2.bedGraph
##Sorting##
sort -k1,1 -k2,2n Day0.bedGraph > Day0.sorted.bedGraph
sort -k1,1 -k2,2n Day1.bedGraph > Day1.sorted.bedGraph
sort -k1,1 -k2,2n Day2.bedGraph > Day2.sorted.bedGraph
##Chromosome Lengh reference##
fetchChromSizes hg38 > hg38.chrom.sizes
##removing 5th column##
awk '{print $1,$2,$3,$4}' Day0.sorted.bedGraph > Day0.bedGraph
awk '{print $1,$2,$3,$4}' Day1.sorted.bedGraph > Day1.bedGraph
awk '{print $1,$2,$3,$4}' Day2.sorted.bedGraph > Day2.bedGraph
##Converting to BigWig##
bedGraphToBigWig Day0.bedGraph hg38.chrom.sizes Day0.bw
bedGraphToBigWig Day1.bedGraph hg38.chrom.sizes Day1.bw
bedGraphToBigWig Day2.bedGraph hg38.chrom.sizes Day2.bw

#Overlap Beetween Peaks/Merging Peaks files in condition order#
bedtools intersect -a R1D0_peaks.narrowPeak -b  R2D0_peaks.narrowPeak R3D0_peaks.narrowPeak -f 0.50 -r > Day0_Overlap.bed
bedtools intersect -a R1D1_peaks.narrowPeak -b  R2D1_peaks.narrowPeak R3D1_peaks.narrowPeak -f 0.50 -r > Day1_Overlap.bed
bedtools intersect -a R1D2_peaks.narrowPeak -b  R2D2_peaks.narrowPeak R3D2_peaks.narrowPeak -f 0.50 -r > Day2_Overlap.bed

##Create gene matrix##
bedops -m Day0_Overlap.bed Day1_Overlap.bed Day2_Overlap.bed

##Compute matrix for plotHeatmap##
computeMatrix scale-regions -S Day0.bw\
                                 Day1.bw \
                                 Day2.bw\
                              -R ALL_Peaks.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                                --skipZeros \
                               -o matrix.gz
## PlotHeatmap ##
plotHeatmap -m matrix.mat.gz \
      -out ExampleHeatmap1.png \

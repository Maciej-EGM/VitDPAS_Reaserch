#!/bin/sh

#SBATCH -J STAR_ATAC
#SBATCH -p batch
#SBATCH -N 4
#SBATCH --ntasks-per-node 24
#SBATCH --mem 64gb
#SBATCH --time 24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=m.jankowski@pan.olsztyn.pl

module load tryton/binutils/2.34
module load tryton/compiler/gcc/7.4
module load tryton/compiler/intel/2019



##Date:20 Jan 2023

## Run using this on command line
## IMPORTANT: DO NOT FORGET TO UPDATE THE SCRIPT NAME BELOW BEFORE RUNNING THIS AFTER MODIFICATIONS!!!
## As an additional safety feature against unintended overwriting of older results, always leave the run
## command below without the concluding '.sh', and only add it manually on the command line
## sh /research/groups/heikkinen/spheikki/scripts/ATACseqPipe_step3_STAR_PE.SE_v1.1_human.THP1s.Carlberg

## Takes as input a folder of quality trimmed fastq read files, and aligns them using STAR.
## After alignment, collects the main output, the sorted .bams, from sample-wise STAR output
## folders into one folder and indexes them. Also runs FastQC on the good sorted.bams and
## finally produces .tdf viz files ("full" as in unfiltered; see step 4 for the final version).

## Note: Make sure all the output folders exist before starting.
## Note: - for pair end analysis, relies on finding the 'R1/2' in file name as the last prefix
##         element before the first period
##          - e.g. 'ATAC_n8_i1_siTCF_R1.ATAC_n8_i1_siTCF_S4_R1_001.fastq.gz'
##       - for single end, the name is simply cut at the first period, so if R1 and R2 exist for
##         a sample, they are aligned separately
## Note: see the bottom of the script for comments on how STAR could be tweaked for RNAseq

## TO DO
## -

## Folder where the subfolders are to be found (do not use concluding slash here or below)
projectFolder="/users/project1/pt01001/VitDPAss_ATAC_seq_Polish_Cohort"


## Common settings
seqfolder=$projectFolder$(echo "/Raw_files/")
##inputSuffix="trimmed.fastq.gz"
inputSuffix=".txt.gz"
STARreadtype="paired" ## possible values 'paired' and 'single'. if neither, script errors out but may still do silly things first
##STARmodule="star/2.7.10b"
STARGenomeIndexPath="/users/scratch1/carlbergc/GRCh38.p14_STAR/"
STARCores="48"
STARreadsubsetSize="-1" ## default is '-1'; anything else makes sense only for testing purposes
STARclip5="0"
STARclip3="0"
samsortCores="16"
STARseedlen="50" ## default is 50
STARendstype="EndToEnd"  ## usually either the default 'Local', or 'EndToEnd' that prevents soft clipping. For ATAC-seq should be 'EndtoEnd' to not soft-clip
SPLICEAwareness="1" ## 1 = off. For ATAC-seq should be set to 1
STARprotrude="0 ConcordantPair" ## default is '0 ConcordantPair'; from the STAR manual (v2.5.4b) it is difficult to understand what this means and how the number would affect things. Therefore, using default for now.
STARbamSortTypes="SortedByCoordinate" ## either Unsorted, SortedByCoordinate, or both  ## 'Unsorted' might be directly useful in some cases, but not really within this pipeline
STARfilterMismatchNmax="10" ## default 10
alnfolder=$projectFolder$(echo "/star_aligned_cohort/")
bamOutputFolder=$projectFolder$(echo "/bam_full_cohort")




### For logging of this pipeline step
logFileBase=$projectFolder$(echo "/log_cohort")
currDateTime=$(date +%y%m%d_%H-%M)           ## don't touch this
logFile="${logFileBase}_${currDateTime}.txt" ## ...or this


{ # this "{" is for logging


echo " "
echo "Loading STAR module"
echo " "


echo " "
echo "Starting alignment"
date
echo " "

echo " "
echo "COMMON SETTINGS"
echo "  Sequence (input) folder            : "${seqfolder}
echo "  Main output folder                 : "${alnfolder}
echo "  STAR genome index base (full path) : "${STARGenomeIndexPath}
echo "  Read type to align                 : "${STARreadtype}
echo "  Cores for alignment                : "${STARCores}
echo "  Cores for filtering and sorting    : "${samsortCores}
echo "  Read subset size (-1 = all reads)  : "${STARreadsubsetSize}
echo "  Clip nts from 5' of each read      : "${STARclip5}
echo "  Clip nts from 3' of each read      : "${STARclip3}
echo "  Seed length                        : "${STARseedlen}
echo "  Align ends type                    : "${STARendstype}
echo "  Align protrude                     : "${STARprotrude}
echo "  BAM output sort type(s)            : "${STARbamSortTypes}
echo "  Log file                           : "${logFile}
echo " "
date
echo " "

cd "${seqfolder}"

SAMPLES="AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane110d0
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane110d1
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane111d0
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane111d1
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane113d0
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane114d0
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane115d1
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane119d1
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane121d0
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane126d1
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane15d0
AAFKGVTM5_ATACP4_260424_24s000684-1-1_Carlberg_lane19d1
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane115d0
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane116d0
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane117d0
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane120d1
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane121d1
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane123d0
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane123d1
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane125d0
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane153d1
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane15d1
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane16d0
AAFKH27M5_ATACP3_260424_24s000683-1-1_Carlberg_lane19d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane118d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane118d1
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane128d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane129d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane130d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane134d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane135d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane135d1
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane143d1
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane153d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane154d0
AAFKH37M5_ATACP2_260424_24s000682-1-1_Carlberg_lane154d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane130d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane133d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane137d0
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane138d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane139d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane140d0
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane144d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane145d0
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane149d0
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane149d1
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane150d0
AAFKH3HM5_ATACP7_260424_24s000687-1-1_Carlberg_lane152d0
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane112d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane114d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane11d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane124d0
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane125d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane126d0
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane129d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane13d0
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane13d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane14d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane16d1
AAFKJJMM5_ATACP6_260424_24s000686-1-1_Carlberg_lane18d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane116d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane117d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane133d0
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane136d0
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane137d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane138d0
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane139d0
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane140d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane146d0
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane150d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane152d1
AAFKKG5M5_ATACP1_260424_24s000681-1-1_Carlberg_lane18d0
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane112d0
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane113d1
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane119d0
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane11d0
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane120d0
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane122d0
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane122d1
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane124d1
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane128d1
AAFM7L7M5_ATACP5_260424_24s000685-1-1_Carlberg_lane14d0
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane134d1
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane136d1
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane143d0
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane144d0
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane145d1
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane146d1
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane147d0
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane147d1
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane148d0
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane148d1
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane151d0
AAFM7T5M5_ATACP8_260424_24s000688-1-1_Carlberg_lane151d1"


echo " "
echo "Processing these sample ids"
echo " "


for SAMPLE in $SAMPLES; do

    echo $SAMPLE
    echo "Processing sample         : ${SAMPLE}"
    echo "full name"
    echo ${seqfolder}${SAMPLE}_1_txt.gz


    if [ "${STARreadtype}" = "paired" ]; then

        infileR1=${seqfolder}${SAMPLE}_1_sequence.txt.gz
        infileR2=${seqfolder}${SAMPLE}_2_sequence.txt.gz


        echo "  Paired input for R1     : ${infileR1}"
        echo "  Paired input for R2     : ${infileR2}"
    fi

    if [ "${STARreadtype}" = "single" ]; then
         infiles=$(ls | grep "$sampleid" )
        echo "  Single-end input string : ${infiles}"
    fi


    outfolderSample=$(echo $alnfolder$SAMPLE$(echo "/"))
    mkdir ${outfolderSample}

    echo "  Output location         : ${outfolderSample}"
    echo " "

	STAR --runThreadN "${STARCores}" \
	--runMode "alignReads" \
	--genomeDir "${STARGenomeIndexPath}" \
	--genomeLoad "LoadAndKeep" \
	--readFilesIn $infileR1 $infileR2 \
	--readMapNumber "${STARreadsubsetSize}" \
	--clip3pNbases "${STARclip3}" \
	--clip5pNbases "${STARclip5}" \
	--seedSearchStartLmax "${STARseedlen}" \
	--alignEndsType "${STARendstype}" \
 	--alignEndsProtrude "${STARprotrude}" \
 	--outFilterMismatchNmax "${STARfilterMismatchNmax}" \
    --alignIntronMax "${SPLICEAwareness}" \
	--outFileNamePrefix "${outfolderSample}" \
	--outSAMtype BAM "${STARbamSortTypes}" \
	--limitBAMsortRAM "60000000000" \
	--outSAMattributes "All" \
	--outBAMsortingThreadN "${samsortCores}" \
	--outStd "Log" \
	--readFilesCommand zcat


###	--readFilesIn "${infileR1}" "${infileR2}" \


  ## Moving the output bam file to bam output folder
  mv -v "${outfolderSample}Aligned.sortedByCoord.out.bam" "${bamOutputFolder}/${SAMPLE}.sorted.bam"

done

STAR --genomeDir "${STARGenomeIndexPath}" --genomeLoad Remove

echo " "
echo "DONE all alignments!"
echo " "
date
echo " "
echo " "


echo " "
echo "Starting bam indexing"
date
echo " "

cd "${bamOutputFolder}"

for bfile in $( ls | grep "\.sorted\.bam$" ); do

  echo " "
  echo "Indexing file: ${bfile}"
  bamsortinputFile="${bamOutputFolder}/${bfile}"

  samtools index "${bamsortinputFile}"
  samtools idxstats "${bamsortinputFile}" > "${bamsortinputFile}.idxstat.txt"
  samtools stats "${bamsortinputFile}" > "${bamsortinputFile}.stats.txt"

  echo "Indexing of ${bfile} complete."
  echo " "

done

echo " "
echo "Bam indexing complete"
date
echo " "





} 2>&1 | tee $logFile # This makes all the outputs that are printed
#+ to terminal between the first { and this one also to be printed in the specified log file


######## BELOW HERE, suggestions for STAR parameter tweaking for RNA-seq from Baruzzo et al, Nat Meth 2017
# "Conclusions
# The default options of STAR achieve one of the best results on the T3 complexity dataset. However, tweaking some parameters is possible to further improve the metrics. The important roles of NUM_FILTER_MISMATCHES (outFilterMismatchNmax), END_ALIGNMENT_TYPE (alignEndsType), OVERHANG (alignSJoverhangMin), NUM_FILTER_SCORE (outFilterScoreMinOverLread), SEED_LENGTH (seedSearchStartLmax) are confirmed on both Human and Malaria. Increasing the number of allowed mismatches and leaving the default value for END_ALIGNMENT_TYPE improve the results. At the same time, increasing OVERHANG and decreasing SEED_LENGTH and NUM_FILTER_SCORE increase the recall at all levels."
# 	- the tested parameters (in the order given in the optimal string)
# 			1. limitOutSJcollapsed (default optimal)
# 			2. limitSjdbInsertNsj (default optimal)
# 		3. outFilterMultimapNmax (100; default 10)
# 		4. outFilterMismatchNmax (33; default 10)
# 			5. outFilterMismatchNoverLmax (default is optimal)
# 		6. seedSearchStartLmax (12; default 50)
# 		7. alignSJoverhangMin (15; default 5)
# 			8. alignEndsType (default optimal)
# 		9. outFilterMatchNminOverLread (0; default 0.66)
# 		10. outFilterScoreMinOverLread (0.3; default 0.66)
# 			11. winAnchorMultimapNmax (default optimal)
# 			12. alignSJDBoverhangMin (default optimal, except for junctions (1))
# 		13. outFilterType (BySJout; default Normal)
# 	- best read level recall     : 1000000-1000000-100-33-0.3-12-15-Local-0-0.3-50-3-BySJout
# 	- best base level recall     : 1000000-1000000-100-33-0.3-12-15-Local-0-0.3-50-3-BySJout
# 	- best junction level recall : 1000000-1000000-100-33-0.3-12-15-Local-0-0.3-50-1-BySJout ## note the difference in the second-to-last parameter
# ==> it is entirely possible that for non-RNA-seq applications the above optima may not apply
# ==> the STAR defaults where shown to be pretty good even for RNAseq, so let's not introduce any particular tweaks to the ATACseq pipe at this point
#
#



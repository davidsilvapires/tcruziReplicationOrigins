#!/usr/bin/env bash

###################################################################################################
#                                        deeptoolsMain                                            #  
#                                                                                                 #
#  This script allows verifying the co-location of the orc1cdc6 peaks in the genes of the         #
# multigenic families, in the features of the genome and also in the origins of replication       #
# identified in this  work. For this analysis we will use the deeptools compute matrix and        #
# plotheatmap tools. Script  also writes another script to run on Slurm.                          #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                   < thiago.franco.esib@esib.butantan.gov.br>                                    #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/deeptoolsMain.slurm <<EOI
#!/usr/bin/env bash

# Creating genome index.
BASENAME_GENOME=`basename \${GENOME%.fa}`
samtools faidx \${GENOME} \\
  -o output/\${BASENAME_GENOME}.fa.fai \\
  2> log/index.err > log/index.out

# Creating chrom.sizes.
# https://genomewiki.ucsc.edu/index.php/GBiB:_From_download_to_BLAT_at_assemply_hubs
faToTwoBit \${GENOME} output/\${BASENAME_GENOME}.2bit
twoBitInfo output/\${BASENAME_GENOME}.2bit stdout | sort -k2nr \\
  > output/\${BASENAME_GENOME}-chromSizes.txt

# Sorting the bed file that represents the region that will be plotted
BASENAME_BEDFILE=`basename \${BEDFILE%.bed}`
sort -k1,1 -k2,2n -k3,3n \${BEDFILE} \\
  > output/\${BASENAME_BEDFILE}-sorted.bed

# Sorting all bed files used in analysis.
for FILE in \`ls input/*Origins.bed | xargs -i basename {}\`
do
  sort -k1,1 -k2,2n -k3,3n input/\${FILE} > \\
output/\${FILE%.bed}-sorted.bed
done

# Converting bed to bedgraph.
for BEDFILE in \`ls output/*Origins-sorted.bed | xargs -i basename {}\`
do
  bedtools genomecov -bg -i output/\${BEDFILE} \\
    -g output/\${BASENAME_GENOME}.fa.fai \\
    > output/\${BEDFILE%-sorted.bed}.bedgraph \\
  2> log/\${BEDFILE%-sorted.bed}.bedgraph.err
done

# Converting bedgraph to bigwig.
for BEDGRAPH_FILE in \`ls output/*.bedgraph | xargs -i basename {}\`
do
  bedGraphToBigWig output/\${BEDGRAPH_FILE} \\
    output/\${BASENAME_GENOME}-chromSizes.txt \\
    output/\${BEDGRAPH_FILE%.bedgraph}.bigwig \\
    2> log/\${BEDGRAPH_FILE%.bedgraph}.bigwig.err
done

# Compute Matrix to predominant, flexible, dormant
for FILE in \`ls output/*bigwig | xargs -i basename {}\`
do
  computeMatrix scale-regions -S output/\${FILE} \\
  -R output/\${BASENAME_BEDFILE}-sorted.bed \\
  -b 3000 -a 3000 \\
  --verbose \\
  -o output/matrix-\${BASENAME_BEDFILE}-\${FILE%.bigwig}.gz \\
  2> log/matrix-\${BASENAME_BEDFILE}-\${FILE%.bigwig}.err \\
  > log/matrix-\${BASENAME_BEDFILE}-\${FILE%.bigwig}.out
done

# Plot Heatmap to dormant, flexible and predominant origins
for FILEANDCOLOR in dormant,blue flexible,orange predominant,forestgreen
do
   Y=\$(echo \${FILEANDCOLOR} | cut -d, -f1)
   C=\$(echo \${FILEANDCOLOR} | cut -d, -f2)
   plotHeatmap -m  output/matrix-\${BASENAME_BEDFILE}-\${Y}Origins.gz \\
      --yMin 0 --yMax 4 \\
      --averageTypeSummaryPlot "mean" \\
      --dpi 700 \\
      --colorList "white,\${C}" \\
      --startLabel 'S' --endLabel 'E' \\
      --zMin 0 --zMax 4 --hclust 3 \\
      --heatmapHeight 56 --heatmapWidth 8 \\
      --plotFileFormat pdf \\
      --samplesLabel "Orc1cd6 \${Y} Origins : \${BASENAME_BEDFILE}" \\
      -x "\${BASENAME_BEDFILE}" \\
      -y "Orc1Cdc6 peak count" \\
      --outFileSortedRegions output/sortedRegions-\${BASENAME_BEDFILE}-\${Y}OriginsCluster3.bed \\
      -o output/heatmap-\${BASENAME_BEDFILE}-\${Y}OriginsCluster3.pdf \\
      2> log/heatmap-\${BASENAME_BEDFILE}-\${Y}OriginsCluster3.err \\
      > log/heatmap-\${BASENAME_BEDFILE}-\${Y}OriginsCluster3.out
done

exit 0
EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/sbatch-deeptoolsMain.time memusg \
  --output-to-file log/sbatch-deeptoolsMain.memusg --time --shell "saveCommand sbatch --nodes 1 \
  --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} -o log/slurm-%A.out -J deeptoolsMain job/deeptoolsMain.slurm \
  2> log/sbatch-deeptoolsMain.err | tee log/sbatch-deeptoolsMain.out"

# Check how are going all the running awk processes. 
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0




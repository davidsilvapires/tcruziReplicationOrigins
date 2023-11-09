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
samtools faidx \${GENOME} \
  -o output/\${BASENAME_GENOME}.fa.fai \
  2> log/index.err > log/index.out

# Creating chrom.sizes.
# https://genomewiki.ucsc.edu/index.php/GBiB:_From_download_to_BLAT_at_assemply_hubs
faToTwoBit \${GENOME} output/\${BASENAME_GENOME}.2bit
twoBitInfo output/\${BASENAME_GENOME}.2bit stdout | sort -k2nr \
  > output/\${BASENAME_GENOME}-chromSizes.txt

# Sorting the bed file that represents the region that will be plotted
BASENAME_BEDFILE=`basename \${BEDFILE%.bed}`
sort -k1,1 -k2,2n -k3,3n \${BEDFILE} \
  > output/\${BASENAME_BEDFILE}-sorted.bed

# Coverting bed to bedgraph
bedtools genomecov -bg -i output/\${BASENAME_BEDFILE}-sorted.bed \
  -g output/\${BASENAME_GENOME}.fa.fai \
  > output/\${BASENAME_BEDFILE}.bedgraph \
  2> log/\${BASENAME_BEDFILE}.bedgraph.err

# Converting bedgraph to bigwig.
bedGraphToBigWig output/\${BASENAME_BEDFILE}.bedgraph \
  output/\${BASENAME_GENOME}-chromSizes.txt \
  output/\${BASENAME_BEDFILE}.bigwig \
  2> log/\${BASENAME_BEDFILE}.bigwig.err

# Sorting all bed files used in analysis.
for FILE in \`ls input/*.bed | xargs -i basename {}\`
do
  sort -k1,1 -k2,2n -k3,3n input/\${FILE} > \
    output/\${FILE%.bed}-sorted.bed
done

# Compute Matrix to predominant and flexible Origins
for FILE in \`ls output/*Origins-sorted.bed | xargs -i basename {} \`
do
  computeMatrix scale-regions -S output/\${BASENAME_BEDFILE}.bigwig \
    -R output/\${FILE} \
    -b 3000 -a 3000 \
    --verbose \
    -o output/matrix-\${FILE%-sorted.bed}-\${BASENAME_BEDFILE}.gz \
    2> log/matrix-\${FILE%-sorted.bed}-\${BASENAME_BEDFILE}.err \
    > log/matrix-\${FILE%-sorted.bed}-\${BASENAME_BEDFILE}.out 
done

# Plot Heatmap to flexible and predominant origins
for FILEANDCOLOR in flexible,orange predominant,forestgreen 
do
  X=\$(echo \${FILEANDCOLOR} | cut -d, -f1)
  C=\$(echo \${FILEANDCOLOR} | cut -d, -f2)
  plotHeatmap -m output/matrix-\${X}Origins-\${BASENAME_BEDFILE}.gz \
    --yMin 0 --yMax 1 \
    --averageTypeSummaryPlot "mean" \
    --dpi 700 \
    --colorList "white,\${C}" \
    --startLabel 'S' --endLabel 'E' \
    --zMin 0 --zMax 2 --hclust 3 \
    --heatmapHeight 56 --heatmapWidth 8 \
    --plotFileFormat pdf \
    --samplesLabel "\${BASENAME_BEDFILE}:\${X} Origins" \
    -x "\${X} Origins" \
    -y "\${BASENAME_BEDFILE} peak count" \
    --outFileSortedRegions output/sortedRegions-\${X}Origins-\${BASENAME_BEDFILE}-Cluster3.bed \
    -o output/heatmap-\${X}Origins-\${BASENAME_BEDFILE}-3kCluster3.pdf \
    2> log/heatmap-\${X}Origins-\${BASENAME_BEDFILE}-3kCluster3.err \
    > log/heatmap-\${X}Origins-\${BASENAME_BEDFILE}-3kCluster3.out

done

# Compute Matrix to genome feature 
for FILE in \`ls output/*Feature-sorted.bed | xargs -i basename {}\`
do
  computeMatrix scale-regions -S output/\${BASENAME_BEDFILE}.bigwig \
    -R output/\${FILE} \
    -b 3000 -a 3000 \
    --verbose \
    -o output/matrix-\${FILE%-Feature-sorted.bed}-\${BASENAME_BEDFILE}.gz \
    2>log/matrix-\${FILE%-Feature-sorted.bed}-\${BASENAME_BEDFILE}.err \
    >log/matrix-\${FILE%-Feature-sorted.bed}-\${BASENAME_BEDFILE}.out 

done

# Plot Heatmap to genomic feature
for FILEANDCOLOR in  snoRNA,green intergenicPolycistronRegions,orange
do
  X=\$(echo \${FILEANDCOLOR} | cut -d, -f1)
  C=\$(echo \${FILEANDCOLOR} | cut -d, -f2)
  plotHeatmap -m output/matrix-\${X}-\${BASENAME_BEDFILE}.gz \
    --yMin 0 --yMax 1 \
    --averageTypeSummaryPlot "mean" \
    --dpi 700 \
    --colorList "white,\${C}" \
    --startLabel 'S' --endLabel 'E' \
    --zMin 0 --zMax 2 --hclust 3 \
    --heatmapHeight 56 --heatmapWidth 8 \
    --plotFileFormat pdf \
    --samplesLabel "\${BASENAME_BEDFILE} : \${X}" \
    -x "\${X}" \
    -y "\${BASENAME_BEDFILE} peak count" \
    --outFileSortedRegions output/sortedRegions-\${X}-\${BASENAME_BEDFILE}-Cluster3.bed \
    -o output/heatmap-\${X}-\${BASENAME_BEDFILE}-3kCluster3.pdf \
    2> log/heatmap-\${X}-\${BASENAME_BEDFILE}-3kCluster3.err \
    > log/heatmap-\${X}-\${BASENAME_BEDFILE}-3kCluster3.out \

done

# Compute Matrix to multigenic families 
for FILE in \`ls output/*Gene-sorted.bed | xargs -i basename {}\`
do
  computeMatrix scale-regions -S output/\${BASENAME_BEDFILE}.bigwig \
    -R output/\${FILE} \
    -b 3000 -a 3000 \
    --verbose \
    -o output/matrix-\${FILE%-Gene-sorted.bed}-\${BASENAME_BEDFILE}.gz \
    2>log/matrix-\${FILE%-Gene-sorted.bed}-\${BASENAME_BEDFILE}.err \
    >log/matrix-\${FILE%-Gene-sorted.bed}-\${BASENAME_BEDFILE}.out 
done

# Plot Heatmap to dormant, flexible and predominant origins
for FILEANDCOLOR in  DGF-1,green MASP,red RHS,blue TS,orange mucin,forestgreen
do
  X=\$(echo \${FILEANDCOLOR} | cut -d, -f1)
  C=\$(echo \${FILEANDCOLOR} | cut -d, -f2)
  plotHeatmap -m output/matrix-\${X}-\${BASENAME_BEDFILE}.gz \
    --yMin 0 --yMax 1 \
    --averageTypeSummaryPlot "mean" \
    --dpi 700 \
    --colorList "white,\${C}" \
    --startLabel 'S' --endLabel 'E' \
    --zMin 0 --zMax 2 --hclust 3 \
    --heatmapHeight 56 --heatmapWidth 8 \
    --plotFileFormat pdf \
    --samplesLabel "\${BASENAME_BEDFILE} : \${X}" \
    -x "\${X}" \
    -y "\${BASENAME_BEDFILE} peak count" \
    --outFileSortedRegions output/sortedRegions-\${X}-\${BASENAME_BEDFILE}-Cluster3.bed \
    -o output/heatmap-\${X}-\${BASENAME_BEDFILE}-3kCluster3.pdf \
    2> log/heatmap-\${X}-\${BASENAME_BEDFILE}-3kCluster3.err \
    > log/heatmap-\${X}-\${BASENAME_BEDFILE}-3kCluster3.out

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

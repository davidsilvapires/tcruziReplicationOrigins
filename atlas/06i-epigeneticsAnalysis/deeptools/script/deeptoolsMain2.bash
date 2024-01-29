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

# Sorting the bed from Orc1Cdc6 
 sort -k1,1 -k2,2n -k3,3n input/predominantOrigins.bed  > output/predominant-sorted.bed

# Compute Matrix to Variant histone H2BV  co-localization with Orc1Cdc6
computeMatrix scale-regions -S input/Epi2-h2bv.bw \\
  -R output/predominant-sorted.bed \\
  -b 10000 -a 10000 \\
  --verbose \\
  --skipZeros \\
  -o output/matrix-HistoneVariant-h2bv-Predominant.gz \\
  2> log/matrix-HistoneVariant-h2bv-Predominant.err\\
  > log/matrix-HistoneVariant-h2bv-Predominant.out


# Plot Heatmap to h2bv  histones variants
  plotHeatmap -m output/matrix-HistoneVariant-h2bv-Predominant.gz  \\
    --yMin 0 \\
    --yMax 10 \\
    --averageTypeSummaryPlot "mean" \\
    --dpi 700 \\
    --colorList 'white,forestgreen' \\
    --zMin 0 \\
    --zMax 10 \\
    --kmeans 1 \\
    --heatmapHeight 56 \\
    --heatmapWidth 8 \\
    --samplesLabel "Epi-H2BV" \\
    --startLabel 'S' \\
    --endLabel 'E' \\
    --plotFileFormat pdf \\
    -x "Predominant Origins" \\
    -y "Histone H2BV peak count" \\
    --outFileSortedRegions output/sortedRegions-H2BV-Predominant-kmeans-10kb.bed \\
    -o output/heatmap-H2BV-Predominant-kmeans-10kb.pdf \\
    2> log/heatmap-H2BV-Predominant-Kmeans-10kb.err \\
    > log/heatmap-H2BV-Predominant-kmeans-10kb.out

exit 0
EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-deeptoolsMain.time memusg \
  --output-to-file log/sbatch-deeptoolsMain.memusg \
  --time \
  --shell "saveCommand sbatch --nodes 1 \
  --ntasks ${NUM_OF_CPUS} \
  --mem ${MEMORY_SIZE} \
  -o log/slurm-%A.out \
  -J deeptoolsMain job/deeptoolsMain.slurm \
  2> log/sbatch-deeptoolsMain.err | 
  tee log/sbatch-deeptoolsMain.out"

# Check how are going all the running awk processes. 
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0




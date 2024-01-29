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
cat> job/deeptoolsMain7.slurm <<EOI
#!/usr/bin/env bash

# Compute Matrix to Variant histone H2BV  co-localization with Orc1Cdc6
computeMatrix scale-regions -S input/Epi2-h2bv.bw \\
  -R output/radon_sequencesFlexible.bed \\
  -b 10000 -a 10000 \\
  --verbose \\
  -o output/matrix-HistoneVariant-h2bv-Flexible-Control.gz \\
  2> log/matrix-HistoneVariant-h2bv-Flexible-Control.err\\
  > log/matrix-HistoneVariant-h2bv-Flexible-Control.out


# Plot Heatmap to h2bv  histones variants
  plotHeatmap -m output/matrix-HistoneVariant-h2bv-Flexible-Control.gz  \\
    --yMin 0 \\
    --yMax 10 \\
    --averageTypeSummaryPlot "mean" \\
    --dpi 700 \\
    --colorList 'white,black' \\
    --zMin 0 --zMax 10 \\
    --kmeans 1 \\
    --heatmapHeight 56  \\
    --heatmapWidth 8 \\
    --samplesLabel "Epi-H2BV" \\
    --startLabel 'S' \\
    --endLabel 'E' \\
    --plotFileFormat pdf \\
    -x "Control Sequences" \\
    -y "Histone H2BV peak count" \\
    --outFileSortedRegions output/sortedRegions-H2BV-Flexible-Kmeans-10kb-Control.bed \\
    -o output/heatmap-H2BV-Flexible-Kmeans-10kb-Control.pdf \\
    2> log/heatmap-H2BV-Flexible-Kmeans-10kb-Control.err \\
    > log/heatmap-H2BV-Flexible-Kmeans-10kb-Control.out

exit 0
EOI

 MAIN_JOB_ID=`cat log/sbatch-deeptoolsMain5.out`
# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose --output=log/sbatch-deeptoolsMain.time memusg \
  --output-to-file log/sbatch-deeptoolsMain.memusg \
  --time \
  --shell "saveCommand sbatch --nodes 1 \
  --dependency=afterany:${MAIN_JOB_ID} \
  --ntasks ${NUM_OF_CPUS} \
  --mem ${MEMORY_SIZE} \
  -o log/slurm-%A.out \
  -J deeptoolsMain job/deeptoolsMain7.slurm \
  2> log/sbatch-deeptoolsMain7.err | 
  tee log/sbatch-deeptoolsMain7.out"

# Check how are going all the running awk processes. 
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0




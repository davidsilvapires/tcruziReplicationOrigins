#mkdir -p output/mergePeaks/histone/unique/{control/individualSample,taggedProtein/individualSample}

# UNIQUE
mergePeaks -d 340 -prefix merged -matrix statsUniqueOrc1Cdc6 -gsize 25827529 -venn vennDiagramUniqueOrc1Cdc6.txt fpUniqueOrc1Cdc6Rep1.txt fpUniqueOrc1Cdc6Rep2.txt fpUniqueOrc1Cdc6Rep3.txt 2> log/mergePeaksHistoneUniqueTaggedProtein1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsUniqueControl -gsize 25827529 -venn vennDiagramUniqueControl.txt fpUniqueControlRep1.txt fpUniqueControlRep2.txt fpUniqueControlRep3.txt 2> log/mergePeaksHistoneUniqueControl1e2e3.erR
mergePeaks -d 340 -prefix merged -matrix statsUniqueOrc1Cdc6Control -gsize 25827529 -venn vennDiagramUniqueOrc1Cdc6Control.txt merged_fpUniqueOrc1Cdc6Rep1.txt_fpUniqueOrc1Cdc6Rep2.txt_fpUniqueOrc1Cdc6Rep3.txt merged_fpUniqueControlRep1.txt_fpUniqueControlRep2.txt_fpUniqueControlRep3.txt  2> log/mergePeaksHistoneUniqueOrc1Cdc6Control.err
pos2bed.pl merged_fpUniqueControlRep1.txt_fpUniqueControlRep2.txt_fpUniqueControlRep3.txt > mergePeaksUniqueControl.bed
pos2bed.pl merged_fpUniqueOrc1Cdc6Rep1.txt_fpUniqueOrc1Cdc6Rep2.txt_fpUniqueOrc1Cdc6Rep3.txt > mergePeaksUniqueOrc1Cdc6.bed
pos2bed.pl merged_merged_fpUniqueOrc1Cdc6Rep1.txt_fpUniqueOrc1Cdc6Rep2.txt_fpUniqueOrc1Cdc6Rep3.txt_merged_fpUniqueControlRep1.txt_fpUniqueControlRep2.txt_fpUniqueControlRep3.txt > mergePeaksUniqueOrc1Cdc6Control.bed

# MULTIMAPPERS
mergePeaks -d 340 -prefix merged -matrix statsMultimappersOrc1Cdc6 -gsize 25827529 -venn vennDiagramMultimappersOrc1Cdc6.txt fpMultimappersOrc1Cdc6Rep1.txt fpMultimappersOrc1Cdc6Rep2.txt fpMultimappersOrc1Cdc6Rep3.txt 2> log/mergePeaksHistoneMultimappersTaggedProtein1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsMultimappersControl -gsize 25827529 -venn vennDiagramMultimappersControl.txt fpMultimappersControlRep1.txt fpMultimappersControlRep2.txt fpMultimappersControlRep3.txt 2> log/mergePeaksHistoneMultimappersControl1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsMultimappersOrc1Cdc6Control -gsize 25827529 -venn vennDiagramMultimappersOrc1Cdc6Control.txt merged_fpMultimappersOrc1Cdc6Rep1.txt_fpMultimappersOrc1Cdc6Rep2.txt_fpMultimappersOrc1Cdc6Rep3.txt merged_fpMultimappersControlRep1.txt_fpMultimappersControlRep2.txt_fpMultimappersControlRep3.txt 2> log/mergePeaksHistoneMultimappersOrc1Cdc6Control.err
pos2bed.pl merged_fpMultimappersControlRep1.txt_fpMultimappersControlRep2.txt_fpMultimappersControlRep3.txt > mergePeaksMultimappersControl.bed
pos2bed.pl merged_fpMultimappersOrc1Cdc6Rep1.txt_fpMultimappersOrc1Cdc6Rep2.txt_fpMultimappersOrc1Cdc6Rep3.txt > mergePeaksMultimappersOrc1Cdc6.bed
pos2bed.pl merged_merged_fpMultimappersOrc1Cdc6Rep1.txt_fpMultimappersOrc1Cdc6Rep2.txt_fpMultimappersOrc1Cdc6Rep3.txt_merged_fpMultimappersControlRep1.txt_fpMultimappersControlRep2.txt_fpMultimappersControlRep3.txt > mergePeaksMultimappersOrc1Cdc6Control.bed


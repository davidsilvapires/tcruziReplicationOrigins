#mkdir -p output/mergePeaks/histone/unique/{control/individualSample,taggedProtein/individualSample}

# UNIQUE
mergePeaks -d 340 -prefix merged -matrix statsFactorUniqueOrc1Cdc6 -gsize 25827529 -venn vennDiagramFactorUniqueOrc1Cdc6.txt fpFactorUniqueOrc1Cdc6Rep1.txt fpFactorUniqueOrc1Cdc6Rep2.txt fpFactorUniqueOrc1Cdc6Rep3.txt 2> log/mergePeaksFactorUniqueTaggedProtein1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsFactorUniqueControl -gsize 25827529 -venn vennDiagramFactorUniqueControl.txt fpFactorUniqueControlRep1.txt fpFactorUniqueControlRep2.txt fpFactorUniqueControlRep3.txt 2> log/mergePeaksFactorUniqueControl1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsFactorUniqueOrc1Cdc6Control -gsize 25827529 -venn vennDiagramFactorUniqueOrc1Cdc6Control.txt merged_fpFactorUniqueOrc1Cdc6Rep1.txt_fpFactorUniqueOrc1Cdc6Rep2.txt_fpFactorUniqueOrc1Cdc6Rep3.txt merged_fpFactorUniqueControlRep1.txt_fpFactorUniqueControlRep2.txt_fpFactorUniqueControlRep3.txt  2> log/mergePeaksFactorUniqueOrc1Cdc6Control.err
pos2bed.pl merged_fpFactorUniqueControlRep1.txt_fpFactorUniqueControlRep2.txt_fpFactorUniqueControlRep3.txt > mergePeaksFactorUniqueControl.bed
pos2bed.pl merged_fpFactorUniqueOrc1Cdc6Rep1.txt_fpFactorUniqueOrc1Cdc6Rep2.txt_fpFactorUniqueOrc1Cdc6Rep3.txt > mergePeaksFactorUniqueOrc1Cdc6.bed
pos2bed.pl merged_merged_fpFactorUniqueOrc1Cdc6Rep1.txt_fpFactorUniqueOrc1Cdc6Rep2.txt_fpFactorUniqueOrc1Cdc6Rep3.txt_merged_fpFactorUniqueControlRep1.txt_fpFactorUniqueControlRep2.txt_fpFactorUniqueControlRep3.txt > mergePeaksFactorUniqueOrc1Cdc6Control.bed

# MULTIMAPPERS
mergePeaks -d 340 -prefix merged -matrix statsFactorMultimappersOrc1Cdc6 -gsize 25827529 -venn vennDiagramFactorMultimappersOrc1Cdc6.txt fpFactorMultimappersOrc1Cdc6Rep1.txt fpFactorMultimappersOrc1Cdc6Rep2.txt fpFactorMultimappersOrc1Cdc6Rep3.txt 2> log/mergePeaksFactorMultimappersTaggedProtein1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsFactorMultimappersControl -gsize 25827529 -venn vennDiagramFactorMultimappersControl.txt fpFactorMultimappersControlRep1.txt fpFactorMultimappersControlRep2.txt fpFactorMultimappersControlRep3.txt 2> log/mergePeaksFactorMultimappersControl1e2e3.err
mergePeaks -d 340 -prefix merged -matrix statsFactorMultimappersOrc1Cdc6Control -gsize 25827529 -venn vennDiagramFactorMultimappersOrc1Cdc6Control.txt merged_fpFactorMultimappersOrc1Cdc6Rep1.txt_fpFactorMultimappersOrc1Cdc6Rep2.txt_fpFactorMultimappersOrc1Cdc6Rep3.txt merged_fpFactorMultimappersControlRep1.txt_fpFactorMultimappersControlRep2.txt_fpFactorMultimappersControlRep3.txt 2> log/mergePeaksFactorMultimappersOrc1Cdc6Control.err
pos2bed.pl merged_fpFactorMultimappersControlRep1.txt_fpFactorMultimappersControlRep2.txt_fpFactorMultimappersControlRep3.txt > mergePeaksFactorMultimappersControl.bed
pos2bed.pl merged_fpFactorMultimappersOrc1Cdc6Rep1.txt_fpFactorMultimappersOrc1Cdc6Rep2.txt_fpFactorMultimappersOrc1Cdc6Rep3.txt > mergePeaksFactorMultimappersOrc1Cdc6.bed
pos2bed.pl merged_merged_fpFactorMultimappersOrc1Cdc6Rep1.txt_fpFactorMultimappersOrc1Cdc6Rep2.txt_fpFactorMultimappersOrc1Cdc6Rep3.txt_merged_fpFactorMultimappersControlRep1.txt_fpFactorMultimappersControlRep2.txt_fpFactorMultimappersControlRep3.txt > mergePeaksFactorMultimappersOrc1Cdc6Control.bed


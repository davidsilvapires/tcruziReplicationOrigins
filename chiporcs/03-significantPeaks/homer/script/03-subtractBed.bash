mkdir -p output/finalPeaks

bedtools subtract -A -a mergePeaksUniqueOrc1Cdc6.bed -b mergePeaksUniqueControl.bed > output/finalPeaks/finalPeaksUniqueOrc1Cdc6allReplicates.bed 2> log/subtractUniqueOrc1Cdc6Control.err
bedtools subtract -A -a mergePeaksMultimappersOrc1Cdc6.bed -b mergePeaksMultimappersControl.bed > output/finalPeaks/finalPeaksMultimappersOrc1Cdc6allReplicates.bed 2> log/subtractMultimappersOrc1Cdc6Control.err

# Generating fasta files for MEME motif finder.

mkdir -p output/findMotifs/getfasta/

bedtools getfasta -fi input/tryCru-clb1.fasta -bed output/finalPeaks/finalPeaksUniqueOrc1Cdc6allReplicates.bed > output/findMotifs/getfasta/finalPeaksUniqueOrc1Cdc6allReplicates.fasta 2> log/getfastaUniqueOrc1cdc6.err
bedtools getfasta -fi input/tryCru-clb1.fasta -bed output/finalPeaks/finalPeaksMultimappersOrc1Cdc6allReplicates.bed > output/findMotifs/getfasta/finalPeaksMultimappersOrc1Cdc6allReplicates.fasta 2> log/getfastaMultimappersOrc1cdc6.err

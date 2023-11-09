# Finding peaks with --style histone
mkdir -p output/findPeaks/{unique,multimappers}

# findPeaks <tag directory> -style histone -size 500 -minDist 1000  -i <control tag directory> -o <path to output file>
# -gsize = genome size. According to the tritrypDB website (Data -> Organisms: Genome Info and Stats -> search for Trypanosoma cruzi CL Brener.
#          Esmeraldo-like = 32.53Mbp = 32530000
#          CL Brener = 36.03Mbp = 36030000

# UNIQUE

findPeaks output/tagDirectory/unique/control/individualSample/chip/2/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/control/individualSample/input/1/ -o fpFactorUniqueControlRep1.txt 2> log/findPeaksFactorUniqueControlIndividualSampleRepĺicate1.err > log/findPeaksFactorUniqueControlIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/unique/control/individualSample/chip/4/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/control/individualSample/input/3/ -o fpFactorUniqueControlRep2.txt 2> log/findPeaksFactorUniqueControlIndividualSampleRepĺicate2.err > log/findPeaksFactorUniqueControlIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/unique/control/individualSample/chip/6/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/control/individualSample/input/5/ -o fpFactorUniqueControlRep3.txt 2> log/findPeaksFactorUniqueControlIndividualSampleRepĺicate3.err > log/findPeaksFactorUniqueControlIndividualSampleRepĺicate3.out
findPeaks output/tagDirectory/unique/taggedProtein/individualSample/chip/8/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/taggedProtein/individualSample/input/7/ -o fpFactorUniqueOrc1Cdc6Rep1.txt 2> log/findPeaksFactorUniqueTaggedProteinIndividualSampleRepĺicate1.err > log/findPeaksFactorUniqueTaggedProteinIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/unique/taggedProtein/individualSample/chip/10/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/taggedProtein/individualSample/input/9/ -o fpFactorUniqueOrc1Cdc6Rep2.txt 2> log/findPeaksFactorUniqueTaggedProteinIndividualSampleRepĺicate2.err > log/findPeaksFactorUniqueTaggedProteinIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/unique/taggedProtein/individualSample/chip/12/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/taggedProtein/individualSample/input/11/ -o fpFactorUniqueOrc1Cdc6Rep3.txt 2> log/findPeaksFactorUniqueTaggedProteinIndividualSampleRepĺicate3.err > log/findPeaksFactorUniqueTaggedProteinIndividualSampleRepĺicate3.out

# MULTIMAPPERS

findPeaks output/tagDirectory/multimappers/control/individualSample/chip/2/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/control/individualSample/input/1/ -o fpFactorMultimappersControlRep1.txt 2> log/findPeaksFactorMultimappersControlIndividualSampleRepĺicate1.err > log/findPeaksFactorMultimappersControlIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/multimappers/control/individualSample/chip/4/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/control/individualSample/input/3/ -o fpFactorMultimappersControlRep2.txt 2> log/findPeaksFactorMultimappersControlIndividualSampleRepĺicate2.err > log/findPeaksFactorMultimappersControlIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/multimappers/control/individualSample/chip/6/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/control/individualSample/input/5/ -o fpFactorMultimappersControlRep3.txt 2> log/findPeaksFactorMultimappersControlIndividualSampleRepĺicate3.err > log/findPeaksFactorMultimappersControlIndividualSampleRepĺicate3.out
findPeaks output/tagDirectory/multimappers/taggedProtein/individualSample/chip/8/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/taggedProtein/individualSample/input/7/ -o fpFactorMultimappersOrc1Cdc6Rep1.txt 2> log/findPeaksFactorMultimappersTaggedProteinIndividualSampleRepĺicate1.err > log/findPeaksFactorMultimappersTaggedProteinIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/multimappers/taggedProtein/individualSample/chip/10/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/taggedProtein/individualSample/input/9/ -o fpFactorMultimappersOrc1Cdc6Rep2.txt 2> log/findPeaksFactorMultimappersTaggedProteinIndividualSampleRepĺicate2.err > log/findPeaksFactorMultimappersTaggedProteinIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/multimappers/taggedProtein/individualSample/chip/12/ -style factor -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/taggedProtein/individualSample/input/11/ -o fpFactorMultimappersOrc1Cdc6Rep3.txt 2> log/findPeaksFactorMultimappersTaggedProteinIndividualSampleRepĺicate3.err > log/findPeaksFactorMultimappersTaggedProteinIndividualSampleRepĺicate3.out



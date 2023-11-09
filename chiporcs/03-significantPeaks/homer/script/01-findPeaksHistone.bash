# Finding peaks with --style histone
mkdir -p output/findPeaks/{unique,multimappers}

# findPeaks <tag directory> -style histone -size 500 -minDist 1000  -i <control tag directory> -o <path to output file>
# -gsize = genome size. According to the tritrypDB website (Data -> Organisms: Genome Info and Stats -> search for Trypanosoma cruzi CL Brener.
#          Esmeraldo-like = 32.53Mbp = 32530000
#          CL Brener = 36.03Mbp = 36030000

# UNIQUE

findPeaks output/tagDirectory/unique/control/individualSample/chip/2/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/control/individualSample/input/1/ -o fpUniqueControlRep1.txt 2> log/findPeaksHistoneUniqueControlIndividualSampleRepĺicate1.err > log/findPeaksHistoneUniqueControlIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/unique/control/individualSample/chip/4/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/control/individualSample/input/3/ -o fpUniqueControlRep2.txt 2> log/findPeaksHistoneUniqueControlIndividualSampleRepĺicate2.err > log/findPeaksHistoneUniqueControlIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/unique/control/individualSample/chip/6/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/control/individualSample/input/5/ -o fpUniqueControlRep3.txt 2> log/findPeaksHistoneUniqueControlIndividualSampleRepĺicate3.err > log/findPeaksHistoneUniqueControlIndividualSampleRepĺicate3.out
findPeaks output/tagDirectory/unique/taggedProtein/individualSample/chip/8/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/taggedProtein/individualSample/input/7/ -o fpUniqueOrc1Cdc6Rep1.txt 2> log/findPeaksHistoneUniqueTaggedProteinIndividualSampleRepĺicate1.err > log/findPeaksHistoneUniqueTaggedProteinIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/unique/taggedProtein/individualSample/chip/10/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/taggedProtein/individualSample/input/9/ -o fpUniqueOrc1Cdc6Rep2.txt 2> log/findPeaksHistoneUniqueTaggedProteinIndividualSampleRepĺicate2.err > log/findPeaksHistoneUniqueTaggedProteinIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/unique/taggedProtein/individualSample/chip/12/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/unique/taggedProtein/individualSample/input/11/ -o fpUniqueOrc1Cdc6Rep3.txt 2> log/findPeaksHistoneUniqueTaggedProteinIndividualSampleRepĺicate3.err > log/findPeaksHistoneUniqueTaggedProteinIndividualSampleRepĺicate3.out

# MULTIMAPPERS

findPeaks output/tagDirectory/multimappers/control/individualSample/chip/2/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/control/individualSample/input/1/ -o fpMultimappersControlRep1.txt 2> log/findPeaksHistoneMultimappersControlIndividualSampleRepĺicate1.err > log/findPeaksHistoneMultimappersControlIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/multimappers/control/individualSample/chip/4/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/control/individualSample/input/3/ -o fpMultimappersControlRep2.txt 2> log/findPeaksHistoneMultimappersControlIndividualSampleRepĺicate2.err > log/findPeaksHistoneMultimappersControlIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/multimappers/control/individualSample/chip/6/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/control/individualSample/input/5/ -o fpMultimappersControlRep3.txt 2> log/findPeaksHistoneMultimappersControlIndividualSampleRepĺicate3.err > log/findPeaksHistoneMultimappersControlIndividualSampleRepĺicate3.out
findPeaks output/tagDirectory/multimappers/taggedProtein/individualSample/chip/8/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/taggedProtein/individualSample/input/7/ -o fpMultimappersOrc1Cdc6Rep1.txt 2> log/findPeaksHistoneMultimappersTaggedProteinIndividualSampleRepĺicate1.err > log/findPeaksHistoneMultimappersTaggedProteinIndividualSampleRepĺicate1.out
findPeaks output/tagDirectory/multimappers/taggedProtein/individualSample/chip/10/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/taggedProtein/individualSample/input/9/ -o fpMultimappersOrc1Cdc6Rep2.txt 2> log/findPeaksHistoneMultimappersTaggedProteinIndividualSampleRepĺicate2.err > log/findPeaksHistoneMultimappersTaggedProteinIndividualSampleRepĺicate2.out
findPeaks output/tagDirectory/multimappers/taggedProtein/individualSample/chip/12/ -style histone -F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i output/tagDirectory/multimappers/taggedProtein/individualSample/input/11/ -o fpMultimappersOrc1Cdc6Rep3.txt 2> log/findPeaksHistoneMultimappersTaggedProteinIndividualSampleRepĺicate3.err > log/findPeaksHistoneMultimappersTaggedProteinIndividualSampleRepĺicate3.out



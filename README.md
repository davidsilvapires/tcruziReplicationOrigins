<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/f/fc/Logo_Instituto_Butantan_horizontal.svg/800px-Logo_Instituto_Butantan_horizontal.svg.png" height="100">

# _Trypanosoma cruzi_ Replication Origins

> "At the replication origins, the molecular symphony of DNA duplication unfolds, weaving the narrative of genetic continuity." - [ChatGPT]

[![Version](https://img.shields.io/badge/version-1.0.0-brightgreen.svg)](https://semver.org)

## Summary

1. [Overview](#overview)
2. [Install](#install)
3. [Pipelines](#pipelines)
4. [Configuration](#configuration)
5. [Usage](#usage)
6. [Examples](#examples)
7. [Contribution](#contribution)
8. [License](#license)
9. [Additional Notes](#notes)

## Overview

_Trypanosoma cruzi_, the parasite responsible for Chagas disease, has a genome containing numerous multigenic families that encode virulent antigens within a dynamically evolving genomic compartment called Disruptive. DNA replication is an accurate biological process, and deviations can lead to mutations and alterations in chromosomal and gene copy numbers. In this manner, understanding the dynamics of  _T. cruzi_ DNA replication can contribute to better elucidating its relationship with the success of the parasite's infection and the evolution of the species. Briefly, we conducted a computational analysis using data from ChIP-seq of Orc1Cdc6, Nanopore sequencing of BrdU-incorporated DNA (D-NAscent), and MFA-seq of epimastigote _T. cruzi_ to locate replication origins. By employing cutting-edge high-throughput and single-molecule analyses, we have identified three distinct categories of origins that compose an atlas: Predominant, Flexible, and Dormant Orc1Cdc6-dependent origins, in addition to Orc1Cdc6-independent origins.

This repository contains scripts relative to experiments published at the following paper: **Integrating high-throughput analysis to create an atlas of replication origins in _Trypanosoma cruzi_ in the context of genome structure and variability** (_to be submitted_).

We present three different pipelines:
1. **Chiporcs:** ChIP-seq pipeline of Orc1Cdc6 in _Trypanosoma cruzi_.
2. **D-NAscent:** Nascent DNA pipeline.
3. **Atlas:** The construction of a dataset that combines the results of the implemented pipelines.

## Install

To install the dependencies, use the following command:

```bash
sudo apt install bedtools hisat2 samtools
curl --output /usr/local/bin/faCount https://github.com/ENCODE-DCC/kentUtils/blob/master/bin/linux.x86_64/faCount
```

Besides the above packages and programs, there are other software that have to be installed according to the instructions at their sites:
* [Homer](http://homer.ucsd.edu/homer/download.html)


## Pipelines

This project is based on three main pipelines to identify and analyze the replication origins identified by different techniques (ChIP-seq, D-NAscent and MFA-seq). Firstly, in this project, we identify the DNA replication origins of _T. cruzi_ through three distinct approaches, as mentioned above. In this manner, a dataset of DNA replication origins identified by the MFA-seq technique was utilized, which had previously been published by our group. Next, Orc1Cdc6 binding sites were identified through the ChIP-seq assay, using the pipeline created by our group and will be detailed subsequently. Finally, DNA replication origins were identified through the methodology of detecting BrDU-incorporated nascent DNA molecules, followed by MinION sequencing and the D-NAscent software, which enable the determination of the direction of the replication fork in individual DNA molecules.

The three pipelines are coded in the following shared files:
1. chiporcs
2. D-NAscent
3. atlas

Each pipeline is composed of steps that should be executed in a linear fashion way. The order of the steps are indicated by numbers. If two steps can be executed in parallel, the number is followed by a letter (e.g.: `5a`, `5b`).


### Chip-seq of ORCs pipeline

Providing further detail on the pipeline for the analysis of Orc1Cdc6 binding sites (ChIP-seq): The sequenced reads were initially subjected to quality visualization using FastQC software, version 0.11.9. Subsequently, adapters were removed, and a quality filter (`LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:15`) was applied using Trimmomatic software, version 0.39. The initial and final bases with quality less than 5 were automatically discarded, an average quality of greater than or equal to 15 was applied to a 4bp window, and reads that did not meet the average were discarded. The filtered reads were mapped to the genomes of _T. cruzi_ CL Brener, _T. cruzi_ CL Brener Esmeraldo-like, and _T. cruzi_ CL Brener Non-Esmeraldo-like, version 32, available on TritrypDB, using `bwa aln` software, version 0.7.17. `bowtie2` (`--local -D 25 -R 4 -N 1 -L 19 -i S1,0.40 --phred33 --n-ceil L,0,0.15 --dpad 15 --gbar --time --no-unal`), version 8.3.1, or `hisat2` (`--n-ceil L,0,0.15 --score-min L,0.6,-0.6, --no-spliced-alignment --no-unal`), version 5.4.0.

The search for enriched regions in ChIP-seq experiments was conducted using the Homer software, version 4.11, employing the following commands: i) `makeTagDirectory` (`-sspe -mapq 10 -tbp 1`, `-unique` or `-keepAll`), for the creation of tag directories containing reads with a quality greater than or equal to 10, removing potential PCR amplification duplications and selecting uniquely mapped or multipappers reads; ii) `findPeaks` (`style histone` or `factor`, `-F 2 -fdr 0.05 -size 50 -minDist 150 -gsize 25827529 -i`), which identifies peaks in broad or narrow styles, with a fold change greater than or equal to 2, a false positive rate of 5%, in a 50bp window with a minimum distance between adjacent peaks of 150bp, normalizing ChIP samples by their respective input; iii) `mergePeaks` (`-d 340 -matrix -gsize 25827529`), to merge common peaks among the three replicates in a 340 bp window, using a hypergeometric distribution that calculates the statistical significance of peak co-occurrence. The effective genome size, disregarding 'N' bases in the reference genome, was used, calculated using the `faCount` tool. Once identified, significant peaks from the control lineage were subtracted from those of the Orc1Cdc6 lineages using the subtract option in the Bedtools software. The genomic coordinates of Orc1Cdc6 binding sites were organized into a bed-format file.

### D-NAscent pipeline

For the assembly of the dataset of origins identified by the D-NAscent method, we employed a pipeline developed by the group using the D-NAscent software. Detailing the steps, the reads generated by MinION sequencing in the FAST5 format, containing unique electrical signals for each nucleotide, were converted to base format using the Guppy software (`--flowcell FLO-MIN106 --kit SQK-LSK109 --qscore_filtering --trim_strategy dna`), version 4.2.2. Subsequently, the reads were mapped to the genome of _T.cruzi_ CL Brener Esmeraldo-like version 30 using the Minimap2 software (`-ax map-ont`), version 2.17. The mapped reads in the SAM format were converted to BAM using the Samtools software, version 1.12. Then, the original reads in FAST5 format, and the mapped reads (BAM), were used as input for BrdU detection through the D-NAscent software, version 2.0.2, using the following commands: i) `detect` (`--quality 20 -–length 1000`) for BrdU signal detection, and ii) `forkSense`, for the detection of the replication fork direction, based on BrdU decay along the read. Next, with the filtered reads, the coordinates of each read containing a fork direction probability greater than or equal to 70% were selected. These data were plotted in a table and used to i) determine the regions of replication initiation or termination, and ii) identify potential regions of conflict between replication and transcription machineries. Thus, for the detection of replication origins, only forks with divergent directions in the same read were considered. The intervals between these forks were determined as the coordinates of the replication start sites. Once the coordinates of the origins in the genome were determined, the midpoint of each origin was calculated, and a 3 kb window was applied to establish a maximum size for the replication origin regions. These coordinates were organized into a BED file containing the genomic coordinates of each identified origin.


### Atlas composition


## Configuration

To configure the pipeline and adapt it to your specific necessities...

You can customize the behaviour of the pipeline by editing the configuration file at `config/pipeline_config.yml`.

## Usage

Our pipeline is basically written in Bash, only the final steps that involve the generation of plots are written in R.

To execute the pipeline, run the following command:

```bash
mkdir log
script/chiporcsDriver.bash 2>&1 | tee log/chiporcsDriver.log
```


## Examples


## Contribution

Suggestions are welcome! Please, send your ideas to the following email address:


## License

This project is licensed under [GPL](LICENSE).


## Additional notes

- [Important note about the pipeline]
- [Other important informations]

<!--
This site was built using [GitHub Pages](https://pages.github.com/).

- George Washington
* John Adams
+ Thomas Jefferson

```bash
echo "Código em Bash."
```

Here is a simple footnote[^1].

A footnote can also have multiple lines[^2].

[^1]: My reference.
[^2]: To add line breaks within a footnote, prefix new lines with 2 spaces.
  This is a second line.

[Main directory](chiporcs/)

> [!NOTE]
> Useful information that users should know, even when skimming content.

> [!TIP]
> Helpful advice for doing things better or more easily.

> [!IMPORTANT]
> Key information users need to know to achieve their goal.

> [!WARNING]
> Urgent info that needs immediate user attention to avoid problems.

> [!CAUTION]
> Advises about risks or negative outcomes of certain actions.
-->

<!-- This content will not appear in the rendered Markdown -->

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

_Trypanosoma cruzi_, the parasite responsible for Chagas disease, has a genome containing numerous multigenic families that encode virulent antigens within a dynamically evolving genomic compartment called Disruptive. DNA replication is an accurate biological process, and deviations can lead to mutations and alterations in chromosomal and gene copy numbers. In this manner, understanding the dynamics of  _T. cruzi_ DNA replication can contribute to better elucidating its relationship with the success of the parasite's infection and the evolution of the species. Briefly, we conducted a computational analysis using data from ChIP-seq of Orc1Cdc6, Nanopore sequencing of BrdU-incorporated DNA (D-Nascent), and MFA-seq of epimastigote _T. cruzi_ to locate replication origins. By employing cutting-edge high-throughput and single-molecule analyses, we have identified three distinct categories of origins that compose an atlas: Predominant, Flexible, and Dormant Orc1Cdc6-dependent origins, in addition to Orc1Cdc6-independent origins.

This repository contains scripts relative to experiments published at the following paper: **Integrating high-throughput analysis to create an atlas of replication origins in _Trypanosoma cruzi_ in the context of genome structure and variability** (_to be submitted_).

We present three different pipelines:
1. **Chiporcs:** ChIP-seq pipeline of Orc1Cdc6 in _Trypanosoma cruzi_.
2. **DNAscent:** Nascent DNA pipeline.
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
Three pipelines are coded in the shared files:
1. chiporcs
2. dnascent
3. atlas

Each pipeline is composed of steps that should be executed in a linear fashion way. The order of the steps are indicated by numbers. If two steps can be executed in parallel, the number is followed by a letter (e.g.: `5a`, `5b`).

### Chip-seq of ORCs pipeline

### DNAscent pipeline

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

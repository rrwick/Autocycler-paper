<p align="center"><img src="https://github.com/rrwick/Autocycler/wiki/images/logo.png" alt="Autocycler logo" width="90%"></p>

This repository contains the supplementary data for our paper describing Autocycler:<br>[Wick RR, Stinear TP. Autocycler: long-read consensus assembly for bacterial genomes. bioRxiv. 2025. doi:10.1101/2025.05.12.653612](https://doi.org/10.1101/2025.05.12.653612)

If you are looking for Autocycler itself (i.e. the software, not the paper), please see the [Autocycler GitHub repository](https://github.com/rrwick/Autocycler).

This repository contains:
* [`figures/`](figures): all main-text and supplementary figures.
* [`tables.xlsx`](tables.xlsx): all supplementary tables for the paper.
* [`commands.md`](commands.md): methods and exact commands used to generate the data shown in Figure 2 and the supplementary tables.
* [`assess_assembly.py`](assess_assembly.py): the Python script used to quantify assembly accuracy.
* [`plots.Rmd`](plots.Rmd): R code for generating the plots in Figure 2 and Figure S1.

The datasets used in this study (reads, assemblies and reference genomes) are too large to host on GitHub but are available from [this figshare repository](https://figshare.unimelb.edu.au/projects/Autocycler/247142). While the five isolates are also available via the SRA (accessions listed in Table S1), the ONT reads were re-basecalled using the latest version of Dorado. The figshare repository provides the exact post-QC read sets used in this paper.

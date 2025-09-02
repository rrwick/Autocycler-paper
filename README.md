<p align="center"><img src="https://github.com/rrwick/Autocycler/wiki/images/logo.png" alt="Autocycler logo" width="90%"></p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16916187.svg)](https://doi.org/10.5281/zenodo.16916187)

This repository contains the supplementary data for our manuscript describing Autocycler:<br>
[Wick RR, Howden BP, Stinear TP (2025). Autocycler: long-read consensus assembly for bacterial genomes. _Bioinformatics_. doi:10.1093/bioinformatics/btaf474.](https://doi.org/10.1093/bioinformatics/btaf474)

If you are looking for Autocycler itself (i.e. the software, not the manuscript), please see the [Autocycler GitHub repository](https://github.com/rrwick/Autocycler).

This repository contains:
* [`manuscript.pdf`](manuscript.pdf): the Autocycler manuscript.
* [`supp_figures.pdf`](supp_figures.pdf): supplementary figures for the manuscript.
* [`supp_tables.xlsx`](supp_tables.xlsx): supplementary tables for the manuscript.
* [`figures/`](figures): directory with individual PDFs for all main-text and supplementary figures.
* [`commands.md`](commands.md): methods and exact commands for all analyses in the manuscript.
* [`assess_assembly.py`](assess_assembly.py): the Python script used to quantify assembly accuracy.
* [`plots.Rmd`](plots.Rmd): R code for generating the figures.
* [`Autocycler-0.5.1.tar.gz`](Autocycler-0.5.1.tar.gz): source code for the version of Autocycler used in the manuscript.

The datasets used in this study (reads, assemblies and reference genomes) are too large to host on GitHub but are available from [this figshare repository](https://figshare.unimelb.edu.au/projects/Autocycler/247142). While the five isolates are also available via the SRA (accessions listed in Table S1), the ONT reads were re-basecalled using the latest version of Dorado. The figshare repository provides the exact post-QC read sets used in this manuscript.

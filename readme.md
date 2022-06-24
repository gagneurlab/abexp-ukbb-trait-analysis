# Splicing Analysis

Check whether spliced variants have a similar effect on traits as pLoF variants.

The same traits are used as in Remo Monti, ... Christoph Lippert  https://www.biorxiv.org/content/10.1101/2021.05.27.444972v1

The traits are listed here: https://docs.google.com/spreadsheets/d/1cqVWj5MSExZUMMHKLlltVqdcjMJNxh9Py2_zKsLqqmA/edit#gid=454465379



currently executing this pipeline using:
```
snakemake -j 5 --cluster "sbatch -N 1 -n 16 --mem=64G" --use-conda
```


# wBuild Docs
## Documentation

Full documentation is available at https://wbuild.readthedocs.io

## Features

* Supports reproducible research
  * Append R-markdown scripts to the Snakemake pipeline
* Render the R scripts to a structured web page
* Snakemake rules written directly in the header of scripts
  * Dependencies updated automatically!

## Installation

Install using pip:

- `pip install wBuild`

See "Installation" tab in the documentation for more details.

## Getting started
  
* Navigate to an empty directory
* Run `wbuild demo`. This will create a wBuild demo project with various examples
* Explore the files in `Scripts/`
* Run `snakemake` to build the project
* Open `Output/html/index.html` in your web browser
  * There you will find useful and understandable examples of features and operations with wBuild

## Usage

* Navigate to the root of your project (new or existing)
* Run `wbuild init`
* Run `snakemake`

## GitHub

wBuild is an open-source software. The source-code is available at https://github.com/gagneurlab/wBuild.

## Credits

Leonhard Wachutka

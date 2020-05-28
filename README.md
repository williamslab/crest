# crest
Classification of RelationShip Types

## Introduction
CREST (**C**lassification of **R**elation**S**hip **T**ypes) is a tool that uses identity-by-descent (IBD) segments to classify second-degree relatives as avuncular, half-siblings, or grandparent/grandchild.

## Quick Start
Extract IBD information from your data with the program of your choice. We recommend using [IBIS](https://github.com/williamslab/ibis).
```
ibis your_data.bed your_data.bim your_data.fam -min_l 7 -mt 500 -er 0.04 -f IBIS_your_data
```
...

Run sex-inference
```
python3 CREST_sex_inference.py -i IBIS_your_data.seg -m your_genetic_map.simmap -b your_data.bim -o sex_inference_output
```
## Thorough Start
### Pre-CREST Data Generation and Curation
### Relationship Type Inference
### Parental Sex Inference
#### Command line arguments:

*  `-i` or `--input`

    * Name of the input file
    * Required argument

* `-o` or `--output`

    * Name of the output file
    * Required argument

* `-m` or `--map`

    * Name of a genetic map in .simmap format
    * Required argument

* `-b` or `--bim`

    * Name of the original .bim file for your data
    * Required argument

* `-w` or `--window`

    * Window size in kilobases
    * Optional argument, defaults to 500 kb
    * Must be an integer

* `-k` or `--keep`

    * List of pairs of samples to keep
    * Optional argument, defaults to `None`

#### Output file format
The format of the output file is
```
id1 id2 number_of_segments LOD_if_GP LOD_if_HS
```
The last two columns are quasi-LOD scores under the HS and GP models, where LOD = log10 (p(maternal) / p(paternal)). In other words,  positive scores indicate a pair is more probably maternal, while negative scores indicate a pair is more probably paternal.

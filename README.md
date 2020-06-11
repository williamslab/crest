# crest

## Introduction
CREST (**C**lassification of **R**elation**S**hip **T**ypes) is a tool that uses identity-by-descent (IBD) segments to classify second-degree relatives as avuncular, half-siblings, or grandparent/grandchild. It can also be used to infer the directionality of relationships and whether they are maternal or paternal related. 

## Quick Start
Follow these steps to get CREST results quickly and easily. All file names and directories in brackets should be replaced with names and directories of your choosing. CREST includes three parts: CREST_ratios to calculate the ratios of IBD sharing with mutual relatives, CREST_relationships to infer relationship types and directionality, and CREST_sex_inference to infer whether they are maternal or paternal related.
### Your Data
CREST, and IBIS, if you choose to use it, requires genotype data in a PLINK binary file format. Note that CREST currently uses autosomal IBD only, so if you have non-autosomal data, be sure to exclude it later on.
### Getting IBD Segments
We recommend using [IBIS](https://github.com/williamslab/ibis) to extract IBD information. It tends to infer contiguous segments, which is especially important for sex inference. 
Before running IBIS, we advise adding a genetic map to your .bim file. See the IBIS documentation [here](https://github.com/williamslab/ibis#Steps-for-running-IBIS):
```
./add-map-plink.pl [your data].bim [map directory]/genetic_map_GRCh37_chr{1..22}.txt > [your new data].bim
```

Then run IBIS itself:
```
ibis [your data].bed [your new data].bim [your data].fam -f [your IBIS data]
```
or if you rename the new .bim, you can supply all three at once:
```
ibis -b [your data] -f [your IBIS data]
```

### Run CREST_ratios

First, complie by running 
```
make
```
Then CREST_ratios takes in the .seg file and .coef file from IBIS as input. 
```
./crest_ratios -i [ibd segment].seg -r [relative list].coef -o [output prefix]
```

CREST_ratios will generate [output prefix].csv file with this format:
```
ID1 ID2 coverage_in_cM ratio1 ratio2
```
Details about other options see below.


### Run CREST_relationships

The CREST_relationships takes in the .csv output of CREST_ratios as the input file. It also needs the total map length in cM to calculate the coverage rate. The basic useage is 
```
./CREST_relationships -i [ratios].csv --total_len [total length of genome] -o [output prefix].csv 
```
Details about other options see below.

The total map length in cM is available in the .bim file. The maplen.awk script in [IBIS](https://github.com/williamslab/ibis) calculates this in the following way:
```
./maplen.awk [bim files ...]
```

CREST_relationships will generate [output prefix].csv file with this format:
```
ID1 ID2 inferred_type prob_gp prob_av prob_hs inferred_direction prob1 prob2
```
For the `inferred_type` column, 1 is for GP, 2 is for AV, and 3 is for HS. For the `inferred_direction` column, 0 means sample1 is genetically older than sample2 and 1 means ample1 is genetically younger than sample2.

### Run CREST_sex_inference

For sex-inference, you will need to convert sex-specific genetic maps of your choosing to a .simmap format file.
Information on how to do this can be found [here](https://github.com/williamslab/ped-sim#map-file):
```bash
bash
wget https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/raw/master/Refined_genetic_map_b37.tar.gz
tar xvzf Refined_genetic_map_b37.tar.gz
printf "#chr\tpos\tmale_cM\tfemale_cM\n" > [your map].simmap
for chr in {1..22}; do
  paste Refined_genetic_map_b37/male_chr$chr.txt Refined_genetic_map_b37/female_chr$chr.txt \
    | awk -v OFS="\t" 'NR > 1 && $2 == $6 {print $1,$2,$4,$8}' \
    | sed 's/^chr//' >> [your map].simmap;
done
```
Now you are ready to run the sex-inference script:
```
CREST_sex_inference.py -i [your IBIS data].seg -m [your map].simmap -b [your (new) data].bim -o [sex inference output]
```
## Thorough Start
### Pre-CREST Data Generation and Curation
The .seg file contains the IBD segments for all samples and .coef file is a list of relative pairs with inferred relatedness. If use IBD segment and relative information from other tools, please make sure to have the same data format. The useful information from .seg file includes sample1, sample2, chromosome, IBD type, genetic start position, genetic end position, genetic length in 1st, 2nd, 3rd, 6th, 7th, 8th, 9th columns. The information sample1, sample2, IBD2_fraction, and degree_of_relatedness in 1st, 2nd, 4th, 6th column of .coef file is used.  

### Relationship Type Inference
CREST_relationships requires python packages including sklearn, numpy. Please make sure you have them installed.
#### Command line arguments for CREST_ratios:
CREST_ratios also have following options:
* `--ibd2 <value between 0 and 1>` : the threshold of IBD2 ratios to exclude relatives from second-degree relatives. The default value is 0.02

* `--max_degree <integer larger than 2>`: the upper bound of degree of relatedness for mutual relatives to the pair. The default value is 6. 

* `--cluster_thres <value>` : the threshold of genetic length in cM to cluster mutual relatives. If shared IBD length between mutual relatives is large than this threshold, then they are considered as relatives to each other too. The default value is 10. 

#### Command line arguments for CREST_relationships:
CREST_relationships have following options:
* `-i` or `--input` : the .csv file contains the ratios information

* `-o` or `--output` :  the .csv file of inferred relationship types and directionality. The default name is out.csv.

* `--total_len <value>` : the total genetic length in cM

* `--models_type` : name of trained models to infer relationship types. The default is type_clf.pickle. It will also be used to store new trained models if `--train` is enabled.

* `--models_direction` : name of trained models to infer relationship directionality. The default is direction_clf.pickle. 

* `--start <value>` : the minimum coverage rate to infer relationships. If one pair has coverage rate smaller than this value, the pair will not be inferred. The default value is 0.025.

* `--end <value>` : the coverage rate to merge models. For pairs with coverage rates larger than this value, they will use the same trained model as the ones with this coverage rate. The default value is 0.2.

* `--inv <value>` : the window size of coverage rate to train different models. It will be used to divide the coverage rate from the start to the end into different windows. The default value is 0.025.

* `--prior [p1 p2 p3]` : the prior probability of three types of relationships. The default is 0.3333 for each type. If you believe the probabilities of three relationship types are equal, specify the prior probability and make sure the sum of these three values equals to 1. If you train new models and believe the testing data has the same distribution with your trained models, then there is no need to change the default.  

* `--train` : train new models with labeled data

* `--labels` : the file of labels for training dataset	

### Parental Sex Inference
#### Command line arguments:

* `-i` or `--input`: name of the input file (required)

* `-o` or `--output`: name of the output file (required)

* `-m` or `--map`: name of a genetic map in .simmap format (required)

* `-b` or `--bim`: name of the .bim file from your IBD data (required)

* `-w` or `--window`: window size in kilobases (optional, integer, defaults to 500 kb)
    * Windows are symmetric about the IBD segment ends, so a 500 kb window extends 250 kb in both directions.

* `-k` or `--keep`: list of sample pairs to keep for sex-inference analysis (optional, defaults to `None`)


#### Output file format
The format of the output file is
```
id1 id2 segment_number GP_lod HS_lod
```
The last two columns are quasi-LOD scores under the HS and GP models, where LOD = log10 (p(maternal) / p(paternal)). In other words,  positive scores indicate a pair is more probably maternal, while negative scores indicate a pair is more probably paternal.

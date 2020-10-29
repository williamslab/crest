# CREST

CREST (**C**lassification of **RE**lation**S**hip **T**ypes) is a tool that uses identical by descent (IBD) segments to classify second degree relatives as grandparent/grandchild (GP), avuncular (AV), or half-siblings (HS). It can also be used to infer the directionality of relationships and whether they are maternally or paternally related.

## Quick Start

### 0. Genetic Data
CREST requires a PLINK format .bim file, and IBIS, if you choose to use it, requires genotype data in PLINK binary format. [PLINK](https://www.cog-genomics.org/plink2/) converts to this format with the `--make-bed` option.

### 1. Getting IBD Segments
We recommend using [IBIS](https://github.com/williamslab/ibis) to extract IBD segments. It tends to infer contiguous segments, which is especially important for sex inference.
Before running IBIS, insert a genetic map to your .bim file with the command below. The [IBIS documentation](https://github.com/williamslab/ibis#Steps-for-running-IBIS) provides links to the HapMap genetic map.
```
./add-map-plink.pl [your data].bim [map directory]/genetic_map_GRCh37_chr{1..22}.txt > [your new data].bim
```

Then run IBIS itself:
```
ibis [your data].bed [your new data].bim [your data].fam -f [your IBIS output] -printCoef
```
or if you rename the new .bim, you can supply all three at once:
```
ibis -b [your data] -f [your IBIS output] -printCoef
```

### 2. Run `crest_ratio`

First, compile by running
```
make
```
Then run `crest_ratio` with the .seg and .coef files from IBIS as input.
```
./crest_ratio -i [your IBIS output].seg -r [your IBIS output].coef -o [ratio output prefix]
```

`crest_ratio` will generate [ratio output prefix].csv file with this format:
```
ID1 ID2 is_pc coverage_in_cM ratio1 ratio2
```
Details about other options are [below](#command-line-arguments-for-CREST_ratios).


### 3. Run `CREST_relationships.py`
`CREST_relationships.py` used Python 3.8 and requires python packages including sklearn 0.23, numpy. Please make sure you have them installed. Other versions may end up with warning or error messages.
`CREST_relationships.py` takes in the .csv from `crest_ratio` as input. It also needs the total map length in cM to calculate the coverage rate. The basic usage is
```
./CREST_relationships.py -i [ratio output prefix].csv --total_len [total length of genome] -o [relationships output prefix]
```

The total map length in cM is available in the .bim file. The maplen.awk script in [IBIS](https://github.com/williamslab/ibis) calculates this in the following way:
```
./maplen.awk [bim files ...]
```
Details about other options to `CREST_relationships.py` are [below](#command-line-arguments-for-CREST_relationships).

`CREST_relationships.py` will generate [relationship output prefix].csv file with this format:
```
ID1 ID2 inferred_class inferred_type prob_gp prob_av prob_hs inferred_direction prob1 prob2
```
The `inferred_class` column is the numeric label for `inferred_type`. Here 0 is for PC, 1 is for GP, 2 is for AV, and 3 is for HS. For the `inferred_direction` column, 0 means sample1 is genetically older than sample2 and 1 means sample1 is genetically younger than sample2.

### 4. Run `CREST_sex_inference.py`

For sex-inference, you will need to put sex-specific genetic maps into a .simmap format file.
Information on how to do this can be found [here](https://github.com/williamslab/ped-sim#map-file); to use the [Bhérer et al. (2017)](http://dx.doi.org/10.1038/ncomms14994) genetic maps, run:

```bash
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
./CREST_sex_inference.py -i [your IBIS output].seg -m [your map].simmap -b [your (new) data].bim -o [sex inference output]
```
The above will compute LOD scores for _all_ pairs. To limit to only second degree relatives, you can use the `-k` (keep) option as:
```
./CREST_sex_inference.py -i [your IBIS output].seg -m [your map].simmap -b [your (new) data].bim -o [sex inference output] -k <(awk '$6 == 2' [your IBIS output].coef)
```
The file redirection `<( ... )` syntax works in bash. Other shells can run `awk '$6 == 2' [your IBIS output].coef > second_degree.coef`, and then run `CREST_sex_inference.py` with `-k second_degree.coef`. If you have already performed relationship inference, you can also use the [relationship output prefix].csv as a keep file, though you may want to extract only the GP or HS pairs first. 

## Thorough Start
### Pre-CREST Data Generation and Curation
The IBIS format .seg file contains the IBD segments for all samples and .coef file is a list of relative pairs with inferred relatedness. If you want to use IBD segments and relative information from other tools, please make sure to have the same data format. The useful information from .seg file includes ID1, ID2, chromosome, IBD type, genetic start position, genetic end position, genetic length in 1st, 2nd, 3rd, 6th, 7th, 8th, 9th columns correspondingly. The information ID1, ID2, IBD2\_fraction, and degree\_of\_relatedness in 1st, 2nd, 4th, 6th columns of .coef file is used.  

### Relationship Type Inference
#### Command line arguments for CREST\_ratios:
CREST_ratios also have following options:
* `--ibd2 <value between 0 and 1>` : the threshold of IBD2 ratios to exclude some relatives from second-degree relatives, since GP, AV, and HS are expected to share 0 ibd2 proportion. The default value is 0.02. 

* `--max_degree <integer larger than 2>`: the upper bound of degree of relatedness for mutual relatives to the pair. The default value is 6. 

* `--pc`: to predict the directionality of PC. 
 
* `--ibd0 <value between 0 and 1>` : the threshold of IBD0 ratios to distinguish PC. The default value is 0.1.

* `--kinship_lw <value between 0 and 1>` : the lower bound of kinship to choose second degree pairs. The default value is 0.0883883.

* `--kinship_up <value between 0 and 1>` : the upper bound of kinship to choose second degree pairs. The default value is 0.1767767. 

* `--kinship_rel_lw <value between 0 and 1>` : the lower bound of kinship to choose mutual relatives. The default value is 0.0055243, corresponding to the lower bound of sixth degree relatives. 

* `--kinship_rel_up <value between 0 and 1>` : the upper bound of kinship to choose mutual relatives. The default value is 0.0883883, corresponding to the upper bound of third degree relatives. 

#### Command line arguments for CREST\_relationships:
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

* `--labels` : the file containing labels for training dataset. Here 1 is for GP, 2 is for AV, and 3 is for HS.	

### Parental Sex Inference
#### Command line arguments:

* `-i` or `--input`: name of the input file (required)
    * Input should be in IBIS format, where rows correspond to segments and columns (separated by whitespace) are: 'id1', 'id2', 'chromosome', 'start (bp)', 'stop (bp)'

* `-o` or `--output`: name of the output file (required)

* `-m` or `--map`: name of a genetic map in .simmap format (required)

* `-b` or `--bim`: name of the .bim file from your IBD data (required)

* `-w` or `--window`: window size in kilobases (optional, integer, defaults to 500 kb)
    * Windows are symmetric about the IBD segment ends, so a 500 kb window extends 250 kb in both directions.

* `-k` or `--keep`: list of sample pairs to keep for sex-inference analysis (optional, defaults to `None`)
    * The fields of this file may be separated by whitespace or commas, as long as the first two columns correspond to: 'id1', 'id2'
'

#### Output file format
The format of the output file is
```
id1 id2 segment_number GP_lod HS_lod
```
The last two columns are quasi-LOD scores under the HS and GP models, where LOD = log10 (p(maternal) / p(paternal)). In other words,  positive scores indicate a pair is more probably maternal, while negative scores indicate a pair is more probably paternal.

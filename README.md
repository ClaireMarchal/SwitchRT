# SwitchRT
R script for repli-seq data to assess replication timing (RT) changes between two conditions with mulitple replicate. For each genomic location, the q-value can be used to call which genomic regions have a RT significantly different between the two conditions (recommanded q-value threshold: 0.01).


## Input
- A tab delimited data file with header, containing the RT values from two or more conditions, with multiple replicates per condition, with each row correspond to a genomic location, the first 3 columns are the chromosome, start and stop positions for the genomic location, and each suplementary column are RT values from one sample.
- the name of the output file to be created
- the first set of indexes corresponding to the column with the first condition
- the second set of indexes corresponding to the column with the second condition

## Output
A tab delimited data file containing the input file with the genomic coordinates (3 first columns), the RT values for each sample selected, several statistics about the comparisons, and the q-value.

## Exemple

`Rscript call_regions.R input_file.txt output_file.txt 4,5,6 7,8,9`
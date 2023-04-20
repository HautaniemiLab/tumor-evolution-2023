# Mutational Signatures

Code and commands to perform signature analysis.

## Data

Data for the most part is acquired using the code provided in this repository,
see below. SigProfiler reference signatures have to be manually downloaded from
Synapse (registration required), and are assumed to be put into
`data/raw/reference_signatures/`.

## Performing analysis

First, prepare the resources once you have the reference genome `GRCh38.d1.vd1.fa`
and its index in `data/raw/reference_genomes/`:

`src/s000_prepare_resources.sh`

You can then perform the analysis by using the commands in `commands.discovery.sh`
(discovery cohort) and `commands.validation.sh` (validation cohort).

Note that you will also need the inputs in the correct paths.

## Directory structure

Code and data are split such that code are under `src/`, data is in `data/`
while results are within `results/`. Each of these have further logical
substructure.


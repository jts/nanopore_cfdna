# Deconvolution Workflow
This repository contains a [nextflow](nextflow.io)  workflow for deconvoluting a complex mixture of DNA into mixture proportions using nanopore methylation calls. A task useful for cell-free DNA analysis. Input files can either start from raw `.fast5` files or basecalled `.bam` files.
![Deconvolution Pipeline](DeconvolutionPipeline_diagram.png)

## External Dependencies
1. [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis)
2. [nanomix](https://github.com/jts/nanopolish)
3. [mbtools](https://github.com/jts/mbtools)


## Run
Make a new directory for your run. Inside that directory make another directory called data:
```
mkdir dev1
cd dev1
mkdir data
```
Then inside the data directory make symbolic links to the directories where the nanopore runs are kept.
```
cd data
ln -s <PATH_TO_DATA>
```
Then whilst in the run directory execute the nextflow script with the following command which takes as input the path to this repository.
```
cd ..
nextflow run <.../nanopore_cfdna/>
```
If you have multiple cores available use the `--threads` parameter to specify how many you would like to use.

## Results
The nextflow pipeline automatically creates a `work` directory and `results` directories. Inside you will find:
- `sample.read_modifications.tsv`: nanopolish methylation calling for reads
- `sample.referenceAtlas.region_modifications.tsv`: nanopolish methylation calling for each region in the deconvolution reference atlas
- `sample.bam`: alignment against reference
- `deconv_output.referenceAtlas.method.tsv`: Methylation deconvolution output. Estimation of proportion of cell type heterogeneity for each sample.
- `deconv_output.referenceAtlas.method.png`: Stacked bar plot of cell type heterogeneity.


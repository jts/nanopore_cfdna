# Simpson Lab cfDNA Analysis Pipeline
This pipeline follows a suite of functions to compute fragmentation and methylation patters of cfDNA samples

## External Dependencies
1. [Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis)
2. [Nanopolish](https://github.com/jts/nanopolish)
3. [mbtools](https://github.com/jts/mbtools)
4. [Reriomodels](https://github.com/nanoporetech/rerio)

## Setup
1. Install the environment with [Conda](https://docs.conda.io/en/latest/miniconda.html)
```
conda env create --file environment.yaml
```
2. Download the archive and the above dependencies with `git clone`
	- when downloading Nanopolish, install the `methylation_bam` version
	```
	git clone --recursive https://github.com/jts/nanopolish.git
	cd nanopolish
	git checkout methylation_bam
	make
	```
3. Edit `nextflow.config` to point the relevant parameters to where you installed them. A reference sequence is also required here. HDF5 plugin is installed with nanopolish. These are the following lines that need to be edited:
```
...
env.HDF5_PLUGIN_PATH="/.mounts/labs/simpsonlab/users/jsimpson/code/cfdna_pipeline/etc/"
...
guppy = "/.mounts/labs/ont/software/guppy-5.0.11/bin/guppy_basecaller"
reference = "/.mounts/labs/simpsonlab/data/references/GRCh38_no_alt_analysis_set.GCA_000001405.15.fna"
nanopolish = "/.mounts/labs/simpsonlab/users/jbroadbent/software/nanopolish/nanopolish"
mbtools = "/.mounts/labs/simpsonlab/users/jsimpson/code/mbtools/target/release/mbtools"
reriomodels = "/.mounts/labs/simpsonlab/users/jbroadbent/software/rerio/basecall_models"
...
```

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
The nextflow pipeline automatically creates a `work` directory and `results` directory for each sample. Inside the results directory you will:
- `sample.fragmentation_ratios.tsv`: Ratio of short to long reads in 5Mb bins. More info [here](https://www.nature.com/articles/s41467-021-24994-w). Fragmentome can be plotted with `plot_fragmentome.r`
- `sample.bamstats.tsv`: *Pomoxis* alignment statistics
- `sample.read_modifications.tsv`: nanopolish methylation calling for reads
- `sample.reference_modifications.tsv`: nanopolish methylation calling for each coordinate on the reference
- `sample.cpgfreq.csv`: Modification frequency of CpG sites uses in methylation atlas deconvolution. Can be used as input [here](https://github.com/nloyfer/meth_atlas).
- `sample.bam`: alignment against reference



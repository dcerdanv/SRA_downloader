# SRA_downloader
Small snakemake pipeline that automates SRA objects downloads; and extract and compresses their FASTQ files  

[Pipeline status](https://github.com/dcerdanv/sra_downloader/commits/main)

## Introduction
**SRA_downloader** is a **[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline** that automates the download of data from the [INSDC Sequence Read Archives](https://www.ncbi.nlm.nih.gov/sra/). It allows to store both sra objects and FASTQ files from a sample. FASTQ files can be saved compressed in GZ format or uncompressed.

Tools recommended by the SRA are used to download and access the data, all belonging to the [sra-tools](https://github.com/ncbi/sra-tools) package. Additionally [pigz](https://zlib.net/pigz/) is used to compress the FASTQ files. It is a parallel implementation of gzip that allows the compression time to be significantly reduced.

The pipeline makes use of Snakemake's integration with the [conda](https://docs.conda.io/en/latest/) package manager, to automatically take care of software requirements and dependencies.

I've built SRA_downloader to improve the speed of the data collection, with the aim of allowing the user to adjust the pipeline to different experiments using configuration parameters.
When in doubt, the default parameters were calculated to offer a good starting configuration.


## Author
 * Daniel Cerdán-Vélez

 
## Setup

The setup of the pipeline consists of the modifying two configuration files, one sets the desired parameters and the location of the output, log, and tmp files; and the other lists which samples you want to download.
A general description of these files follows. See the *Usage* section for more details.

### Configuration files

* **config.yaml** contains all pipeline parameters.
* **samples.tsv** contains the samples list.


## Usage 

### 1. Set up the environment 

SRA_downloader requires the conda package manager in order to work. I recommend using [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge), a conda installer that includes Mamba in the base environment. In addition, of course, it is essential to install Snakemake; following the steps in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 


### 2. Download SRA_downloader repository from GitHub.
Use git clone command to create a local copy. 

    git clone https://github.com/dcerdanv/sra_downloader.git


### 3. Configure the pipeline.

Before executing the pipeline, the users must configure it according to their samples. To do this, they must fill these files:

#### **a. config.yaml**

This is the pipeline configuration file, where you can tune all the available parameters to customise your downloads. It allows you to adjust the resources assigned to each rule, as well as add optional parameters to any of the tools used: prefetch, vdb-validate, fasterq-dump or pigz.

It also gives the user the possibility of decide whether or not they want to compress the FASTQ file or whether they want to delete the folder with the SRA object once the fasta has been dumped.

>**No default parameter** has been set for **fasterq-dump**. This would be equivalent to using fastq-dump with the parameters `--split-3` and `--skip-technical`. That is:
>```bash
>fastq-dump SRRXXXXXX --split-3 --skip-technical
>fasterq-dump SRRXXXXXX
>```
>You can find more information about fasterq-dump at [sra-tools documentation](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump), including the equivalences with fastq-dump.

Rename the example file ([`template_config.yaml`](https://github.com/dcerdanv/sra_downloader/blob/main/template_config.yaml)) to `config.yaml` and edit its contents to use it during your execution.


#### **b. samples.tsv**

This single column TSV is just a list with the accession numbers of the samples that want to be downloaded.

Rename the example file ([`template_samples.tsv`](https://github.com/dcerdanv/sra_downloader/blob/main/template_samples.tsv)) to `samples.tsv` and edit its contents to use it during your execution.


### 4. Create the Conda environments.

To run the pipeline, the user can create the conda environment first, which will take some minutes.
This step is done automatically using this command:
```bash
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba
```

If the environment has not been created before, it will be created before the first execution of the pipeline. 


### 5. Run the pipeline.

Once the pipeline is configured the user just needs to run SRA_downloader.
```bash
snakemake --use-conda -j 10
```

The mandatory arguments are:
* **--use-conda**: to install and use the conda environemnts.
* **-j**: number of jobs provided to snakemake.

## Pipeline steps

Here you can find a general description of the main steps of the pipeline.

1. **Prefetch sra (prefetch)**. Prefetch tool downloads SRA files and their dependencies.
    
        prefetch SRRXXXXXX

2. **Validate sra (vdb-validate)**. Vdb-validate validates the integrity of downloaded SRA data.

        vdb-validate SRRXXXXXX

    This step, like the one that removes the SRA object, will generate an empty file as a flag to confirm that the step has completed correctly. It is necessary for the pipe to works correctly, since this step does not modify or create any new file. This {sample}_validate.done file is stored in the 'flags' folder in the logs directory.

3. **Dump fastq (fasterq-dump)**. The fasterq-dump tool extracts data in FASTQ or FASTA format from SRA-accessions. It allows parallelization.

    This tool can return one or more fastq files, depending on the type of sequencing of the experiment (single-end, paired-end) or fasterq-dump's own configuration parameters. The compression tool, pigz, is prepared to deal with both situations 

        fasterq-dump -v -L info SRRXXXXXX

4. **Compress fastq (pigz)** *[Optional]*. Compress FASTQ sequence to gzip format. It can be enabled through the 'compress' parameter in `config.yaml`.

        pigz SRRXXXXXX.fastq

5. **Remove sra** *[Optional]*. Remove the SRA object folder. It can be enabled through the 'remove_sra' parameter in `config.yaml`.

        rm -r SRRXXXXXX

    This step, like the one that validates the SRA object, will generate an empty file as a flag to confirm that the step has completed correctly. It is necessary for the pipe to works correctly, since this step does not modify or create any new file. This {sample}_remove.done/{sample}_not_remove.done file is stored in the 'flags' folder in the logs directory.

## Cluster Usage

Although this pipeline can be launched on any machine, it is designed to be executed in environments with high computational capacity. If they use Slurm as a task manager, you can add the parameter `--slurm` which gives snakemake specific support for it. So the command would be something like: 

    snakemake --use-conda --slurm -j 32 


## Configuration of computation resources

The user can configure SRA_downloader to optimise the available computational resources in order to reduce the computational time. The optimization is achieved thanks to Snakemake's ability to parallelize the jobs and to assign specific resources to each rule in a simple way. Resources configuration is done through the configuration file (`config.yaml`). This file has a field called *resources*, where the user can define the RAM memory usage, the number of threads (if the rule admits multithreading) available to a specific step and the maximum execution time. Additionally, if the user does not provide any value for some of these fields, the pipeline will use the default values.
___



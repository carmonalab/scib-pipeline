# Pipeline for benchmarking supervised integration of single-cell RNA-seq atlases

This repository contains the snakemake pipeline for our benchmarking
analysis of scRNA-seq data integration tools. It is based on the package
[`scib`](https://github.com/theislab/scib.git) and the previous
[scib-pipeline](https://github.com/theislab/scib-pipeline) from the
Theis Lab [Luecken et al,
2020](https://doi.org/10.1038/s41592-021-01336-8).

Compared to this previous study our benchmark focus on (semi) supervised
tools with the addition of our new version of
[STACAS](https://github.com/carmonalab/STACAS). It includes our new
integration metric CiLISI we computed with our
[scIngrationMetrics](https://github.com/carmonalab/scIntegrationMetrics)
R package. It also assesses how well (semi)-supervised integration tools
are robust to noise we introduce with shuffling of cell type labels and
partial annotations. It also tests the capacity of (semi) supervised
tools to separate cell type when they are guided with a broader
annotation.

## Major modifications regarding the original pipeline

-   Adding STACAS and semi-supervised STACAS
-   Using embedding output (PCA computed on scaled integrated data with Seurat) for R based methods
-   Embedding/latent space for integration with a size fixed (e.g. 30, 50) for all tools (reduced space or bottleneck layer of autoencoders)
-   Testing noisy and missing cell type labels to guide integration (i.e. partially removing and shuffling cell type labels)
-   New batch-correction metric CiLISI (cell type-aware LISI) computed with [scIngrationMetrics](https://github.com/carmonalab/scIntegrationMetrics)
-   Packages of integration tools updated


## Installation

As in the original pipeline, to reproduce the results from this study,
two separate conda environments are needed for python and R operations.
Please make sure you have either
[`mambaforge`](https://github.com/conda-forge/miniforge) or
[`conda`](https://conda.io/projects/conda) installed on your system to
be able to use the pipeline. We recommend using
[`mamba`](https://mamba.readthedocs.io), which is also available for
conda, for faster package installations with a smaller memory footprint.

We provide python and R environment YAML files in `envs/`, together with
an installation script for setting up the correct environments in a
single command. based on the R version you want to use. Our new pipeline
currently only supports R 4.1 Call the script as follows

``` shell
bash envs/create_conda_environments.sh -r 4.1
```

Once installation is successful, you will have the python environment
`scib-pipeline-R<version>` and the R environment `scib-R<version>` that
you must specify in the [config file](#setup-configuration-file).

| R version | Python environment name | R environment name | Test data config YAML file   |
|------------------|------------------|------------------|-------------------|
| 4.1       | `scib-pipeline-R4.1`    | `scib-R4.1`        | `configs/test_data-R4.1.yml` |

## Running the Pipeline

This repository contains a
[snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to run
integration methods and metrics reproducibly for different data
scenarios preprocessing setups.

### Setup Configuration File {#setup-configuration-file}

The parameters and input files are specified in config files. A
description of the config formats and example files can found in
`configs/`. You can use the example config that use the test data to get
the pipeline running quickly, and then modify a copy of it to work with
your own data.

### Pipeline Commands

To call the pipeline on the test data e.g. using R 4.1 to reproduce our
benchmarking with the original annotations to guide supervised tools:

``` shell
snakemake --configfile configs/test_original_annotations-R4.1.yaml -n
```

This gives you an overview of the jobs that will be run. In order to
execute these jobs with up to 10 cores, call

``` shell
snakemake --configfile configs/test_original_annotations-R4.1.yaml --cores 10
```

We strongly recommand to use this snakemake on a HPC cluster e.g. using slurm 
and the config file `configs/cluster.yml` you can run the workflow as follow:

``` shell
mkdir -p cluster/snakemake/; \
snakemake -j 100 --configfile configs/test_original_annotations-R4.1.yaml \
--cluster-config configs/cluster.yml \
--cluster "sbatch -A {cluster.account} \
    -p {cluster.partition} \
    -N {cluster.N} \
    -t {cluster.time} \
    --job-name {cluster.name} \
    --mem {cluster.mem} \
    --cpus-per-task {cluster.cpus-per-task}\
    --output {cluster.output} \
    --error {cluster.error}"
```


Then you can generate a table gathering the snakemake benchmark files
(cpu time, memory usage...)

``` shell
snakemake --configfile configs/test_original_annotations-R4.1.yaml --cores 1 benchmarks
```

More snakemake commands can be found in the
[documentation](snakemake.readthedocs.io/).

## Reproduce/Visualize our results

We provide the config files to reproduce our 3 different analyses
together with the rmarkdown we used to generate the figures from the
results of the pipeline that you can find on the results directory

| Analysis             | config YAML file                                                                 | Rmarkdown file                                                                                   |
|------------------|-------------------------|-----------------------------|
| original annotations | [test_original_annotations-R4.1.yml](configs/test_original_annotations-R4.1.yml) | [originalAnnotationAnalysis.Rmd](results/original_annotations/originalAnnotationAnalysis.Rmd) |
| robustness to noise  | [test_supervised_methods-R4.1.yml](configs/test_supervised_methods-R4.1.yml)     | [SupervisedToolAnalysis.Rmd](results/supervised_tools_analysis/SupervisedToolAnalysis.Rmd)       |
| final benchmark      | [test_final_benchmark-R4.1.yml](configs/test_final_benchmark-R4.1.yml)           | [finalBenchmarkAnalysis.Rmd](results/final_benchmark/finalBenchmarkAnalysis.Rmd)                 |



## Tools

Tools that are compared include: -
[STACAS](https://github.com/carmonalab/STACAS) -
[Scanorama](https://github.com/brianhie/scanorama) - -
[scANVI](https://github.com/chenlingantelope/HarmonizationSCANVI) - -
[FastMNN](https://bioconductor.org/packages/batchelor/) - -
[scGen](https://github.com/theislab/scgen) - -
[scVI](https://github.com/YosefLab/scVI) - [Seurat v4 (CCA and
RPCA)](https://github.com/satijalab/seurat) - -
[Harmony](https://github.com/immunogenomics/harmony)


## References

[Benchmarking atlas-level data integration in single-cell genomics. Luecken et al, 2020](https://doi.org/10.1038/s41592-021-01336-8)

[Semi-supervised integration of single-cell transcriptomics data. Andreatta et al, 2023](https://doi.org/10.1101/2023.07.07.548105)

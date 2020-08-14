# cofrag: cfDNA co-fragmentation analysis

## Getting Started

* Input and output: cofrag heavily relies on UNIX pipes as IO.
* Coordinates: unless otherwise stated, cofrag IO assumes UCSC BED convention, i.e. 0-based half-open coordinates.

Examples:

Read fragmentation data from stdin and calculate the contact matrix (the output is a genome contact matrix in BEDPE format):

```
gzip -dc frag.tsv.gz | ./cofrag.R --verbose contact > contact_matrix.bed
```

Perform the analysis using multiple (four) CPU cores:

```
gzip -dc frag.tsv.gz | ./cofrag.R -v contact -n 4 > contact_matrix.bed
```

Only use fragments no longer than 500bp and with MAPQ >= 30:

```
gzip -dc frag.tsv.gz | ./cofrag.R contact -q 30 --max-frag-size 500 > contact_matrix.bed
```

Specify the bin size (500k bp):

```
gzip -dc frag.tsv.gz | ./cofrag.R contact --bin-size 500k > contact_matrix.bed
```

Process a subset of the fragment file (need tabix index file):

```
tabix frag.tsv.gz 22:16000001-28000000 | ./cofrag.R contact > contact_matrix.bed
```

Call A/B compartments from a contact matrix file:

```
cat contact_matrix.bed | ./cofrag.R compartment --genes gene_density.bed
```

## Prerequisites

We sugguest using conda to mananage the environment:

```
conda config --add channels bioconda && conda config --add channels conda-forge
conda install -y r-base=3.6.3
conda install -y r-essentials=3.6
conda install -y r-optparse=1.6.6
conda install -y r-here=0.1
conda install -y r-doparallel=1.0.15
conda install -y r-logging=0.10
conda install -y -c bioconda bioconductor-genomicranges=1.38.0
```

## Containerization

We have created a Docker image for cofrag (https://hub.docker.com/repository/docker/epifluidlab/cofrag):

```
docker pull epifluidlab/cofrag:v0.1.1
```

You can find the `Dockerfile` at: https://github.com/epifluidlab/cofrag/blob/master/Dockerfile

Apart from Docker, we also encourage you to run the container through singularity, since it provides better support for pipes, file system binding, and doesn't require root to run:

```
singularity pull docker://epifluidlab/cofrag:v0.1.1
```

An example of running cofrag through singularity:

```
gzip -dc frag.tsv.gz | singularity exec PROJECT_PATH/cofrag.R contact > contact_matrix.bed
```

Since cofrag analysis can run across multiple CPU cores, prior to running the container, you may want to check the Docker/Singularity CPU resource management strategies.

## Usage

Main command arguments:

* `-v`, `--verbose`: verbose logging

Sub commands: `contact` and `compartment`:

### Contact matrix

Input: read fragment data from stdin. cofrag assume the data is organized in a BED-like format:

```
22	16050006	16050156	27	-
22	16050029	16050203	9	-
22	16050032	16050200	21	-
22	16050036	16050203	0	-
```

Meanings of the five columns:

1. chromosome name
2. fragment start (inclusive)
3. fragment end (exclusive)
4. MAPQ
5. strand

Output: contact matrix in BEDPE format.

Arguments:

* `-s`, `--bin-size`: specify the size of each bin. Can be a plain integer, or notations such as 200k or 1m. Default: `500k`.
* `-b`, `--block-size`: in order to facilitate multi-thread processing, cofrag comes with a schedule which first divides the genome range into multiple *blocks*, then creates a fleet of workers and assigns blocks to each worker. After all workers have done their jobs, all results are merged together. This argument specifies how large the block is. The size should be multiples of `bin_size`. Default: `5m`.
* `-m`, `--metrics`: specify which statistic test to use for calculating the distance. `ks` for Kolmogorov–Smirnov test, `cucconi` for Cucconi test. Default: `ks`
* `-q`, `--min-mapq`: exclude all fragments with lower MAPQ scores. Default: `0`.
* `--max-frag-size`: exclude exceptionally long fragments. Default: `1000`.
* `-n`, `--num-cores`: number of CPU cores. Default: `1`.

### Calling A/B compartment

Input:

* Read contact matrix from stdin. The data should be organized in BEDPE format, for example:

```
22      16000000        16500000        22      16000000        16500000        15.698970004336
22      16000000        16500000        22      16500000        17000000        12.8406235365063
22      16500000        17000000        22      16500000        17000000        15.698970004336
22      16000000        16500000        22      17000000        17500000        0
```

* A gene density file (BED format). This will be used to determine whether the compartment is open or closed. For example:

```
22      16062157        16063236        1
22      16076052        16076172        1
22      16084249        16084826        1
22      16100517        16124973        1
```

Output: compartment data in BED format. For example:

```
22	16000001	16500001	0.130025497104604
22	16500001	17000001	0.131736981249446
22	17000001	17500001	-0.127187986290912
22	17500001	18000001	-0.0618555647760052
```

Arguments:

* `-g`, `--genes`: the gene density file. Can be plain text, .gz or .bzip2.

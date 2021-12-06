![Build Status](https://github.com/iqbal-lab-org/viridian_workflow/actions/workflows/build.yaml/badge.svg)

# Viridian Workflow

## Installation

Clone this repository and then from its root, build either a Singularity
or Docker container.

To build a singularity container:
```
singularity build viridian_workflow.img Singularity.def
```

To build a docker container:
```
docker build --network=host .
```
(without `--network=host` you will likely get `pip install` timing out and
the build failing).

Both the Docker and Singularity container will have the main script
`viridian_workflow` installed.

## Usage

To run on paired Illumina reads:
```
viridian_workflow run_one_sample \
  --tech illumina
  --ref_fasta data/MN908947.fasta \
  --amplicon_json data/covid-artic-v3.json \
  --reads1 reads_1.fastq.gz \
  --reads2 reads_2.fastq.gz \
  --outdir OUT
```
To run on unpaired nanopore reads:
```
viridian_workflow run_one_sample \
  --tech ont
  --ref_fasta data/MN908947.fasta \
  --amplicon_json data/covid-artic-v3.json \
  --reads reads.fastq.gz \
   --outdir OUT
```

The FASTA and JSON files in those commands can be found in the `data/`
directory of this repository.

Other options:
* `--sample_name MY_NAME`: use this to change the sample name
  (default is "sample") that is put in the final FASTA file, BAM file, and
  VCF file.
* `--keep_bam`: use this option to keep the BAM file of original input reads
  mapped to the reference genome.
* `--force`: use with caution - it will overwrite the output directory if
  it already exists.



## Output files

The default files in the output directory are:

* `consensus.fa`: a FASTA file of the consensus sequence.
* `variants.vcf`: a VCF file of the identified variants between the consensus
  sequence and the reference genome.
* `log.json`: contains logging information for the viridian workflow run.

If the option `--keep_bam` is used, then a sorted BAM file of the reads mapped
to the reference will also be present, called
`reference_mapped.bam` (and its index file `reference_mapped.bam.bai`).

## Configuration

Quality control thresholds are configured in `viridian_workflow/config.ini`:

```INI
[minimap2]
threads = 1

[viridian]
illumina_read_lengths = 50

[qcovid]
min_depth = 50
min_template_coverage_75 = 0.80
freq_threshold = 0.80
```

Viridian workflow applies quality control metrics at the amplicon and base level, as follows:

* Reads are aligned to the reference (default `MN908947`). This removes reads that do not sufficiently match this sequence.
* Individual reads must be sufficiently long. For Illumina reads, `illumina_read_lengths` must be, by default, at least 50.
* Templates are reconstructed from either paired or single end reads that have been aligned. The start and end positions are used to infer which amplicon they belong to (e.g. `nCoV-artic-v3.bed`). "off-target" templates are discarded.
* Enough templates must cover an amplicon or else all calls inside the region covered solely by that amplicon are voided with `N`s: `min_depth`, default 50.
* Enough templates must span at least 75% of the amplicon they belong to, excluding internal indels, where template length is the `(alignment_end - alignment_start) / amplicon_length`. This threshold is controlled by `min_template_coverage_75`, default 80%. This filters out short fragments caused by PCR artefacts. If this threshold fails the entire amplicon is discarded.
* After the assembly is constructed the original reads are mapped to it. Individual bases are voided with `N` if there is less than `freq_threshold` agreement. Default is 80%.

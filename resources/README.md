# Reference files

The AMBER processing scripts need a genome FASTA with a bwa-mem2 index and a RepeatMasker bed. These are large or genome-specific, so they are built on your machine rather than stored here. `2.mapUG.sh` reads their locations from environment variables.

| Variable | File | Used for |
|---|---|---|
| `HG38_FASTA` / `MM10_FASTA` / `HGMM_FASTA` | Genome FASTA (bwa-mem2 index files next to it) | Alignment |
| `HG38_RMSK` / `MM10_RMSK` / `HGMM_RMSK` | RepeatMasker bed, milliDiv < 140 | Removing reads in repeats |

`HGMM_*` is only needed for the human+mouse species-mixing (barnyard) experiment.

## 1. Genome FASTA and bwa-mem2 index

Download the assemblies from UCSC (use the same one throughout):

```
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip hg38.fa.gz mm10.fa.gz
```

Build the bwa-mem2 index (needs roughly 64 GB RAM for a human or mouse genome, so run it on a compute node):

```
bwa-mem2 index hg38.fa
bwa-mem2 index mm10.fa
```

The scripts keep standard chromosomes only (`chr1`..`chrN`, `chrX`, `chrY`) and drop alternate contigs and `chrM`.

### Human + mouse combined genome (species mixing only)

Prefix the chromosome names with `hg_` and `mm_`, concatenate, and index:

```
sed 's/^>/>hg_/' hg38.fa >  hg38_mm10_combined.fa
sed 's/^>/>mm_/' mm10.fa >> hg38_mm10_combined.fa
bwa-mem2 index hg38_mm10_combined.fa
```

## 2. RepeatMasker bed

```
bash build_repeatmasker_bed.sh hg38 hg38.repeatmasker.milliDiv_lt140.bed
bash build_repeatmasker_bed.sh mm10 mm10.repeatmasker.milliDiv_lt140.bed
```

This downloads the UCSC `rmsk` table and keeps repeats with `milliDiv < 140` (the combined human+mouse bed command is in the header of the script).

## 3. Point the pipeline at them

```
export HG38_FASTA=/path/to/hg38.fa
export HG38_RMSK=/path/to/hg38.repeatmasker.milliDiv_lt140.bed
export MM10_FASTA=/path/to/mm10.fa
export MM10_RMSK=/path/to/mm10.repeatmasker.milliDiv_lt140.bed
```

Run these in the shell before `sbatch`, and the job inherits them.

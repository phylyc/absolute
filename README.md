# ABSOLUTE quantification of somatic DNA alterations in cancer

## Installation

Create a virtual environment and use the `install.sh` script to install all dependencies, or use [this docker](https://hub.docker.com/repository/docker/phylyc/absolute/general).

## Usage

1. Run library/scripts/run_absolute.R on your segmented copy ratio (and SNV/INDEL) data. 
2. Inspect the generated output plots and write down the best-matching solution.
3. Run library/scripts/extract_solution.R on the output from the first step and specify the solution number.


### Input for library/scripts/run_absolute.R

#### Required

```--sample```: Sample name

```--seg_dat_fn```:
File path to segmentation table, which is assumed to be in AllelicCapSeg format (column order does not matter):
```
Chromosome	Start.bp	End.bp	n_probes	length	n_hets	f	tau	sigma.tau	mu.minor	sigma.minor	mu.major	sigma.major	SegLabelCNLOH
```

- Columns:
  - `n_probes`: Number of target intervals in each segment
  - `n_hets`: Number of heterozygous SNPs in each segment
  - `f`: minor allele fraction
  - `tau`: total copy number
  - `mu.minor`: minor allelic copy number (f * tau)
  - `mu.major`: major allelic copy number ((1 - f) * tau)
  - `sigma.xxx`: standard error of `xxx`
  - `SegLabelCNLOH`: copy-neutral loss of heterozygosity label: 0 is flanked on both sides, 1 is one side, 2 is no cn.loh

See [this script](https://github.com/phylyc/somatic_workflow/blob/master/python/acs_conversion.py) on how to convert GATK ModelSegments output to that format.

#### Optional

```--results_dir```: Directory into which the result files will be written

```--ssnv_skew```: ~ 2 / (1 + ref_bias), where the reference bias skews the observed distribution of alternate allele read counts towards f / (f + (1 - f) * ref_bias) for true minor allele fraction f. 

``--maf``: somatic SNV table with functional annotations (e.g. output of GATK Oncotator or Funcotator). Note that no entry in the table cells can be longer than ~1000 characters.

Alternatively, ``--maf`` accepts a **VCF file** (`.vcf`, `.vcf.gz`, or `.vcf.bgz`). The VCF is parsed directly into the internal mutation table, so a separate ``--indel_maf`` is not needed — SNVs and indels are detected automatically from the REF/ALT alleles. Requirements:
- A `#CHROM` header line and a per-sample `FORMAT` genotype column.
- Tumor allelic depths read from `FORMAT/AD` (`ref,alt`), with `RO`+`AO` (freebayes) and `AF`+`DP` as fallbacks.
- Only `PASS` / unfiltered (`.`) records are kept.
- For multi-sample VCFs (e.g. tumor/normal), the tumor genotype column is selected from the `##tumor_sample=` header (Mutect2) if present, otherwise from the ``--sample`` value; a single-sample VCF is used as-is.
- Chromosome names must match the segmentation table (e.g. `1` vs `chr1`), otherwise mutations are silently dropped during segment mapping.

Functional annotations (gene, variant classification, protein change, COSMIC counts) are not available from a bare VCF, so the affected report fields are left blank; the purity/ploidy/CCF fit is unaffected.

```--indel_maf```: somatic INDEL table with functional annotations (output of GATK Oncotator or Funcotator). Ignored when ``--maf`` is a VCF.

```--gender```: biological sex for ploidy assumptions; {"F", "Female", "female", "XX", "M", "Male", "male", "XY"}

```--alpha```: purity for force-calling

```--tau```: ploidy for force-calling

```--refit_alpha_tau```: flag (off by default). When set together with ``--alpha`` and ``--tau``, the given purity/ploidy are treated as *seeds* rather than taken at face value: ABSOLUTE locally optimizes alpha/tau in a small window around the seeds to snap them onto the nearest well-fitting mode. Useful when the force-called values are approximate (e.g. carried over from another tool or a previous run). Without this flag, ``--alpha``/``--tau`` are used exactly as provided.

```--refit_window_alpha```: half-width of the purity search window used when ``--refit_alpha_tau`` is set (default 0.05, i.e. alpha ± 0.05).

```--refit_window_tau```: half-width of the ploidy search window used when ``--refit_alpha_tau`` is set. If unset, it is auto-derived from ``--refit_window_alpha`` as ``(2 / alpha^2) * refit_window_alpha`` (capped at 1.0) to match the copy-ratio comb spacing.

```--copy_num_type```: {"allelic", "total"}; determines purity/ploidy based on allelic or total copy ratios. Both modes are supported. This value is also embedded into the output file names (see below) so allelic and total runs of the same sample do not overwrite each other.

```--genome_build```: This package currently supports human {hg18, hg19, hg38} and mouse {mm9, mm10, mm39}

```--pkg_dir```: path to folder in which the library folder lies.



### Output of library/scripts/run_absolute.R

1. A PDF containing plots to choose the best fitting solution from.
2. An RData object containing all input and model data.


### Input for library/scripts/extract_solution.R

#### Required

```--sample```: Sample name

```--rdata```: RData object from first step.

```--solution_num```: Ordinal number of the picked solution.

```--results_dir```: Directory into which a `results` folder will be created

```--analyst_id```: Initials of the analyst who picked the solution (for future blame :P)

#### Optional

```--copy_num_type```: {"allelic", "total"}; determines purity/ploidy based on allelic or total copy ratios. Both modes are supported.

```--genome_build```: This package currently supports human {hg18, hg19, hg38} and mouse {mm9, mm10, mm39}

```--pkg_dir```: path to folder in which the library folder lies.


### Output of library/scripts/extract_solution.R

The `copy_num_type` (`allelic` or `total`) is embedded into each output file name so that allelic and total runs of the same sample can coexist in the same results directory.

1. `.../reviewed/SEG_MAF/*.ABS_MAF.<copy_num_type>.txt`: annotated SNV/INDEL MAF file.
2. `.../reviewed/SEG_MAF/*.segtab.<copy_num_type>.txt`: Absolute copy number estimates for the segmentation table
3. `.../reviewed/SEG_MAF/*.IGV.seg.<copy_num_type>.txt`: IGV-compatible segmentation table of absolute copy number
4. `.../reviewed/*.<analyst_id>.ABSOLUTE.table.<copy_num_type>.txt`: purity/ploidy table

Consider running [this script](https://github.com/phylyc/somatic_workflow/blob/master/python/map_to_absolute_copy_number.py) on the output to rescue any dropped segments.


## Cite

If you use this package, please consider citing: 

- Carter, S., Cibulskis, K., Helman, E. et al. Absolute quantification of somatic DNA alterations in human cancer. Nat Biotechnol 30, 413–421 (2012). https://doi.org/10.1038/nbt.2203

If you used any of the suggested pre- or post-processing scripts, please consider citing:

- https://github.com/phylyc/somatic_workflow/
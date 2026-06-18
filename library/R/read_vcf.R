## --------------------------------------------------------------------------
## VCF input support for ABSOLUTE
##
## ABSOLUTE natively consumes MAF files to fit SSNV modes. This module lets a
## user pass a VCF (.vcf / .vcf.gz / .vcf.bgz) to the --maf input instead. The
## VCF is parsed into MAF-shaped data frames that classic_CreateMutCnDat()
## already understands, so nothing downstream has to change.
##
## Only four pieces of information are actually required to fit SSNVs:
##   Chromosome, position, tumor alt read depth, tumor ref read depth.
## These come from CHROM, POS and the per-sample AD field (with RO/AO and
## AF+DP fallbacks). Everything else a MAF carries (Hugo_Symbol,
## Variant_Classification, Protein_Change, dbSNP_Val_Status, ...) is only used
## for annotation / reporting and is filled with neutral placeholders here.
##
## SNVs and indels are separated by REF/ALT allele length and returned as two
## frames. classic_CreateMutCnDat() applies indel-specific filtering to the
## indel frame; the beta-binomial error model itself treats both identically
## (single global SSNV_skew / rho / epsilon).
## --------------------------------------------------------------------------

is_vcf_file <- function(fn) {
  if (is.null(fn) || length(fn) == 0 || is.na(fn)) return(FALSE)
  grepl("\\.vcf(\\.gz|\\.bgz)?$", fn, ignore.case = TRUE)
}

## MAF columns produced for the ABSOLUTE pipeline. Start_Position / End_Position
## use the standard MAF capitalization so rename_maf_colnames() converts them.
.vcf_maf_template <- function() {
  data.frame(
    Hugo_Symbol            = character(0),
    Chromosome             = character(0),
    Start_Position         = integer(0),
    End_Position           = integer(0),
    Variant_Classification = character(0),
    Variant_Type           = character(0),
    Reference_Allele       = character(0),
    Tumor_Seq_Allele2      = character(0),
    Tumor_Sample_Barcode   = character(0),
    Protein_Change         = character(0),
    UniProt_AApos          = character(0),
    t_ref_count            = integer(0),
    t_alt_count            = integer(0),
    stringsAsFactors       = FALSE
  )
}

## SNP / DNP / ... vs INS / DEL purely from allele lengths.
classify_variant <- function(ref, alt) {
  rl <- nchar(ref); al <- nchar(alt)
  if (rl == 1L && al == 1L) return(list(type = "SNP", is.indel = FALSE))
  if (rl == al)             return(list(type = if (rl == 2L) "DNP" else if (rl == 3L) "TNP" else "ONP",
                                        is.indel = FALSE))
  if (rl < al)              return(list(type = "INS", is.indel = TRUE))
  return(list(type = "DEL", is.indel = TRUE))
}

## Extract tumor ref/alt read depths from one FORMAT + sample genotype pair.
## Returns list(ref = <int>, alt = <int vector, one per ALT allele>) or NULL.
parse_vcf_counts <- function(format.str, sample.str) {
  if (is.na(format.str) || is.na(sample.str)) return(NULL)
  fkeys <- strsplit(format.str, ":", fixed = TRUE)[[1]]
  svals <- strsplit(sample.str, ":", fixed = TRUE)[[1]]
  get <- function(key) {
    i <- match(key, fkeys)
    if (is.na(i) || i > length(svals)) return(NA_character_)
    svals[i]
  }

  ## Preferred: AD = "ref,alt1,alt2,..." (GATK/Mutect2, etc.)
  ad <- get("AD")
  if (!is.na(ad) && ad != ".") {
    v <- suppressWarnings(as.integer(strsplit(ad, ",", fixed = TRUE)[[1]]))
    if (length(v) >= 2L && !is.na(v[1])) return(list(ref = v[1], alt = v[-1]))
  }

  ## freebayes: RO (ref obs) + AO (alt obs, one per ALT)
  ro <- get("RO"); ao <- get("AO")
  if (!is.na(ro) && !is.na(ao) && ro != "." && ao != ".") {
    return(list(ref = suppressWarnings(as.integer(ro)),
                alt = suppressWarnings(as.integer(strsplit(ao, ",", fixed = TRUE)[[1]]))))
  }

  ## Last resort: reconstruct from AF + DP.
  af <- get("AF"); dp <- get("DP")
  if (!is.na(af) && !is.na(dp) && af != "." && dp != ".") {
    dpi <- suppressWarnings(as.integer(dp))
    afs <- suppressWarnings(as.numeric(strsplit(af, ",", fixed = TRUE)[[1]]))
    if (!is.na(dpi) && !any(is.na(afs))) {
      alt <- as.integer(round(afs * dpi))
      return(list(ref = as.integer(dpi - sum(alt)), alt = alt))
    }
  }

  return(NULL)
}

## Decide which genotype column is the tumor. Honors an explicit name, then a
## Mutect2 "##tumor_sample=" header, then a sole sample; otherwise errors.
pick_tumor_sample <- function(meta.lines, sample.cols, tumor.sample = NA) {
  if (!is.na(tumor.sample)) {
    if (!tumor.sample %in% sample.cols) {
      stop("Requested tumor sample '", tumor.sample, "' not found in VCF samples: ",
           paste(sample.cols, collapse = ", "))
    }
    return(tumor.sample)
  }
  ts <- grep("^##tumor_sample=", meta.lines, value = TRUE)
  if (length(ts) >= 1L) {
    name <- sub("^##tumor_sample=", "", ts[length(ts)])
    if (name %in% sample.cols) return(name)
  }
  if (length(sample.cols) == 1L) return(sample.cols[1])
  stop("VCF has multiple samples (", paste(sample.cols, collapse = ", "),
       ") and no usable ##tumor_sample header. ",
       "Specify the tumor sample via --sample.")
}

read_vcf_as_mafs <- function(vcf.fn, tumor.sample = NA, pass.only = TRUE, verbose = FALSE) {
  con <- if (grepl("\\.(gz|bgz)$", vcf.fn, ignore.case = TRUE)) gzfile(vcf.fn, "rt") else file(vcf.fn, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)

  chrom.ix <- grep("^#CHROM", lines)
  if (length(chrom.ix) == 0L) stop("VCF has no #CHROM header line: ", vcf.fn)
  chrom.ix <- chrom.ix[1]
  meta.lines <- lines[grep("^##", lines)]
  col.names  <- strsplit(sub("^#", "", lines[chrom.ix]), "\t", fixed = TRUE)[[1]]

  if (!"FORMAT" %in% col.names) {
    stop("VCF lacks a FORMAT column; per-sample allelic depths (AD/RO+AO/AF+DP) are required.")
  }
  fmt.pos <- match("FORMAT", col.names)
  sample.cols <- col.names[(fmt.pos + 1L):length(col.names)]
  if (length(sample.cols) == 0L) stop("VCF has no sample / genotype columns.")
  tumor.sample <- pick_tumor_sample(meta.lines, sample.cols, tumor.sample)
  if (verbose) print(paste("VCF: using tumor sample '", tumor.sample, "'", sep = ""))

  body.lines <- lines[(chrom.ix + 1L):length(lines)]
  body.lines <- body.lines[nchar(body.lines) > 0L]
  if (length(body.lines) == 0L) {
    if (verbose) print("VCF contains no variant records.")
    empty <- .vcf_maf_template()
    return(list(snv = empty, indel = empty))
  }

  fields <- strsplit(body.lines, "\t", fixed = TRUE)
  n.expected <- length(col.names)
  bad <- vapply(fields, length, integer(1)) != n.expected
  if (any(bad)) {
    print(paste("VCF: dropping", sum(bad), "malformed records with unexpected column count"))
    fields <- fields[!bad]
  }
  vcf <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
  colnames(vcf) <- col.names

  if (pass.only && "FILTER" %in% col.names) {
    keep <- vcf[["FILTER"]] %in% c("PASS", ".", "")
    if (verbose) print(paste("VCF: keeping", sum(keep), "of", length(keep), "PASS/unfiltered records"))
    vcf <- vcf[keep, , drop = FALSE]
  }
  if (nrow(vcf) == 0L) {
    empty <- .vcf_maf_template()
    return(list(snv = empty, indel = empty))
  }

  alt.list <- strsplit(vcf[["ALT"]], ",", fixed = TRUE)
  is.symbolic <- function(a) a %in% c("*", ".", "<NON_REF>") | grepl("^<", a)
  total <- sum(vapply(alt.list, function(a) sum(!is.symbolic(a)), integer(1)))

  chrom <- character(total); start <- integer(total); endp <- integer(total)
  refA  <- character(total); altA  <- character(total)
  vtype <- character(total); isindel <- logical(total)
  rc <- integer(total); ac <- integer(total)

  fmt.col <- vcf[["FORMAT"]]
  smp.col <- vcf[[tumor.sample]]
  k <- 0L
  for (i in seq_len(nrow(vcf))) {
    alts   <- alt.list[[i]]
    counts <- parse_vcf_counts(fmt.col[i], smp.col[i])
    pos    <- suppressWarnings(as.integer(vcf[["POS"]][i]))
    refal  <- vcf[["REF"]][i]
    for (j in seq_along(alts)) {   ## j indexes ALT alleles == AD alt index
      a <- alts[j]
      if (is.symbolic(a)) next
      k <- k + 1L
      cls <- classify_variant(refal, a)
      chrom[k]   <- vcf[["CHROM"]][i]
      start[k]   <- pos
      endp[k]    <- pos + max(0L, nchar(refal) - 1L)
      refA[k]    <- refal
      altA[k]    <- a
      vtype[k]   <- cls$type
      isindel[k] <- cls$is.indel
      rc[k] <- if (!is.null(counts)) counts$ref else NA_integer_
      ac[k] <- if (!is.null(counts) && j <= length(counts$alt)) counts$alt[j] else NA_integer_
    }
  }

  maf <- data.frame(
    Hugo_Symbol            = "Unknown",
    Chromosome             = chrom,
    Start_Position         = start,
    End_Position           = endp,
    Variant_Classification = "Unknown",
    Variant_Type           = vtype,
    Reference_Allele       = refA,
    Tumor_Seq_Allele2      = altA,
    Tumor_Sample_Barcode   = tumor.sample,
    Protein_Change         = NA_character_,
    UniProt_AApos          = NA_character_,
    t_ref_count            = rc,
    t_alt_count            = ac,
    stringsAsFactors       = FALSE
  )

  snv   <- maf[!isindel, , drop = FALSE]
  indel <- maf[ isindel, , drop = FALSE]
  if (verbose) {
    print(paste("VCF: parsed", nrow(snv), "SNV/MNP and", nrow(indel), "indel records",
                "across chromosomes:", paste(utils::head(unique(chrom), 5), collapse = ", "),
                if (length(unique(chrom)) > 5) "..." else ""))
    print(paste("NOTE: Chromosome names must match the segmentation data",
                "(e.g. '1' vs 'chr1'); mismatched names silently drop mutations."))
  }
  list(snv = snv, indel = indel)
}

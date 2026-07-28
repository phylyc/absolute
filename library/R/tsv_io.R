## Robust tab-delimited input / output for MAF and seg tables.
##
## Two read.delim() defaults are wrong for MAF and both fail silently:
##
##  * quote="\"" -- MAF has no quoting convention, so a double quote is a literal
##    data character. Funcotator emits HGNC alias / previous-name values that are
##    quoted comma-separated lists (`"plexin 2", "plexin-A2"`). An odd number of
##    quotes on a line opens a quoted region that swallows every following tab
##    *and newline* until the next quote, merging records into one giant field and
##    dropping every absorbed row. Seen in production: a truncated HGNC alias
##    `"transforming growth factor-&` (the `&beta;` entity was cut at the `;`) in
##    one SSNV MAF line reduced 454 parsed mutations to 271 with nothing but an
##    "EOF within quoted string" warning, and the surviving record carried 331 kB
##    of raw file text in HGNC_Alias_names -- which write.table(quote=FALSE) then
##    emitted verbatim, producing a non-rectangular ABS_MAF.
##
##  * comment.char="#" -- R truncates a line at a "#" found *anywhere* in it, not
##    just at the start, silently NA-filling the rest of the row. MAF uses "#" as
##    a leading header-block marker but also carries it inside free-text
##    annotation, so the leading block has to be counted and skipped explicitly.
##
## Reading is therefore always quote-blind and comment-blind, with the leading
## comment block located by an explicit pre-scan, and the parsed row count checked
## against the number of data lines on disk so this class of corruption can never
## again pass as a successful read.


.tsv_connection = function(fn) {
  if (grepl("\\.(gz|bgz)$", fn, ignore.case = TRUE)) gzfile(fn, "rt") else file(fn, "rt")
}


## One pass over the file, reporting everything the reader needs to know:
##   skip      - size of the leading comment/blank block, for read.delim(skip=)
##   n_data    - non-blank lines after that block, minus the header line
##   bad_quote - line numbers whose double quotes are unbalanced; harmless to us
##               now, but they break any quote-aware consumer downstream, so they
##               are worth naming rather than swallowing.
.scan_tsv_layout = function(fn, comment = "#", header = TRUE, max_report = 5L,
                            chunk_size = 50000L) {
  con = .tsv_connection(fn)
  on.exit(close(con))

  skip = 0L
  in_lead = TRUE
  n_content = 0L
  bad_quote = integer(0)
  offset = 0L

  repeat {
    chunk = readLines(con, n = chunk_size, warn = FALSE)
    if (length(chunk) == 0L) break

    ## nzchar() and useBytes=TRUE keep this safe on lines that are not valid in
    ## the current locale's encoding, which annotation columns occasionally are.
    blank = !nzchar(chunk)
    keep = !blank

    if (in_lead) {
      skippable = blank | (nzchar(comment) & startsWith(chunk, comment))
      first_content = which(!skippable)[1]
      if (is.na(first_content)) {
        skip = skip + length(chunk)
        offset = offset + length(chunk)
        next
      }
      skip = skip + (first_content - 1L)
      keep[seq_len(first_content - 1L)] = FALSE
      in_lead = FALSE
    }

    n_content = n_content + sum(keep)

    if (length(bad_quote) < max_report) {
      n_quotes = nchar(chunk, type = "bytes") -
                 nchar(gsub('"', "", chunk, fixed = TRUE, useBytes = TRUE), type = "bytes")
      hit = which(keep & (n_quotes %% 2L == 1L))
      if (length(hit) > 0L) {
        bad_quote = head(c(bad_quote, offset + hit), max_report)
      }
    }

    offset = offset + length(chunk)
  }

  list(skip = skip,
       n_data = max(n_content - as.integer(isTRUE(header)), 0L),
       bad_quote = bad_quote)
}


## read.delim() for MAF / seg tables: quote-blind, comment-blind, and checked in the
## two ways this format gets silently mangled.
.read_tsv = function(fn, comment = "#", header = TRUE, check_rows = TRUE,
                     label = basename(fn), ...) {
  layout = .scan_tsv_layout(fn, comment = comment, header = header)

  if (length(layout$bad_quote) > 0L) {
    print(paste("WARNING: ", label, ": unbalanced double quote on line(s) ",
                paste(layout$bad_quote, collapse = ", "),
                " (truncated annotation value?). Parsed here regardless, but any",
                " quote-aware reader will misparse this file.", sep = ""))
  }

  ## fill=FALSE, unlike the read.delim() default: a short line means the writer put a
  ## tab or a newline inside a field, and filling it with NA turns that corruption
  ## into a plausible-looking row. R's own error names the offending line.
  dat = read.delim(fn, header = header, skip = layout$skip,
                   quote = "", comment.char = "", fill = FALSE,
                   row.names = NULL, stringsAsFactors = FALSE,
                   check.names = FALSE, blank.lines.skip = TRUE, ...)

  ## An explicit nrows= is a deliberate partial read, not a corrupt one.
  if ("nrows" %in% names(list(...))) { check_rows = FALSE }

  ## Catches the complementary failure: rows that parsed at full width but got
  ## absorbed into a neighbour, which leaves the count short with no ragged line.
  if (check_rows && nrow(dat) != layout$n_data) {
    stop(label, ": parsed ", nrow(dat), " rows but the file holds ",
         layout$n_data, " data lines. Refusing to continue on a partial read.")
  }

  return(dat)
}


## MAF is written unquoted, so a tab or newline inside a field would silently
## destroy the row structure. Replace them and name the columns, so a bad
## annotation source gets reported instead of producing a corrupt output file.
.sanitize_tsv_fields = function(dat, label = "table") {
  dat = as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  offenders = character(0)

  for (j in seq_along(dat)) {
    col = dat[[j]]
    if (!is.character(col) && !is.factor(col)) { next }
    col = as.character(col)
    hit = !is.na(col) & grepl("[\t\r\n]", col, useBytes = TRUE)
    if (!any(hit)) { next }

    col[hit] = gsub("[\t\r\n]+", " ", col[hit], useBytes = TRUE)
    dat[[j]] = col
    offenders = c(offenders, paste(colnames(dat)[j], "=", sum(hit), sep = ""))
  }

  if (length(offenders) > 0L) {
    print(paste("WARNING: ", label, ": replaced embedded tab/newline characters",
                " before writing: ", paste(offenders, collapse = ", "), sep = ""))
  }

  return(dat)
}


## write.table() for MAF / seg tables: plain, unquoted, guaranteed rectangular.
.write_tsv = function(dat, out_fn, label = basename(out_fn), row.names = FALSE,
                      col.names = TRUE, append = FALSE) {
  dat = .sanitize_tsv_fields(dat, label = label)
  write.table(dat, file = out_fn, sep = "\t", quote = FALSE,
              row.names = row.names, col.names = col.names, append = append)
}

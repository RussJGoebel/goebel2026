#' Recode a land-cover code vector to names
#'
#' @param x Integer/character vector of land-cover codes (e.g., 1..8).
#' @param key Named character vector mapping codes -> labels
#'   (e.g., c("1"="urban", "2"="croplands", ...)).
#' @param as_factor Logical; return a factor with levels in key order? Default TRUE.
#'
#' @return A character (or factor) vector of labels; unknown codes become NA.
#' @examples
#' landcover_key <- c(
#'   "1"="urban","2"="croplands","3"="forest","4"="shrub & scrub",
#'   "5"="grass","6"="water","7"="wetlands","8"="bare surface"
#' )
#' x <- c(1,3,6,NA,9)
#' recode_landcover_codes(x, landcover_key)
#' @export
recode_landcover_codes <- function(x, key, as_factor = TRUE) {
  out <- unname(key[as.character(x)])
  if (as_factor) {
    out <- factor(out, levels = unname(key))
  }
  out
}


#' Rename proportion_* columns from numeric suffixes to land-cover names and fill NAs with 0
#'
#' @param df A data.frame/tibble with wide columns like `proportion_1`, ..., `proportion_8`.
#' @param key Named character vector mapping codes -> labels (names must be "1","2",...).
#' @param prefix String prefix used by the wide columns (default "proportion_").
#' @param keep_prefix If TRUE (default), keep the prefix in the new names
#'   (e.g., "proportion_urban"); if FALSE, use just the label (e.g., "urban").
#'
#' @return `df` with those columns renamed and any NAs in them replaced by 0.
#' @examples
#' landcover_key <- c(
#'   "1"="urban","2"="croplands","3"="forest","4"="shrub & scrub",
#'   "5"="grass","6"="water","7"="wetlands","8"="bare surface"
#' )
#' # df2 <- rename_and_fill_proportions(df, landcover_key)
#' @export
rename_and_fill_proportions <- function(df,
                                        key,
                                        prefix = "proportion_",
                                        keep_prefix = TRUE) {
  stopifnot(is.data.frame(df))

  pat <- paste0("^", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}])", "\\\\\\1", prefix), "([0-9]+)$")
  m    <- regexec(pat, names(df))
  hits <- regmatches(names(df), m)
  idx  <- vapply(hits, length, integer(1)) == 2

  if (!any(idx)) return(df)

  codes  <- vapply(hits[idx], function(z) z[2], character(1))
  labels <- unname(key[codes])

  keep <- !is.na(labels)
  if (!all(keep)) {
    idx[idx][!keep] <- FALSE
    codes  <- codes[keep]
    labels <- labels[keep]
  }

  new_names <- if (keep_prefix) paste0(prefix, labels) else labels
  names(df)[which(idx)] <- new_names

  fill_cols <- which(names(df) %in% new_names)
  for (j in fill_cols) {
    x <- df[[j]]
    if (is.numeric(x)) {
      x[is.na(x)] <- 0
    } else if (is.logical(x)) {
      x[is.na(x)] <- FALSE
    } else {
      x[is.na(x)] <- 0
      x <- as.numeric(x)
    }
    df[[j]] <- x
  }

  df
}

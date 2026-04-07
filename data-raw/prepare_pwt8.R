## Download and prepare Jan Ditzen's xtdcce2 example dataset
## Source: https://github.com/JanDitzen/xtdcce2

url <- "https://github.com/JanDitzen/xtdcce2/raw/master/xtdcce2_sample_dataset.dta"
tmp <- tempfile(fileext = ".dta")
download.file(url, destfile = tmp, mode = "wb")
pwt_raw <- haven::read_dta(tmp)

# Inspect
cat("Variables:", paste(names(pwt_raw), collapse = ", "), "\n")
cat("Countries:", length(unique(pwt_raw$country)), "\n")
cat("Years:", paste(range(pwt_raw$year), collapse = "-"), "\n")
cat("Rows:", nrow(pwt_raw), "\n")

# Construct d_log_rgdpo (first difference of log_rgdpo)
pwt8 <- as.data.frame(pwt_raw)

# The dataset may use 'id' instead of 'country'
if ("id" %in% names(pwt8) && !"country" %in% names(pwt8)) {
  pwt8$country <- as.character(pwt8$id)
} else {
  pwt8$country <- as.character(pwt8$country)
}
pwt8$year <- as.integer(pwt8$year)

pwt8 <- pwt8[order(pwt8$country, pwt8$year), ]
rownames(pwt8) <- NULL

# Panel-aware first difference
units <- unique(pwt8$country)
pwt8$d_log_rgdpo <- NA_real_
for (u in units) {
  idx <- which(pwt8$country == u)
  if (length(idx) < 2) next
  pwt8$d_log_rgdpo[idx[-1]] <- diff(pwt8$log_rgdpo[idx])
}

cat("d_log_rgdpo NAs:", sum(is.na(pwt8$d_log_rgdpo)), "\n")
cat("d_log_rgdpo non-NA:", sum(!is.na(pwt8$d_log_rgdpo)), "\n")

usethis::use_data(pwt8, overwrite = TRUE)

row_tibble <- function(x, col_names) {
  tibble::as_tibble(rbind(setNames(x, col_names)))
}

parse_format <- function(format, sample) {
  purrr::map2(
    strsplit(sample, ":", fixed = TRUE),
    strsplit(format, ":", fixed = TRUE),
    row_tibble)
}

parse_info <- function(info) {
  strsplit(info, ";", fixed = TRUE) %>%
    purrr::map(~row_tibble(sub("^.*=(.*)", "\\1", .x), sub("^(.*)=.*", "\\1", .x)))
}

split_tiers <- function(x) {
  name <- names(x)
  x$`1` <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[1])) %>%
    purrr::flatten_int()
  x$`2` <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[2])) %>%
    purrr::flatten_int()
  names(x)[2:3] <- paste0(name, "_", names(x)[2:3])
  x[2:3]
}

calc_counts <- function(x) {
  T_REF_COUNT_KEY <- paste0("T_", x$REF, "U_", x$TQSS_NT)
  T_ALT_COUNT_KEY <- paste0("T_", x$ALT, "U_", x$TQSS_NT)
  N_REF_COUNT_KEY <- paste0("N_", x$REF, "U_", x$TQSS_NT)
  N_ALT_COUNT_KEY <- paste0("N_", x$ALT, "U_", x$TQSS_NT)
  T_REF_COUNT <- x[[T_REF_COUNT_KEY]]
  T_ALT_COUNT <- x[[T_ALT_COUNT_KEY]]
  N_REF_COUNT <- x[[N_REF_COUNT_KEY]]
  N_ALT_COUNT <- x[[N_ALT_COUNT_KEY]]
  T_VAF <- T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT)
  N_VAF <- N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT)
  data.frame(T_VAF = T_VAF, N_VAF = N_VAF)
}

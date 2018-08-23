parse_vcf_alt1 <- function(x) {
    vcf <- readr::read_tsv(x, comment = "##")
    vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
    vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                         info = parse_info(INFO),
                         normal = parse_format(FORMAT, NORMAL),
                         tumour = parse_format(FORMAT, TUMOR))
    vcf <- tidyr::unnest(vcf, info)
    vcf <- tidyr::unnest(vcf, N = normal, T = tumour, .sep = "_")
    vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_[ACTG]U", vars = names(vcf)), split_tiers)
    vcf <- readr::type_convert(vcf)
    vcf <- tidyr::gather(vcf, allele, count, dplyr::matches("[TN]_[ACTG]U_[12]"))
    vcf <- tidyr::separate(vcf, allele, c("sample", "base", "tier"), sep = "[_U]+", remove = FALSE)
    vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
    vcf <- dplyr::mutate(vcf, T_REF_COUNT = count[REF == base & sample == "T" & tier == TQSS_NT],
                         T_ALT_COUNT = count[ALT == base & sample == "T" & tier == TQSS_NT],
                         N_REF_COUNT = count[REF == base & sample == "N" & tier == TQSS_NT],
                         N_ALT_COUNT = count[ALT == base & sample == "N" & tier == TQSS_NT],
                         T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                         N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
    vcf <- dplyr::select(vcf, -sample, -base, -tier)
    vcf <- tidyr::spread(vcf, allele, count)  # Takes up most of the time
    vcf
}

parse_vcf_alt2 <- function(x) {
    vcf <- readr::read_tsv(x, comment = "##")
    vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
    vcf <- dplyr::mutate(vcf, ALT = purrr::map_chr(strsplit(ALT, ","), 1),
                         info = parse_info(INFO),
                         normal = parse_format(FORMAT, NORMAL),
                         tumour = parse_format(FORMAT, TUMOR))
    vcf <- tidyr::unnest(vcf, info)
    vcf <- tidyr::unnest(vcf, N = normal, T = tumour, .sep = "_")
    vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_[ACTG]U", vars = names(vcf)), split_tiers)
    vcf <- readr::type_convert(vcf)
    vcf <- purrr::by_row(vcf, calc_counts)
    vcf <- tidyr::unnest(vcf, .out)
    vcf
}

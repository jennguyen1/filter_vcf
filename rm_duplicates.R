# Ensure matches on chr:pos, a1, a2 of dose/info with the snp file
# Jennifer N Nguyen

library(scriptR)
scriptrR::lib()

# command line args
option_list <- list(
	optparse::make_option("--prefix"),
	optparse::make_option("--vcf_dose"),
	optparse::make_option("--vcf_info"),
	optparse::make_option("--snp_file_name"),
	optparse::make_option("--output_dose"),
	optparse::make_option("--output_info")
)

opt <- parse_args(optparse::OptionParser(option_list=option_list))
start_logging()

# open snp file
snp <- read.in(opt$snp_file_name, read.table, header = FALSE, colClasses = "character")
colnames(snp)[1:3] <- c("id", "effect", "noneffect")
snp <- dplyr::select(snp, id, effect, noneffect)

# open data file (VCF dose/info)
dose <- read.in(opt$vcf_dose, read.table, header = FALSE, colClasses = "character")
colnames(dose)[3:5] <- c("id", "a1", "a2") # dose template: columns 3-5 (ID, REF, ALT)

info <- read.in(opt$vcf_info, read.table, header = FALSE, colClasses = "character")
colnames(info)[1:3] <- c("id", "a1", "a2") # info template: columns 1-3 (SNP, REF, ALT)

# create index for filter
dose$index <- 1:nrow(dose)
info$index <- 1:nrow(info)

# match with alleles (via merge), combine, rm dups, and filter data
do_merge <- function(d){
    assert_cols_in(d, c("id", "a1", "a2", "index"))
    d1 <- merge(d %>% mutate(effect = a1, noneffect = a2, a1 = NULL, a2 = NULL), snp, by = c("id", "effect", "noneffect"))
    d2 <- merge(d %>% mutate(effect = a2, noneffect = a1, a1 = NULL, a2 = NULL), snp, by = c("id", "effect", "noneffect"))

    c <- rbind.data.frame(d1, d2)
    c <- dplyr::distinct(c)

    keep_index <- c$index
    r <- subset(d, index %in% keep_index)
    r$index <- NULL
    return(r)
}

filter_dose <- do_merge(dose)
filter_info <- do_merge(info)

# checks on dimensions before writing out
assert_dim(function(x) x, c(nrow(dose), NA))(info) -> trash
assert_dim(function(x) x, c(nrow(filter_dose), NA))(filter_info) -> trash
assert_dim(function(x) x, c(NA, ncol(info) - 1))(filter_info) -> trash
assert_dim(function(x) x, c(NA, ncol(dose) - 1))(filter_dose) -> trash

suppressWarnings(write.table(filter_dose, file = opt$output_dose, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"))
suppressWarnings(write.table(filter_info, file = opt$output_info, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"))

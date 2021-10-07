library(ggplot2)
library(tools)

args <- (commandArgs(TRUE))
if (length(args) == 1) {
    sex.filename <- args[1]
} else {
    stop("Script accepts exactly 1 argument: filename with the sex calling results (e.g. 'output/test/gender.csv'")
}
extension <- file_ext(sex.filename)
output <- sub(extension, "png", sex.filename, fixed=TRUE)
palette <- c("XY" = "#1F78B4", "XX" = "#E31A1C", "XYY" = "#4DAF4A", "XXX" = "#FF7F00", "XXY" = "#984EA3", "X0" = "#F781BF")

data <- read.delim(sex.filename)

pl <- ggplot(data) + geom_point(aes(x=XAutoRatio, y=YAutoRatio, colour=Karyotype), shape=1, alpha=0.6) +
    theme_bw() + theme(aspect.ratio=1) + scale_colour_manual(values=palette)

ggsave(output, pl, width = 10, height = 10, units = "cm",
           antialias="gray", dpi=300, type="cairo")


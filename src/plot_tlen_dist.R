library(tidyverse)
args <- commandArgs(trailingOnly=T)
input <- args[1]
output <- args[2]

#  input <- "qc/TLEN_dist/SRR12362020_MappedOn_tair10_nuclear_sort.md.tlen_dist.txt"
#  output <- "qc/TLEN_dist/plot/sdfsdfsdf.txt"

output_log10 <-  str_replace(output, ".pdf$", "_log10.pdf")
output_png <-  str_replace(output, ".pdf$", ".png")
output_png_log10 <-  str_replace(output, ".pdf$", "_log10.png")

tlen_count <- read_delim(file=input, delim=" ", col_names=c("count", "tlen"))

p.tlen_dist <- ggplot(data=tlen_count, aes(x=tlen, y=count)) +
    geom_col() +
    theme(text=element_text(size=8, family="Helvetica"))

p.tlen_dist.log10 <- ggplot(data=tlen_count, aes(x=tlen, y=count)) +
    geom_col() +
    scale_y_continuous(trans="log10") +
    theme(text=element_text(size=8, family="Helvetica"))

pdf(file=output, width=2, height=2)
plot(p.tlen_dist)
dev.off()

pdf(file=output_log10, width=2, height=2)
plot(p.tlen_dist.log10)
dev.off()

png(file=output_png, width=2, height=2, units="in", res=300)
plot(p.tlen_dist)
dev.off()

png(file=output_png_log10, width=2, height=2, units="in", res=300)
plot(p.tlen_dist.log10)
dev.off()

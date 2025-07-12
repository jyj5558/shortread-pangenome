library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b <- a + theme_light()
b_lim <- a + theme_light() + xlim(0,2500)
ggsave("var_qual.png", plot = b, width = 8, height = 6)
ggsave("var_qual_lim.png", plot = b_lim, width = 8, height = 6)
capture.output(summary(var_qual$qual))

   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b <- a + theme_light()   
b_lim <- a + theme_light() + xlim(0,100)      
ggsave("var_depth.png", plot = b, width = 8, height = 6)
ggsave("var_depth_lim.png", plot = b_lim, width = 8, height = 6)
capture.output(summary(var_depth$mean_depth))


var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b <- a + theme_light()
ggsave("var_miss.png", plot = b, width = 8, height = 6)
capture.output(summary(var_miss$fmiss))

                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b <- a + theme_light()
ggsave("ind_depth.png", plot = b, width = 8, height = 6)
capture.output(summary(ind_depth$depth))


ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b <- a + theme_light()
ggsave("ind_miss.png", plot = b, width = 8, height = 6)
capture.output(summary(ind_miss$fmiss))

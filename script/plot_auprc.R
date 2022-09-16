#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--auprc", required = TRUE, help = "the path to the ggplot.tsv file where the first column is sample, the second column is tool, and the third column is auprc.")
parser$add_argument("--labels", required = FALSE, help = "the path to the label.lst file. If specified, will be used to label each figure")
parser$add_argument("--sig", action = "store_true", help = "if specified, add the significance label to the box/bar plot")
parser$add_argument("--box", action = "store_true", help = "if specified, generate a box plot")
parser$add_argument("--bar", action = "store_true", help = "if specified, generate a bar plot")
parser$add_argument("-f", "--xtickfont", default = 8, required = FALSE, help = "the font size for xticks")
parser$add_argument("-o", "--output", required = TRUE, help = "output file prefix")

args <- parser$parse_args()

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
library(gridExtra)
suppressPackageStartupMessages(library("rstatix"))

auprc_file = args$auprc
label_list = args$labels
produce_sig = args$sig
produce_boxplot = args$box
produce_barplot = args$bar
prefix = args$output
x_tick_font_size = args$xtickfont

df <- read.table(file = auprc_file, sep = '\t', header = TRUE)
labels <- readLines(label_list) 
# add ast4risk to show p-value
stat.test <- compare_means(auprc ~ tool, data = df, group.by = "sample", method = "t.test")
stat.test <- add_y_position(stat.test, fun = "max", data = df, formula = auprc ~ tool)
stat.test <- add_x_position(stat.test)
# stat.test <- df %>% group_by(sample) %>% t_test(auprc ~ tool) %>% add_significance("p")
# don't know what this step does
# stat.test <- stat.test %>% add_xy_position(fun = "mean_sd", x = "sample", dodge = 0.8)
# produce boxplot
if (produce_boxplot == TRUE)
{
    bp <- ggbarplot(df, x = "sample", y = "auprc", add = "mean_sd", color= "tool", position = position_dodge(0.8)) + scale_x_discrete(labels = labels)
    if (produce_sig == TRUE)
    {
        bp <- bp + stat_pvalue_manual(stat.test, x = "sample", label = "p.signif", tip.length = 0.01, bracket.nudge.y = 0.02, hide.ns = TRUE) + coord_cartesian(ylim = c(0, 1)) +
 theme(axis.text = element_text(size = x_tick_font_size))
    }
    ggsave(filename = paste(prefix, ".barplot.png", sep = ''), width = 5, height = 5, dpi = 320, plot = bp)
}
# produce bar plot
if (produce_barplot == TRUE)
{
    bxp <- ggboxplot(df, x = "sample", y = "auprc", color = "tool") + scale_x_discrete(labels = labels)
    if (produce_sig == TRUE)
    {
        bxp <- bxp + stat_pvalue_manual(stat.test, x = "sample", label = "p.signif", tip.length = 0.01, bracket.nudge.y = 0.03, hide.ns = TRUE) + coord_cartesian(ylim = c(0, 1)) + theme(axis.text = element_text(size = x_tick_font_size))
    }
    ggsave(filename = paste(prefix, ".boxplot.png", sep = ''), width = 5, height = 5, dpi = 320, plot = bxp)
    # ggexport(bxp, filename = paste(prefix, ".boxplot.png", sep = ''))
}



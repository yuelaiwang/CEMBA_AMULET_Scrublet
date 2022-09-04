#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--dir", required = FALSE, help = "the directory that contains the feather files")
parser$add_argument("--samples", required = TRUE, help = "the path to the sample.lst file")
parser$add_argument("--labels", required = FALSE, help = "the path to the label.lst file. If specified, will be used to label each figure")
# parser$add_argument("--labelx", default = 0.1, "the x position for the label of each sub figure")
parser$add_argument("--suffix", required = FALSE, help = "the suffix of each feather file")
parser$add_argument("-n", "--number_or_datasets", default = 10, type = "integer", help="the number of simulated datasets in a feather file")
parser$add_argument("--row", default = 2, required = FALSE, type = "integer", help = "the number of rows in the figure panel")
parser$add_argument("--column", default = 4, type = "integer", help="the number of columns in the figure panel")
parser$add_argument("-o", "--output", required = TRUE, help = "output file name")

args <- parser$parse_args()

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
library(feather)
library(gridExtra)
suppressPackageStartupMessages(library("tidyverse"))
# library("tictoc")

# feathers <- args$feather
input_dir = args$dir
if (!endsWith(input_dir, '/')) input_dir <- paste(input_dir, '/')
sample_list = args$samples
label_list = args$labels
# labelx = args$labelx
suffix = args$suffix
number_of_datasets = args$number_or_datasets
row_number = args$row
column_number = args$column
outF = args$output

samples <- readLines(sample_list)
if (length(label_list) != 0) labels <- readLines(label_list) else labels <- samples
# save every PRC to this list
plots = vector("list", length(samples))
# plots = list()
if (length(samples) > row_number * column_number) stop("Please reshape the image row/column numbers! Sample size larger than image size.")

for (i in 1 : length(samples))
{
    feather_file <- paste(input_dir, samples[i], suffix, sep = '')
    df <- read_feather(feather_file)
    number_of_data_points <- nrow(df) / number_of_datasets
    # create a new dataframe with an extra column specifying the id of the dataset
    unique_df <- data.frame(matrix(ncol = 4, nrow = 0))
    x <- c("recall", "tool", "precision", "id")
    colnames(unique_df) <- x
    # add id to each dataset
    for (j in 1 : number_of_datasets)
    {
        sub_df <- df[as.integer(((j - 1) * number_of_data_points + 1)) : as.integer(j * number_of_data_points), ]
        sub_df <- sub_df %>% distinct()
        sub_df$id = j
        unique_df <- rbind(unique_df, sub_df)
        # sub_plot <- sub_plot + geom_line(sub_df, aes(x = recall, y = precision, color = tool))
    }
    plots[[i]] <- local({
        i <- i
        sub_plot <- ggplot(unique_df, aes(x = recall, y = precision, color = tool, group = interaction(tool, id))) + geom_line(aes(color = tool)) + theme(axis.title = element_text(size = 50), axis.text = element_text(size = 50), legend.text = element_text(size = 120), legend.key.size = unit(2, 'cm')) + coord_cartesian(ylim = c(0, 1)) + theme_classic()
        print(sub_plot) 
     })
    # sub_plot <- ggplot(unique_df, aes(x = recall, y = precision, color = tool, group = interaction(tool, id))) + geom_line(aes(color = tool)) + theme_classic() 
    # append(plots, sub_plot)
}
# arrange all the PRCs 
# final_plot <- do.call("grid.arrange", c(plots, nrow = row_number, ncol = column_number)) 
final_plot <- ggarrange(plotlist = plots, nrow = row_number, ncol = column_number, common.legend = TRUE, legend = "right", labels = labels, font.label = list(size = 10), label.x = 0.4)  

# save figure
ggexport(final_plot, width = column_number * 3.2, height = row_number * 3, filename = outF)

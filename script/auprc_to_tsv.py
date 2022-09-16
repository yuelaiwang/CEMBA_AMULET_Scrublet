import pandas as pd
import numpy as np
import seaborn as sns
from os.path import exists
import argparse

parser = argparse.ArgumentParser(description = "convert the cemba.mop.prc.auc.txt file to a tsv file where the first column is sample, and the second/third columns are the names of the tool. The tsv file would be further converted to an Excel file to make box plot.")
parser.add_argument("auc", type = str, help = "the path to the prc.auc.txt file. The number of lines in this file should be the product of that of sample.lst and that of tools.lst")
parser.add_argument("samples", type = str, help = "the path to the sample.lst file")
parser.add_argument("tools", type = str, help = "the path to the tools.lst file")
parser.add_argument("output", type = str, help = "the name of the output file")
parser.add_argument("--ggplot", action = "store_true", help = "whether to generate a tsv for ggplot") 
args = parser.parse_args()
with open(args.auc, 'r') as file3:
    aucs = list()
    for line in file3:
        aucs.append(line.strip())
with open(args.samples, 'r') as file1:
    samples = list()
    for line in file1:
        samples.append(line.strip())
with open(args.tools, 'r') as file2:
    tools = list()
    for line in file2:
        tools.append(line.strip())

def main1():
    if len(tools) == 2:
        matrix = pd.DataFrame({"sample": [], tools[0] : [], tools[1] : []})
    else:
        matrix = pd.DataFrame({"sample": [], tools[0] : []})
    '''
    fill the tsv table
    '''
    for i in range(len(samples)):
        sample = samples[i]
        if len(tools) == 2:
            auprc1 = aucs[i * 2].split(',')
            auprc2 = aucs[i * 2 + 1].split(',')
            for k in range(len(auprc1)):
                matrix.loc[len(matrix)] = [sample, auprc1[k], auprc2[k]]
        else:
            auprc1 = aucs[i]
            for k in range(len(auprc1)):
                matrix.loc[len(matrix)] = [sample, auprc1[k]]

    matrix.to_csv(path_or_buf = args.output, sep = '\t', index = False)

'''
output a tsv that has three columns: sample, tool, auprc
'''
def main2():
    matrix = pd.DataFrame({"sample": [], "tool" : [], "auprc" : []})
    for i in range(len(samples)):
        sample = samples[i]
        if len(tools) == 2:
            auprc1 = aucs[i * 2].split(',')
            auprc2 = aucs[i * 2 + 1].split(',')
            for k in range(len(auprc1)):
                matrix.loc[len(matrix)] = [sample, tools[0], auprc1[k]]
            for k in range(len(auprc2)):
                matrix.loc[len(matrix)] = [sample, tools[1], auprc2[k]]
        else:
            auprc1 = aucs[i]
            for k in range(len(auprc1)):
                matrix.loc[len(matrix)] = [sample, tools[0], auprc1[k]]
    matrix.to_csv(path_or_buf = args.output, sep = '\t', index = False)

if not args.ggplot:
    main1()
else:
    main2()


import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import feather
from os.path import exists
import argparse

parser = argparse.ArgumentParser(description = "Take a set of simulated datasets (replicates of each other), calculate coordinates of the precision-recall curve for each dataset depicting the performance of AMULET and Scrublet, respectively, in order to compare them. Options can be specified to make a panel of preliminary PRCs or store the coordinates into a feather file as an input to R.")
parser.add_argument("datasets", type = str, help = "the path to the simulated dataset")
parser.add_argument("--prefix", type = str, required = False, help = "the prefix of where to store the image file")
parser.add_argument("--tools", type = int, default = 2, help = "determine how many tools are applied on each simulated dataset, 2 means both, and 1 means AMULET only")
parser.add_argument("--end", type = int, default = 9, help = "the last dataset to compare the performance of AMULET and Scrublet")
parser.add_argument("--start", type = int, default = 1, help = "the first dataset to compare the performance of AMULET and Scrublet.")
parser.add_argument("--AMULETorder", type = int, default = -20, help = "a negative integer representing the starting order of magnitude for q-value of AMULET results")
parser.add_argument("--highqvalue", type = float, default = 1.01, help = "the upper bound of the q-value for AMULET results")
parser.add_argument("--scrhigh", type = float, default = 0.8, help = "a float between zero and one representing the highest doublet score for Scrublet results")
parser.add_argument("--scrlow", type = float, default = 0.001, help = "a float between zero and one representing the lowest doublet score for Scrublet results")
parser.add_argument("--scrstep", type = float, default = -0.002, help = "Scrublet step")
parser.add_argument("-f", "--feather", type = str, required = False, help = "if specified, output the dataframe as feather format to this path")

args = parser.parse_args()
path = args.datasets
prefix = args.prefix
numberOfTools = args.tools
endingDataset = args.end
startingDataset = args.start
AMULETstart = args.AMULETorder
AMULEThigh = args.highqvalue
ScrubletHigh = args.scrhigh
ScrubletLow = args.scrlow
ScrubletStep = args.scrstep
feather_path = args.feather

numberOfDatasets = endingDataset - startingDataset + 1
numberOfDatasetsToInclude = numberOfDatasets
rowNumber = math.floor(math.sqrt(numberOfDatasetsToInclude))
columnNumber = math.floor(math.sqrt(numberOfDatasetsToInclude))
if rowNumber * columnNumber < numberOfDatasets:
    columnNumber += 1
if rowNumber * columnNumber < numberOfDatasets:
    rowNumber += 1

def GradientCutoffs(pathToSimulatedDatasets):
    #if not pathToSimulatedDatasets.endswith('/'):
    #    pathToSimulatedDatasets = pathToSimulatedDatasets + '/'
    f = plt.figure(figsize =(20, 15))
    raw_data_10_datasets = pd.DataFrame({"recall" : [], "tool": [], "precision" : []})
    datasetSkipped = 0
    '''
    Use dataset 1-9 (by default) in simulated_datasetsID
    '''
    for j in range(numberOfDatasets):
        sns.set(font_scale = 1.35)
        f.add_subplot(rowNumber, columnNumber, j + 1)
        groundTruth = pd.read_csv(pathToSimulatedDatasets + str(j + startingDataset) + "/ground.truth.tsv", sep = '\t')
        if exists(pathToSimulatedDatasets + str(j + startingDataset) + "/MultipletProbabilities.txt"):
            multipletProbabilities = pd.read_csv(pathToSimulatedDatasets + str(j + startingDataset) + "/MultipletProbabilities.txt", sep = '\t')
        else:
            datasetSkipped += 1
            continue
        if numberOfTools == 2:
            scrubletScores = pd.read_csv(pathToSimulatedDatasets + str(j + startingDataset) + "/simulated.gmat.fitDoublets.txt", sep = '\t')
        doublets = set(groundTruth.loc[:, "doubletBarcode"])
        table = {"recall" : [], "tool": [], "precision" : []}
        rawData = pd.DataFrame(table)
        '''
        add AMULET data with a default start equal to -20, i.e. 10 ** -20
        '''
        for i in range(AMULETstart, -2):
            for k in range(1, 10):
                cutoff = 10 ** i * k
                lines = multipletProbabilities[multipletProbabilities["q-value"] <= cutoff]
                positives = set(lines.loc[:, "cell_id"])
                truePositives = doublets.intersection(positives)
                recall = len(truePositives) / len(doublets)
                if len(positives) == 0:
                    precision = 1
                else:
                    precision = len(truePositives) / len(positives)
                rawData.loc[len(rawData)] = [recall, "AMULET", precision]
                raw_data_10_datasets.loc[len(raw_data_10_datasets)] = [recall, "AMULET", precision]
        for cutoff in np.arange(0.01, AMULEThigh, 0.01):
            lines = multipletProbabilities[multipletProbabilities["q-value"] <= cutoff]
            positives = set(lines.loc[:, "cell_id"])
            truePositives = doublets.intersection(positives)
            recall = len(truePositives) / len(doublets)
            precision = len(truePositives) / len(positives)
            rawData.loc[len(rawData)] = [recall, "AMULET", precision]
            raw_data_10_datasets.loc[len(raw_data_10_datasets)] = [recall, "AMULET", precision]
        '''
        add Scrublet data; default = np.arange(0.8, 0.001, -0.002)
        '''
        if numberOfTools == 2:
            for score in np.arange(ScrubletHigh, ScrubletLow, ScrubletStep):
                lines = scrubletScores[scrubletScores["doublet_scores"] >= score]
                positives = set(lines.loc[:, "barcode"])
                truePositives = doublets.intersection(positives)
                recall = len(truePositives) / len(doublets)
                if len(positives) == 0:
                    precision = 1
                else:
                    precision = len(truePositives) / len(positives)
                rawData.loc[len(rawData)] = [recall, "Scrublet", precision]
                raw_data_10_datasets.loc[len(raw_data_10_datasets)] = [recall, "Scrublet", precision]
        rawData = rawData.drop_duplicates()
        ax = sns.lineplot(data = rawData, x = "recall", y = "precision", hue = "tool")
        ax.set_ylim([0, 1])
        ax.set_title("dataset" + str(j + startingDataset))
    plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
    '''
    if pathToSimulatedDatasets.endswith('_'):
        plt.savefig(pathToSimulatedDatasets + "performance_prc.png")
    else: 
        plt.savefig('/'.join(pathToSimulatedDatasets.split('/')[:-1]) + "/performance_prc.png")    
    '''
    if prefix != None:
        plt.savefig(prefix + "performance_prc.png")
    if feather_path != None:
        feather.write_dataframe(raw_data_10_datasets, feather_path)
    return raw_data_10_datasets

_ = GradientCutoffs(path)

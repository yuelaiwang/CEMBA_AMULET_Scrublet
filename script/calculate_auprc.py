import math
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description = "Take a set of simulated datasets (technical replicates of each other), calculate the area under the PRC for each dataset depicting the performance of AMULET and Scrublet, respectively, and write the values of the auc to a text file")
parser.add_argument("datasets", type = str, help = "the prefix of a simulated dataset, like \"$cemba_mop/simulated_dataset_rep1/dataset11_\"")
parser.add_argument("prefix", type = str, help = "the prefix of where to store the resulting text file")
parser.add_argument("--tools", type = int, default = 2, help = "determine how many tools are applied on each simulated dataset, 2 means both, and 1 means AMULET only")
parser.add_argument("--end", type = int, default = 10, help = "the last dataset to compare the performance of AMULET and Scrublet")
parser.add_argument("--start", type = int, default = 1, help = "the first dataset to compare the performance of AMULET and Scrublet.")
parser.add_argument("--AMULETorder", type = int, default = -20, help = "a negative integer representing the starting order of magnitude for q-value of AMULET results")
parser.add_argument("--highqvalue", type = float, default = 1.01, help = "the upper bound of the q-value for AMULET results")
parser.add_argument("--scrhigh", type = float, default = 0.8, help = "a float between zero and one representing the highest doublet score for Scrublet results")
parser.add_argument("--scrlow", type = float, default = 0.001, help = "a float between zero and one representing the lowest doublet score for Scrublet results")
parser.add_argument("--scrstep", type = float, default = -0.002, help = "Scrublet step")

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

numberOfDatasets = endingDataset - startingDataset + 1

def main():
    tables = GradientCutoffs(path)
    AMULETauc = list()
    if numberOfTools == 2:
        Scrubletauc = list()
    for table in tables:
        auc1, auc2 = CalculateAUC(table)
        AMULETauc.append(auc1)
        if numberOfTools == 2:
            Scrubletauc.append(auc2)
    '''
    append to a text file where every two lines represent a set of simulated datasets with the first line representing the AUCs for AMULET and second for Scrublet.
    '''
    file1 = open(prefix + ".auc.txt", 'a')
    for i in range(len(AMULETauc) - 1):
        file1.write(str(AMULETauc[i]))
        file1.write(',')
    file1.write(str(AMULETauc[-1]))
    file1.write('\n')
    if numberOfTools == 2: 
        for i in range(len(Scrubletauc) - 1):
            file1.write(str(Scrubletauc[i]))
            file1.write(',')
        file1.write(str(Scrubletauc[-1]))
        file1.write('\n')
    file1.close()

# takes simulated dataset 1-10 (by default), returns a list of pandas tables, each containing scatter points for the PRC for that simulated dataset 
def GradientCutoffs(pathToSimulatedDatasets) -> list:
    #if not pathToSimulatedDatasets.endswith('/'):
    #    pathToSimulatedDatasets = pathToSimulatedDatasets + '/'
    '''
    Use dataset 1-9 in simulated_datasetsID
    '''
    rv = list()
    for j in range(numberOfDatasets):
        groundTruth = pd.read_csv(pathToSimulatedDatasets + str(j + startingDataset) + "/ground.truth.tsv", sep = '\t')
        multipletProbabilities = pd.read_csv(pathToSimulatedDatasets + str(j + startingDataset) + "/MultipletProbabilities.txt", sep = '\t')
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
                # rawData.loc[len(rawData)] = [recall, "precision", falsePositive]
        for cutoff in np.arange(0.01, AMULEThigh, 0.01):
            lines = multipletProbabilities[multipletProbabilities["q-value"] <= cutoff]
            positives = set(lines.loc[:, "cell_id"])
            truePositives = doublets.intersection(positives)
            recall = len(truePositives) / len(doublets)
            precision = len(truePositives) / len(positives)
            rawData.loc[len(rawData)] = [recall, "AMULET", precision]
        '''
        add Scrublet data; default = np.arange(0.5, 0.006, -0.002)
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
            
        rawData = rawData.drop_duplicates()
        rv.append(rawData)
    # plt.savefig(pathToSimulatedDatasets + pathToSimulatedDatasets[(-14 - 1):-1] + "_performance_prc.png")
    
    return rv

# takes a pandas table with three columns (recall, tool, & precision), calculate the area under the curve for AMULET & Scrublet, respectively
def CalculateAUC(table):
    AMULETAUC = 0
    ScrubletAUC = 0
    AMULET = table[table["tool"] == "AMULET"]
    if numberOfTools == 2:
        Scrublet = table[table["tool"] == "Scrublet"]
    '''
    iterate through every scatter point and compute the area of the trapezium formed by adjacent points
    '''
    x0 = AMULET.loc[AMULET.index[0]]["recall"]
    y0 = AMULET.loc[AMULET.index[0]]["precision"]
    AMULETAUC += (1 + y0) * x0 / 2
    for i in range(len(AMULET.index) - 1):
        x1 = AMULET.loc[AMULET.index[i]]["recall"]
        y1 = AMULET.loc[AMULET.index[i]]["precision"]
        x2 = AMULET.loc[AMULET.index[i + 1]]["recall"]
        y2 = AMULET.loc[AMULET.index[i + 1]]["precision"]
        area = (y1 + y2) * (x2 - x1) / 2
        AMULETAUC += area
    AMULETAUC += y2 * (1 - x2) / 2
    '''
    iterate through every scatter point and compute the area of the trapezium formed by adjacent points
    '''   
    if numberOfTools == 2:
        x0 = Scrublet.loc[Scrublet.index[0]]["recall"]
        y0 = Scrublet.loc[Scrublet.index[0]]["precision"]
        ScrubletAUC += (1 + y0) * x0 / 2
        for i in range(len(Scrublet.index) - 1):
            x1 = Scrublet.loc[Scrublet.index[i]]["recall"]
            y1 = Scrublet.loc[Scrublet.index[i]]["precision"]
            x2 = Scrublet.loc[Scrublet.index[i + 1]]["recall"]
            y2 = Scrublet.loc[Scrublet.index[i + 1]]["precision"]
            area = (y1 + y2) * (x2 - x1) / 2
            ScrubletAUC += area
        ScrubletAUC += y2 * (1 - x2) / 2    
    return AMULETAUC, ScrubletAUC    
main()

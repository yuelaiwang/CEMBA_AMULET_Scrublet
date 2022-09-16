import random
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description = "Generate the stat.txt needed by the script /projects/ps-renlab/yangli/projects/CEMBA/00.data/szu/bin/snapATAC.qc.filter.R. Theoretically there is no need to generate stat.txt to run the script again, because barcodes in the simulated datasets all pass the TSSE and UQ thresholds. However, I decide to just follow the whole pipeline.")
parser.add_argument('stat', type=str, help='the path to the stat.txt file for this region')
parser.add_argument('doublets', type=str, help='the path to the \"ground.truth.tsv\" file')
parser.add_argument("--location", type = str, default = "None", help = "the path to the directory that stores the resulting stat.txt file. If not specified, store the resulting file under the same directory as the ground truth file")
args = parser.parse_args()

def main():
    #stat = args.stat
        #location = location + '/'
    #doubletsFile = open(args.doublets, 'r')
    groundTruth = pd.read_csv(args.doublets, sep = '\t')
    doublets = list(groundTruth.loc[:, "doubletBarcode"])
    # check where to write the resulting file
    if args.location == "None":
        location = '/'.join(args.doublets.split('/')[:-1]) + '/'
    else:
        location = args.location
        if not location.endswith('/'):
            location = location + '/'
    statFile = open(location + "stat.txt", 'w')
    
    stats = open(args.stat, 'r').readlines()
    for stat in stats:
        statFile.write(stat)
    
    for doublet in doublets:
        doublet = doublet.strip()
        randomTSSE = random.uniform(10,15)
        newRecord = doublet + "\t0\t0\t" + str(randomTSSE) + "\t1001\t0\n"
        statFile.write(newRecord)
    statFile.close()
    # stats.close()

main()

import numpy
import pandas as pd
import seaborn as sns
import argparse
parser = argparse.ArgumentParser(description = "Find barcodes that are identified as singlets by both AMULET and Scrublet in a cerebrum region.")
parser.add_argument('--meta', type=str, dest="meta", help='the path to the filtered meta file of the cerebrum region')
parser.add_argument('--scrublet', type=str, dest="scrublet", help='the path to the file \"gmat.fitDoublets.txt\" ')
parser.add_argument('--AMULET', type=str, dest="AMULET", help='the path to the file \"MultipletBarcodes_01.txt\"')
parser.add_argument('--output', type=str, dest= "outfile", help = "the path to the file that stores the singlet barcodes.")
args = parser.parse_args()

def main():
    meta = args.meta
    scrublet = args.scrublet
    AMULET = args.AMULET
    outfile = args.outfile

    cells = pd.read_table(meta)
    scrubletResult = pd.read_table(scrublet)

    AMULETDoublets = open(AMULET, 'r').readlines()
    AMULETDoublets = list(i.strip() for i in AMULETDoublets)

    scrubletDoublets = ScrubletDoublets(scrubletResult)
    allBarcodes = AllBarcodes(cells)

    commonSinglets = open(outfile, 'w')
    for i in allBarcodes:
        if i not in scrubletDoublets and i not in AMULETDoublets:
            commonSinglets.write(i)
            commonSinglets.write('\n')
    commonSinglets.close()        
# return a list of barcodes of doublets found by Scrublet for this region
def ScrubletDoublets(scrubletResult) -> list:
    scrubletDoublets = list()
    for i in range(len(scrubletResult)):
        if scrubletResult.loc[i, "predicted_doublets"] == True:
            scrubletDoublets.append(scrubletResult.loc[i, "barcode"])
    return scrubletDoublets

# return a list of all barcodes of this region (after filtering)
def AllBarcodes(cells) -> list:
    return list(cells.loc[:, "barcode"])

main()

import pandas as pd
import numpy as np
import pysam
import random
import argparse
parser = argparse.ArgumentParser(description = "Generate bam files each containing singlets and simulated doublets merged from singlets. For example, if there are 8800 singlets, and the proportion of doublets is set to be 0.1, then 1600 random singlets would be selected to generate 800 doublets, and the final bam file would contain 8000 cells in total with 7200 singlets and 800 doublets.")
parser.add_argument('bam', type = str, help='the path to the parent bam file')
parser.add_argument('annotation', type = str, help = "the path to an tsv annotation file, here is metatable.tsv")
parser.add_argument("classifier", type = str, help = "the column name in the annotation file used to classify cell type for doublet annotation as homo/hetero")
parser.add_argument("singlets", type = str, help='the path to the \"SingletBarcodes_01.txt\" file')
parser.add_argument("summary", type = str, help = 'the path to the \"OverlapSummary.txt\" file given by running AMULET on the original bam file of this region/replication. This is used to record the number of uniquely aligned fragments for each singlet')
parser.add_argument("--fclower", type = int, dest = "low", default = 0, help = "an integer number that the script would only pick cells with a fragment counts above")
parser.add_argument("--fchigher", type = int, dest = "high", default = 100000000, help = "an integer number that the script would only pick cells with a fragment counts below")
parser.add_argument("--poolsize", type = int, dest = "size", default = 100000, help = "the number of singlets to keep in the pool (when downsampling, need to keep the pool size the same)")
parser.add_argument('--proportion', type=float, dest="proportion", default = 0.1, help='the proportion of doublets in a simulated dataset')
parser.add_argument('--homoratio', type = str, dest = "homoratio", default = "False", help = "whether to generate constant homotypic doublet ratio datasets, TRUE means yes and FALSE means no")
parser.add_argument('location', type=str, help='the location to store generated bam files and corresponding ground truth file of simulated doublet barcodes. The groung truth file will be a 6-column tsv file whose columns are doublet_barcode, singlet1_barcode, singlet2_barcode, singlet1_majortype, singlet2_majortype, and doublet_type.')
parser.add_argument("--cellID", type = str, default = "Barcode", help = "the column in the CEMBA 1.0 metatable that is used here")
parser.add_argument('--region', type=str, dest = "region", default = "None", help = 'the region name, e.g.CEMBA171213_4B. This is totally optional, specifying it would speed up the execution of the script a bit.')

args = parser.parse_args()
cellID = args.cellID
def main():
    bam = args.bam
    annotation = args.annotation
    singletsFile = args.singlets
    summaryFile = args.summary
    FCLow = args.low
    FCHigh = args.high
    singletPoolSize = args.size
    proportion = args.proportion
    region = args.region
    classifier = args.classifier
    homoratio = args.homoratio
    location = args.location
    # cellID = args.cellID
    if not location.endswith('/'):
        location = location + '/'

    lines = open(singletsFile, 'r').readlines()
    singlets2 = list(singlet.strip() for singlet in lines)
    overlapSummary = pd.read_csv(summaryFile, sep = '\t')
    # make a dictionary whose keys are barcodes and values are read counts so that when filtering cells with low read counts, there is no need to go through the pandas table every time.
    barcodeReadCounts = dict()
    for index in range(len(overlapSummary)):
        barcode = overlapSummary.loc[index, "Cell Id"]
        fragmentCounts = overlapSummary.loc[index, "Number of Valid Read Pairs"]
        barcodeReadCounts[barcode] = fragmentCounts
    # use the metatable from CEMBA 1.0, restrict the singlet pool to cells in the table only.
    cellType = annotateCellType(annotation, classifier, region)
    # use annotated cells and cells with a read counts/fragment counts higher than a threshold. This step is important: if a threshold is specified, cells with low read/fragment counts are not included as cells anymore, and they won't affect the background of AMULET in downstream analysis
    singlets = list() 
    for singlet in singlets2:
        if singlet in cellType and barcodeReadCounts[singlet] >= FCLow and barcodeReadCounts[singlet] <= FCHigh:
            singlets.append(singlet)
    # code added above to filter cells below the specified read count threshold

    # check if user wants to downsample the size of the  singlets
    if singletPoolSize < 100000:
        if singletPoolSize < len(singlets):
            print("downsampling to " + str(singletPoolSize) + " number of singlets")
            indexToKeep = random.sample(range(0, len(singlets)), singletPoolSize)
            downsampledSinglets = list()
            for index in indexToKeep:
                downsampledSinglets.append(singlets[index])
            singlets = downsampledSinglets
        # if the user-specified size of singlet pool size is even larger than the number of singlets in this region
        else:
            print("trying to down sample to " + str(singletPoolSize) + " number of singlets, but there are only " + str(len(singlets))  + " singlets. Stop downsampling. ")
    else:
        print("no downsampling required by user.")
    artificialBarcodes = GenerateBarcodes(singlets, cellType, proportion)
    artificialBarcodes.to_csv(path_or_buf = location + "ground.truth.tsv", sep = '\t')
    '''
    regenerate the singlecell.csv file used in AMULET
    '''
    # make a dictionary in which keys are chosen singlets and values are artifical doublets. This would reduce runtime by avoiding going over the whole artificialBarcodes pandas table every time.
    singletsToDoublets = dict()
    for index in range(len(artificialBarcodes)):
        singletsToDoublets[artificialBarcodes.loc[index, "singlet1Barcode"]] = artificialBarcodes.loc[index, "doubletBarcode"]
        singletsToDoublets[artificialBarcodes.loc[index, "singlet2Barcode"]] = artificialBarcodes.loc[index, "doubletBarcode"]
    file2 = open(location + "simulated.singlecell.csv", 'w')
    # uniqueSet = WholeSetOfBarcodes(artificialBarcodes, singlets)
    for barcode in singlets:
        if barcode not in singletsToDoublets:
            line = ",".join([barcode, barcode, "1"]) + '\n'
        else:
            line = ",".join([barcode, singletsToDoublets[barcode],"1" ]) + '\n'
        file2.write(line)
    file2.close()
    # generate the bam file for the simulated dataset for efficient conversion to a snap file
    samfile = pysam.AlignmentFile(bam, "rb")    
    simulation = pysam.AlignmentFile(location + "simulated.bam", "wb", template=samfile)
    GenerateSimulatedDataset(samfile, simulation, singlets, artificialBarcodes)
    simulation.close()
    samfile.close()

# take the path to the annotation file and return a dictionary that contains the annotation of majortype of the cells in this region. This function would collect information of the type of simulated doublets, homo/hetero. The type would be reflected in the ground.truth.tsv file.
def annotateCellType(annotationFile, classifier, region) -> dict:
    cellTypeAnnotation = dict() # keys are barcodes, and values are major type
    annotationTable = pd.read_csv(annotationFile, sep = '\t')
    if region != "None":
        regionalTable = annotationTable[annotationTable["Sample"] == region] # the name "Sample" corresponds to the column name in the metatable
    else:
        regionalTable = annotationTable
    for i in regionalTable.index:
        barcode = regionalTable.loc[i, cellID] # cellID is customizable
        cellType = regionalTable.loc[i, classifier]
        cellTypeAnnotation[barcode] = cellType
    return cellTypeAnnotation

# This function generates artificial barcodes. It stores the information into a 6-column pandas dataframe, where the columns are, in order, artificial doublet barcode, singlet1 barcode, singlet2 barcode, singlet1 major type, singlet2 major type, and doublet type. The number of artificial barcodes to be generated is N(singlets) * proportion / (1 + proportion). The dataframe would be stored as a tsv file called ground.truth.tsv after running the script. 
def GenerateBarcodes(singlets : list, cellType : dict, proportion) -> dict:
    multiDict = {"doubletBarcode" : [], "singlet1Barcode" : [], "singlet2Barcode": [], "singlet1MajorType" : [], "singlet2MajorType" : [], "doubletType" : []}
    artificialBarcodes = pd.DataFrame(multiDict)
    nucleotides = ['A', 'C', 'G', 'T']
    lengthOfBarcodes = len(singlets[0])
    numberOfArtificialBarcodes = int(proportion * len(singlets) / (1 + proportion))
    # generate a list of random indices in singlets. These singlets will be replaced by doublets
    singletsToReplace = random.sample(range(0, len(singlets)), 2 * numberOfArtificialBarcodes)
    for j in range(numberOfArtificialBarcodes):
        artificialBarcode = ""
        ''' generate a random barcode of the same length as in singlets '''
        for i in range(lengthOfBarcodes):
            index = random.randint(0, 3)
            artificialBarcode = artificialBarcode + nucleotides[index]
        # assign values of two real barcodes to be an artificial barcode
        singlet1 = singlets[singletsToReplace[j * 2]]
        singlet2 = singlets[singletsToReplace[j * 2 + 1]]
        singlet1type = cellType[singlet1]
        singlet2type = cellType[singlet2]
        doubletType = "heterotypic"
        if singlet1type == singlet2type:
            doubletType = "homotypic"
        artificialBarcodes.loc[len(artificialBarcodes.index)] = [artificialBarcode, singlet1, singlet2, singlet1type, singlet2type, doubletType]

    return artificialBarcodes

# This function generates a list of barcodes in this region after replacing some barcodes with artificial barcodes. This function is used for regenerating the singlecell.csv file for running AMULET.
def WholeSetOfBarcodes(artificialBarcodes, singlets : list) -> list:
    wholeset = list()
    singlet1Barcodes = list(artificialBarcodes["singlet1Barcode"])
    singlet2Barcodes = list(artificialBarcodes["singlet2Barcode"])
    for singlet in singlets:
        if singlet in singlet1Barcodes:
            index = singlet1Barcodes.index(singlet)
            wholeset.append(artificialBarcodes.loc[index, "doubletBarcode"])
        elif singlet in singlet2Barcodes:
            index = singlet2Barcodes.index(singlet)
            wholeset.append(artificialBarcodes.loc[index, "doubletBarcode"])
            # wholeset.append(artificialBarcodes[singlet])
        else:
            wholeset.append(singlet)
    uniqueSet = list(set(wholeset))
    return uniqueSet

# generate a bam file that contains artificial doublet barcodes. These barcodes replace singlet barcodes that are picked to form the artificial doublets.
def GenerateSimulatedDataset(samfile, simulation, singlets, artificialBarcodes):
    singlet1Barcodes = list(artificialBarcodes["singlet1Barcode"])
    singlet2Barcodes = list(artificialBarcodes["singlet2Barcode"])
    for record in samfile.fetch(until_eof = True): # don't have to make index file for the bam file
        barcode = record.query_name.split(':')[0]
        newRecord = record
        if barcode in singlet1Barcodes: # if the barcode is a singlet to be merged into a doublet
            index = singlet1Barcodes.index(barcode)
            artificialBarcode = artificialBarcodes.loc[index, "doubletBarcode"] # obtain the artificial doublet barcode
            newRecord.set_tag("CB", artificialBarcode) # set CB attribute to be the artificial barcode
            queryName = newRecord.query_name
            queryNameSuffix = queryName[len(barcode):]
            newRecord.query_name = artificialBarcode + queryNameSuffix # fix the query name, the first field of a bam file, to have proper input for Scrublet
            simulation.write(newRecord)
        elif barcode in singlet2Barcodes:
            index = singlet2Barcodes.index(barcode)
            artificialBarcode = artificialBarcodes.loc[index, "doubletBarcode"]
            newRecord.set_tag("CB", artificialBarcode)
            queryName = newRecord.query_name
            queryNameSuffix = queryName[len(barcode):]
            newRecord.query_name = artificialBarcode + queryNameSuffix
            simulation.write(newRecord)
        elif barcode in singlets:
            simulation.write(record)
main()

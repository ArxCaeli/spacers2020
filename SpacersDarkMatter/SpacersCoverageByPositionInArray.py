import os
import random

CRISPRInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_info_known_021718_len.tsv"
ArraysToSpacersQtyDict = dict()
for Line in open(CRISPRInfoFileName):
    LineValues = Line[:-1].split("\t")
    ArraysToSpacersQtyDict["_".join(LineValues[1:4])] = int(LineValues[8])

GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/Linked.tsv"

ArraysInfoDict = dict()
for Line in open(GenomeToArraysFileName):
    Arrays = Line[:-1].split("\t")[1]
    for Array in Arrays.split(","):
        ArrayID = "_".join(Array.split("_")[:-1])
        Type = Array.split("_")[-1]
        ArraysInfoDict[ArrayID] = [Type, ArraysToSpacersQtyDict[ArrayID]]

# print(len(ArraysDict))
# print(ArraysDict)

FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"
ArrayHitsInfo = dict()
for FileName in os.listdir(FolderName):
    if not "KPletsCoverageHitsLinked_" in FileName:
        continue

    for Line in open(FolderName + FileName):
        #CP001891.1_998980_1000411_21_spacer_1000229_32  GCCAGTCCGCATCTCGCCAAAA  1       Genome          Viral   NC_004775.2:34952,-
        LineValues = Line[:-1].split("\t")
        ArrayID = "_".join(LineValues[0].split("_")[0:3])

        if not ArrayID in ArrayHitsInfo:
            ArrayHitsInfo[ArrayID] = dict()

        KPletSize = len(LineValues[1])
        if not KPletSize in ArrayHitsInfo[ArrayID]:
            ArrayHitsInfo[ArrayID][KPletSize] = [set(), set(), set(), set()] # Spacers with hits into genome / viral / rand genome / rand / viral

        if LineValues[4] != "":
            if not "_Random" in LineValues[0]:
                ArrayHitsInfo[ArrayID][KPletSize][0].add(LineValues[0])
            else:
                ArrayHitsInfo[ArrayID][KPletSize][2].add(LineValues[0])

        if LineValues[6] != "":
            if not "_Random" in LineValues[0]:
                ArrayHitsInfo[ArrayID][KPletSize][1].add(LineValues[0])
            else:
                ArrayHitsInfo[ArrayID][KPletSize][3].add(LineValues[0])

## arrays coverage
ArrayDirectionsFileName = "LinkedCRISPRsDirection.tsv"
ArrayDirections = dict()
for Line in open(ArrayDirectionsFileName):
    LineValues = Line[:-1].split("\t")
    ArrayDirections[LineValues[0]] = LineValues[1]

for ArrayID in ArraysInfoDict:
    for KPlet in range(8, 23):
        if (ArrayID in ArrayHitsInfo) and (KPlet in ArrayHitsInfo[ArrayID]):
        #    Counts = "\t".join([str(len(x)) for x in ArrayHitsInfo[ArrayID][KPlet]])
            ArrayInfo = ArrayID + "\t" + ArraysInfoDict[ArrayID][0] + "\t" + str(KPlet) + "\t" + str(ArraysInfoDict[ArrayID][1])

            for SpacerID in ArrayHitsInfo[ArrayID][KPlet][0]: # hits into genomes
                SpacerPosition = int(SpacerID.split("_")[3])
                if ArrayDirections[ArrayID] == "-":
                    SpacerPosition = ArraysInfoDict[ArrayID][1] - SpacerPosition + 1

                print(ArrayInfo + "\t" + "GenomeHits" + "\t" + str(SpacerPosition)) # spacer no
            for SpacerID in ArrayHitsInfo[ArrayID][KPlet][1]: # hits into viruses
                SpacerPosition = int(SpacerID.split("_")[3])
                if ArrayDirections[ArrayID] == "-":
                    SpacerPosition = ArraysInfoDict[ArrayID][1] - SpacerPosition + 1

                print(ArrayInfo + "\t" + "ViralHits" + "\t" + str(SpacerPosition)) # spacer no
            for SpacerID in ArrayHitsInfo[ArrayID][KPlet][0]:  # random hits positions
                print(ArrayInfo + "\t" + "RandomHits" + "\t" + str(random.randint(1, ArraysInfoDict[ArrayID][1])))  # spacer no

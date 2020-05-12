import os

CRISPRInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_info_known_021718_len.tsv"
ArraysToSpacersQtyDict = dict()
for Line in open(CRISPRInfoFileName):
    LineValues = Line[:-1].split("\t")
    ArraysToSpacersQtyDict["_".join(LineValues[1:4])] = int(LineValues[8])


#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/Linked.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteLinked.tsv"
GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/EnterobacterialesLinked.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelected.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelectedEscherichia.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelectedPseudomonas.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelectedStreptococcus.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelectedBacillus.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelectedSulfolobus.tsv"
#GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelectedMycobacterium.tsv"

ArraysInfoDict = dict()
for Line in open(GenomeToArraysFileName):
    Arrays = Line[:-1].split("\t")[1]
    for Array in Arrays.split(","):
        ArrayID = "_".join(Array.split("_")[:-1])
        Type = Array.split("_")[-1]
        ArraysInfoDict[ArrayID] = [Type, ArraysToSpacersQtyDict[ArrayID]]

# print(len(ArraysDict))
# print(ArraysDict)

#FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"
FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResultsFixed/"
#FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/Hits/"
#FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/PermutatedResultsFixed/"
#FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResultsFixedComplete/"
ArrayHitsInfo = dict()
for FileName in os.listdir(FolderName):
    if not "KPletsCoverageHitsLinked_" in FileName:
        continue

    for Line in open(FolderName + FileName):
        LineValues = Line[:-1].split("\t")
        ArrayID = "_".join(LineValues[0].split("_")[0:3])

        if not ArrayID in ArraysInfoDict:
            break # file not in selected genomes

        if not ArrayID in ArrayHitsInfo:
            ArrayHitsInfo[ArrayID] = dict()

        KPletSize = len(LineValues[1])
        if not KPletSize in ArrayHitsInfo[ArrayID]:
            ArrayHitsInfo[ArrayID][KPletSize] = [set(), set(), set(), set()] # Spacers with hits into genome / viral / rand genome / rand / viral

        #print(Line)
        if len(LineValues) < 4:
            #print(Line)
            continue

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

# ## arrays coverage
# for ArrayID in ArraysInfoDict:
#     for KPlet in range(8, 23):
#         if (ArrayID in ArrayHitsInfo) and (KPlet in ArrayHitsInfo[ArrayID]):
#             Counts = "\t".join([str(len(x)) for x in ArrayHitsInfo[ArrayID][KPlet]])
#         else:
#             Counts = "0\t0\t0\t0"
#
#         print(ArrayID + "\t" + ArraysInfoDict[ArrayID][0] + "\t" + str(KPlet) + "\t" + str(ArraysInfoDict[ArrayID][1]) + "\t" + Counts)

## genomes coverage
for Line in open(GenomeToArraysFileName):
    Arrays = Line[:-1].split("\t")[1]
    for KPlet in range(8, 23):
        GenomeCounts = [0, 0, 0, 0]
        SpacersQty = 0
        for Array in Arrays.split(","):
            ArrayID = "_".join(Array.split("_")[:-1])
            SpacersQty += ArraysInfoDict[ArrayID][1]
            if (ArrayID in ArrayHitsInfo) and (KPlet in ArrayHitsInfo[ArrayID]):
                for I in range(0,4):
                    GenomeCounts[I] += len(ArrayHitsInfo[ArrayID][KPlet][I])

        Counts = "\t".join([str(x) for x in GenomeCounts])
        print(Arrays + "\t" + str(KPlet) + "\t" + str(SpacersQty) + "\t" + Counts)


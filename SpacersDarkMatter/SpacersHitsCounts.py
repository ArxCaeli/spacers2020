import os

# import sys
# sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
# import Helper1603
#
# # CRISPRInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_info_known_021718_len.tsv"
# # ArraysToSpacersQtyDict = dict()
# # for Line in open(CRISPRInfoFileName):
# #     LineValues = Line[:-1].split("\t")
# #     ArraysToSpacersQtyDict["_".join(LineValues[1:4])] = int(LineValues[8])
# #
# #
# # GenomeToArraysFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/Linked.tsv"
# #
# # ArraysInfoDict = dict()
# # for Line in open(GenomeToArraysFileName):
# #     Arrays = Line[:-1].split("\t")[1]
# #     for Array in Arrays.split(","):
# #         SpacerID = "_".join(Array.split("_")[:-1])
# #         Type = Array.split("_")[-1]
# #         ArraysInfoDict[SpacerID] = [Type, ArraysToSpacersQtyDict[SpacerID]]
# #
# # # print(len(ArraysDict))
# # # print(ArraysDict)
# # SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"
# # SpacersFNA = Helper1603.ReadFasta(SpacersFNAFileName)
# #
# # Spacers = set()
# # for SpacerID in SpacersFNA:
# #     if any(CRISPRID in SpacerID for CRISPRID in ArraysInfoDict):
# #         if (len(SpacersFNA[SpacerID]) > 22) and (len(SpacersFNA[SpacerID]) < 50):
# #             Spacers.add(Spacers)


FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"
SpacerHitsInfo = dict()
for FileName in os.listdir(FolderName):
    if not "KPletsCoverageHitsLinked_" in FileName:
        continue

    for Line in open(FolderName + FileName):
        LineValues = Line[:-1].split("\t")
        #ArrayID = "_".join(LineValues[0].split("_")[0:3])
        SpacerID = LineValues[0]

        if not SpacerID in SpacerHitsInfo:
            SpacerHitsInfo[SpacerID] = dict()

        KPletSize = len(LineValues[1])
        if not KPletSize in SpacerHitsInfo[SpacerID]:
            SpacerHitsInfo[SpacerID][KPletSize] = [set(), set(), set(), set()] # Spacers with hits into genome / viral / rand genome / rand / viral

        if LineValues[4] != "":
            if not "_Random" in LineValues[0]:
                SpacerHitsInfo[SpacerID][KPletSize][0].add(LineValues[1])
            else:
                SpacerHitsInfo[SpacerID][KPletSize][2].add(LineValues[1])

        if LineValues[6] != "":
            if not "_Random" in LineValues[0]:
                SpacerHitsInfo[SpacerID][KPletSize][1].add(LineValues[1])
            else:
                SpacerHitsInfo[SpacerID][KPletSize][3].add(LineValues[1])


for SpacerID in SpacerHitsInfo:
    if "_Random" in SpacerID:
        continue
    SpacerSize = int(SpacerID.split("_")[-1])
    #for KPletSize in SpacerHitsInfo[SpacerID]:
    for KPletSize in range(8, 23):
        ResStr = SpacerID + "\t" + str(KPletSize) + "\t"
        KPletCount = SpacerSize - KPletSize + 1

        if KPletSize in SpacerHitsInfo[SpacerID]:
            ResStr += str(len(SpacerHitsInfo[SpacerID][KPletSize][0]) / KPletCount) + "\t" + str(len(SpacerHitsInfo[SpacerID][KPletSize][1]) / KPletCount) + "\t"
        else:
            ResStr += "0\t0\t"

        if KPletSize in SpacerHitsInfo[SpacerID + "_Random"]:
            ResStr += str(len(SpacerHitsInfo[SpacerID + "_Random"][KPletSize][2]) / KPletCount) + "\t" + str(len(SpacerHitsInfo[SpacerID + "_Random"][KPletSize][3]) / KPletCount)
        else:
            ResStr += "0\t0"

        print(ResStr)


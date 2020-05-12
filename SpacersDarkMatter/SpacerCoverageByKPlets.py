import os

def FillKPletInSpacer(SpacerCoverage, StartPositions, KPletSize, SpacerLength, FixedLength, SpacerDirection):
    for Position in StartPositions.split(";"):
        Start = int(Position)
        AdjustedStart = round((Start / SpacerLength) * (FixedLength - 1))
        AdjustedEnd = round(((Start + KPletSize) / SpacerLength) * (FixedLength - 1))
        #print(Position + "\t" + str(Start) + "\t" + str(KPletSize) + "\t" + str(AdjustedStart) + "\t" + str(AdjustedEnd))

        if SpacerDirection == "-":
            tmp = AdjustedStart
            AdjustedStart = FixedLength - AdjustedEnd
            AdjustedEnd = FixedLength - tmp
        for I in range(AdjustedStart, AdjustedEnd):
            SpacerCoverage[I] = 1
    return SpacerCoverage


ArrayDirectionsFileName = "LinkedCRISPRsDirection.tsv"
ArrayDirections = dict()
for Line in open(ArrayDirectionsFileName):
    LineValues = Line[:-1].split("\t")
    ArrayDirections[LineValues[0]] = LineValues[1]


FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"
#FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/PermutatedResults/"
SpacerCoverageDict = dict()
FixedLength = 32

for FileName in os.listdir(FolderName):
    # if not "KPletsCoverageHitsLinked_" in FileName:
    if not "PermutatedKPletsCoverageHitsLinked_" in FileName:
        continue

    #print(FileName)
    for Line in open(FolderName + FileName):
        #CP001891.1_998980_1000411_21_spacer_1000229_32  GCCAGTCCGCATCTCGCCAAAA  1       Genome          Viral   NC_004775.2:34952,-
        LineValues = Line[:-1].split("\t")
        ArrayID = "_".join(LineValues[0].split("_")[0:3])

        KPletSize = len(LineValues[1])

        if len(LineValues[0].split("_")) > 6:
            SpacerLength = int(LineValues[0].split("_")[6])
        else:
            SpacerLength = int(LineValues[0].split("_")[4])

        if not LineValues[0] in SpacerCoverageDict:
            SpacerCoverageDict[LineValues[0]] = dict()

        if not KPletSize in SpacerCoverageDict[LineValues[0]]:
            SpacerCoverageDict[LineValues[0]][KPletSize] = [[0] * FixedLength, [0] * FixedLength, [0] * FixedLength, [0] * FixedLength] # genome / rand genome / vir / randviral

        SpacerDirection = ArrayDirections[ArrayID]
        #print(LineValues[2])
        if LineValues[4] != "":
            if not "_Random" in LineValues[0]:
                SpacerCoverageDict[LineValues[0]][KPletSize][0] = FillKPletInSpacer(SpacerCoverageDict[LineValues[0]][KPletSize][0],
                                                                         LineValues[2], KPletSize, SpacerLength, FixedLength, SpacerDirection)
            else:
                SpacerCoverageDict[LineValues[0]][KPletSize][1] = FillKPletInSpacer(SpacerCoverageDict[LineValues[0]][KPletSize][1],
                                                                         LineValues[2], KPletSize, SpacerLength, FixedLength, SpacerDirection)

        if LineValues[6] != "":
            if not "_Random" in LineValues[0]:
                SpacerCoverageDict[LineValues[0]][KPletSize][2] = FillKPletInSpacer(SpacerCoverageDict[LineValues[0]][KPletSize][2],
                                                                         LineValues[2], KPletSize, SpacerLength, FixedLength, SpacerDirection)
            else:
                SpacerCoverageDict[LineValues[0]][KPletSize][3] = FillKPletInSpacer(SpacerCoverageDict[LineValues[0]][KPletSize][3],
                                                                         LineValues[2], KPletSize, SpacerLength, FixedLength, SpacerDirection)
        # if KPletSize > 8:
        #     break
    #break

# #print(SpacerCoverageDict)
# for SpacerID in SpacerCoverageDict:
#     for KPletSize in SpacerCoverageDict[SpacerID]:
#         #print(KPletSize)
#         for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][0])):
#             if SpacerCoverageDict[SpacerID][KPletSize][0][I] == 1:
#                 print("GenomeHit" + "\t" + str(KPletSize) + "\t" + str(I + 1)) # SpacerID + "\t" +
#         for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][1])):
#             if SpacerCoverageDict[SpacerID][KPletSize][1][I] == 1:
#                 print("RandomGenomeHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))
#         for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][2])):
#             if SpacerCoverageDict[SpacerID][KPletSize][2][I] == 1:
#                 print("ViralHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))
#         for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][3])):
#             if SpacerCoverageDict[SpacerID][KPletSize][3][I] == 1:
#                 print("RandomViralHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))


## to get 1 - 0 difference between true and random
## problem - true can give random hit -???
#print(SpacerCoverageDict)
for SpacerID in SpacerCoverageDict:
    #print(SpacerID)
    if "Random" in SpacerID:
        continue

    for KPletSize in SpacerCoverageDict[SpacerID]:
        #print(KPletSize)
        for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][0])):
            RandomSpacerPositionCovered = False
            if (SpacerID + "_Random" in SpacerCoverageDict) and (KPletSize in SpacerCoverageDict[SpacerID + "_Random"]):
                if SpacerCoverageDict[SpacerID + "_Random"][KPletSize][1][I] == 1:
                    RandomSpacerPositionCovered = True
            if (SpacerCoverageDict[SpacerID][KPletSize][0][I] == 1) and not RandomSpacerPositionCovered:
                print("GenomeHit" + "\t" + str(KPletSize) + "\t" + str(I + 1)) # SpacerID + "\t" +
        # for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][1])):
        #     if SpacerCoverageDict[SpacerID][KPletSize][1][I] == 1:
        #         print("RandomGenomeHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))
        for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][2])):
            RandomSpacerPositionCovered = False

            if (SpacerID + "_Random" in SpacerCoverageDict) and (KPletSize in SpacerCoverageDict[SpacerID + "_Random"]):
                if SpacerCoverageDict[SpacerID + "_Random"][KPletSize][3][I] == 1:
                    RandomSpacerPositionCovered = True

            if (SpacerCoverageDict[SpacerID][KPletSize][2][I] == 1) and not RandomSpacerPositionCovered:
                print("ViralHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))
        # for I in range(0, len(SpacerCoverageDict[SpacerID][KPletSize][3])):
        #     if SpacerCoverageDict[SpacerID][KPletSize][3][I] == 1:
        #         print("RandomViralHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))
    #break
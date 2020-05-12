import os

def AddKPletHits(KPletHits, ContigsHits):
    if ContigsHits == "":
        return KPletHits

    for ContigHits in ContigsHits.split("|"):
        ContigID = ContigHits.split(":")[0]
        Coordinates = ContigHits.split(":")[1].split(";")
        if not ContigID in KPletHits:
            KPletHits[ContigID] = [set(), set()] # +, - hits
        for Coordinate in Coordinates:
            HitInfo = Coordinate.split(",")
            if HitInfo[1] == "+":
                KPletHits[ContigID][0].add(int(HitInfo[0]))
            else:
                KPletHits[ContigID][1].add(int(HitInfo[0]))

    return KPletHits

def GetLocalisationMetric(KPletHits, LookAheadDistance, KPletSize):
    HitDensities = []
    for Contig in KPletHits:
        Coordinates = sorted(KPletHits[Contig][0])
        "bug - need to check - strand as well"
        LastCheckedCoord = 0
        for I in range(0, len(Coordinates)):
            Coord = Coordinates[I]
            if Coord <= LastCheckedCoord:
                continue

            Coveredhits = []
            for J in range(I, len(Coordinates)):
                if Coordinates[J] > Coord + LookAheadDistance - KPletSize:
                    break
                if Coordinates[J] >= Coord and Coordinates[J] <= Coord + LookAheadDistance - KPletSize:
                    Coveredhits.append(Coordinates[J])
            #Coveredhits = [x for x in Coordinates if x >= Coord and x <= Coord + LookAheadDistance - KPletSize]

            if len(Coveredhits) > 0:
                LastCheckedCoord = Coveredhits[-1]
                HitDensities.append(len(Coveredhits))
            else:
                LastCheckedCoord = Coord
                HitDensities.append(1)

    return HitDensities


## hits in contig not separated by spacers
# FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"
#
# LookAheadDistance = 32
#
# for FileName in os.listdir(FolderName):
#     if not "KPletsCoverageHitsLinked_" in FileName:
#         continue
#
#     Hits = dict()
#     ResFileName = FileName.split(".")[0] + "_densities.tsv"
#     print(ResFileName)
#     with open(ResFileName, "w") as ResFile:
#         for Line in open(FolderName + FileName):
#             #CP001891.1_998980_1000411_21_spacer_1000229_32  GCCAGTCCGCATCTCGCCAAAA  1       Genome          Viral   NC_004775.2:34952,-
#             LineValues = Line[:-1].split("\t")
#             ArrayID = "_".join(LineValues[0].split("_")[0:3])
#
#             KPletSize = len(LineValues[1])
#
#             if not KPletSize in Hits:
#                 Hits[KPletSize] = [dict(), dict(), dict(), dict()] # Spacers with hits into genome / rand genome / viral / rand viral
#
#             if LineValues[4] != "":
#                 if not "_Random" in LineValues[0]:
#                     Hits[KPletSize][0] = AddKPletHits(Hits[KPletSize][0], LineValues[4])
#                 else:
#                     Hits[KPletSize][1] = AddKPletHits(Hits[KPletSize][1], LineValues[4])
#
#             if LineValues[6] != "":
#                 if not "_Random" in LineValues[0]:
#                     Hits[KPletSize][2] = AddKPletHits(Hits[KPletSize][2], LineValues[6])
#                 else:
#                     Hits[KPletSize][3] = AddKPletHits(Hits[KPletSize][3], LineValues[6])
#             #break
#         for KPletSize in Hits:
#             HitDensities = GetLocalisationMetric(Hits[KPletSize][0], LookAheadDistance, KPletSize)
#             for Hit in HitDensities:
#                 ResFile.write("Genome" + "\t" + str(KPletSize) + "\t" + str(LookAheadDistance) + "\t" + str(Hit) + "\n")
#             HitDensities = GetLocalisationMetric(Hits[KPletSize][1], LookAheadDistance, KPletSize)
#             for Hit in HitDensities:
#                 ResFile.write("RandomGenome" + "\t" + str(KPletSize) + "\t" + str(LookAheadDistance) + "\t" + str(Hit) + "\n")
#             HitDensities = GetLocalisationMetric(Hits[KPletSize][2], LookAheadDistance, KPletSize)
#             for Hit in HitDensities:
#                 ResFile.write("Viral" + "\t" + str(KPletSize) + "\t" + str(LookAheadDistance) + "\t" + str(Hit) + "\n")
#             HitDensities = GetLocalisationMetric(Hits[KPletSize][3], LookAheadDistance, KPletSize)
#             for Hit in HitDensities:
#                 ResFile.write("RandomViral" + "\t" + str(KPletSize) + "\t" + str(LookAheadDistance) + "\t" + str(Hit) + "\n")
#
#     #break


#FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"
FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResultsFixed/"

#LookAheadDistance = 32

for FileName in os.listdir(FolderName):
    if not "KPletsCoverageHitsLinked_" in FileName:
        continue

    Hits = dict()
    ResFileName = FileName.split(".")[0] + "_densities.tsv"
    print(ResFileName)
    with open(ResFileName, "w") as ResFile:
        for Line in open(FolderName + FileName):
            #CP001891.1_998980_1000411_21_spacer_1000229_32  GCCAGTCCGCATCTCGCCAAAA  1       Genome          Viral   NC_004775.2:34952,-
            LineValues = Line[:-1].split("\t")
            ArrayID = "_".join(LineValues[0].split("_")[0:3])

            KPletSize = len(LineValues[1])

            if not KPletSize in Hits:
                Hits[KPletSize] = dict() # for spacers

            if not LineValues[0] in Hits[KPletSize]:
                Hits[KPletSize][LineValues[0]] = [dict(), dict()] # Spacers with hits into genome / viral

            if LineValues[4] != "":
                Hits[KPletSize][LineValues[0]][0] = AddKPletHits(Hits[KPletSize][LineValues[0]][0], LineValues[4])

            if LineValues[6] != "":
                Hits[KPletSize][LineValues[0]][1] = AddKPletHits(Hits[KPletSize][LineValues[0]][1], LineValues[6])

            #break
        for KPletSize in Hits:
            for SpacerID in Hits[KPletSize]:
                HitTypePrefix = ""
                if "_Random" in SpacerID:
                    HitTypePrefix = "Random"
                    LookAheadDistance = int(SpacerID.split("_")[-2])
                else:
                    LookAheadDistance = int(SpacerID.split("_")[-1]) # spacer length
                HitDensities = GetLocalisationMetric(Hits[KPletSize][SpacerID][0], LookAheadDistance, KPletSize)
                for Hit in HitDensities:
                    ResFile.write(HitTypePrefix + "Genome" + "\t" + str(KPletSize) + "\t" + str(LookAheadDistance) + "\t" + str(Hit) + "\n")
                HitDensities = GetLocalisationMetric(Hits[KPletSize][SpacerID][1], LookAheadDistance, KPletSize)
                for Hit in HitDensities:
                    ResFile.write(HitTypePrefix + "Viral" + "\t" + str(KPletSize) + "\t" + str(LookAheadDistance) + "\t" + str(Hit) + "\n")

    #break


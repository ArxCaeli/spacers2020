# def GetCRISPRID(SpacerId):
#     return "_".join(SpacerId.split("_")[0:3])
#
#
# GenomeArrays = dict()
# for Line in open("LinkedResultsFixed/KPletsCoverageLinked_50.tsv"):
#     LineValues = Line[:-1].split("\t")
#     GenomeArrays[LineValues[2]] = dict()
#     for I in range(8, 23):
#         GenomeArrays[LineValues[2]][I] = [set(), set(), set(), set()]
#
# for Line in open("LinkedResultsFixed/KPletsCoverageHitsLinked_50.tsv"):
#     LineValues = Line[:-1].split("\t")
#
#     for ArraySet in GenomeArrays:
#         CRISPRID = GetCRISPRID(LineValues[0])
#
#         if CRISPRID in ArraySet:
#             if not "_Random" in LineValues[0]:
#                 if LineValues[4] != "":
#                     GenomeArrays[ArraySet][len(LineValues[1])][0].add(LineValues[0])
#                 if LineValues[6] != "":
#                     GenomeArrays[ArraySet][len(LineValues[1])][1].add(LineValues[0])
#             else:
#                 if LineValues[4] != "":
#                     GenomeArrays[ArraySet][len(LineValues[1])][2].add(LineValues[0])
#                 if LineValues[6] != "":
#                     GenomeArrays[ArraySet][len(LineValues[1])][3].add(LineValues[0])


## get experimental spacers coverage
import os

# import sys
# sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
# import Helper1603
#
# #SpacersFasta = Helper1603.ReadFasta("Data/T5-spacers_uniq.fasta")
# SpacersFasta = Helper1603.ReadFasta("Data/lambda-spacers_uniq.fasta")

## to get raw numbers for kmers in initial set of spacers
import random

Spacers = dict()
#for Line in open("Data/T5new-spacers.fasta"):
for Line in open("Data/lambda-spacers.fasta"):
    if Line[0] == ">":
        continue

    if Line[:-1] in Spacers:
        Spacers[Line[:-1]] += 1
    else:
        Spacers[Line[:-1]] = 1

SpacerCounts = dict()
#for Line in open("Data/T5-spacers_uniq.fasta"):
for Line in open("Data/lambda-spacers_uniq.fasta"):
    if Line[0] == ">":
        SpacerID = Line[1:-1]
    else:
        SpacerCounts[SpacerID] = Spacers[Line[:-1]]

HitsFolder = "HitsFixed/"
ArrayHitsInfo = dict()
#ArrayID = "T5Spacers"
ArrayID = "LaSpacers"
ArrayHitsInfo[ArrayID] = dict()

for FileName in os.listdir(HitsFolder):
    #if not "KPletsCoverageHitsLinked_T5_uniq_" in FileName:
    if not "KPletsCoverageHitsLinked_La_uniq_" in FileName:
        continue

    KPletSize = int(FileName.split("_")[-1].split(".")[0])
    # if KPletSize > 22:
    #     continue
    if KPletSize != 22:
        continue

    # if KPletSize != 16:
    #     continue

    if not KPletSize in ArrayHitsInfo[ArrayID]:
        ArrayHitsInfo[ArrayID][KPletSize] = [set(), set(), set(), set()] # Spacers with hits into genome / viral / rand genome / rand / viral

    for Line in open(HitsFolder + FileName):
        LineValues = Line[:-1].split("\t")


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

# for SpacerID in ArrayHitsInfo[ArrayID][16][1]:
#     print(SpacerID)

# res files
# LaSpacersCoverage.tsv
# T5SpacersCoverage.tsv

## genomes coverage
# SpacersQty = 0
# for Spacer in Spacers:
#     SpacersQty += Spacers[Spacer]
SpacersQty = 0
for Spacer in SpacerCounts:
    SpacersQty += SpacerCounts[Spacer]

#for KPlet in range(8, 23):
for KPlet in range(22, 23):
    GenomeCounts = [0, 0, 0, 0]
    #SpacersQty = 7266 #len(SpacersFasta)

    #if (ArrayID in ArrayHitsInfo) and (KPlet in ArrayHitsInfo[ArrayID]):
    #for I in range(0,4):
    #    GenomeCounts[I] += sum([SpacerCounts[x] for x in ArrayHitsInfo[ArrayID][KPlet][I]])
    GenomeCounts[0] += sum([SpacerCounts[x] for x in ArrayHitsInfo[ArrayID][KPlet][0]])
    GenomeCounts[1] += sum([SpacerCounts[x] for x in ArrayHitsInfo[ArrayID][KPlet][1]])
    # GenomeCounts[2] += sum([SpacerCounts["_".join(x.split("_")[:-1])] for x in ArrayHitsInfo[ArrayID][KPlet][2]])
    # GenomeCounts[3] += sum([SpacerCounts["_".join(x.split("_")[:-1])] for x in ArrayHitsInfo[ArrayID][KPlet][3]])
    # GenomeCounts[2] += len(ArrayHitsInfo[ArrayID][KPlet][2])
    # GenomeCounts[3] += len(ArrayHitsInfo[ArrayID][KPlet][3])
    BootstrapQty = len(SpacerCounts)
    for I in range(0, BootstrapQty):
        GenomeCounts[2] += sum(random.sample(list(SpacerCounts.values()), len(ArrayHitsInfo[ArrayID][KPlet][2]) + 1))
        GenomeCounts[3] += sum(random.sample(list(SpacerCounts.values()), len(ArrayHitsInfo[ArrayID][KPlet][3]) + 1))
    GenomeCounts[2] = GenomeCounts[2] / BootstrapQty
    GenomeCounts[3] = GenomeCounts[3] / BootstrapQty


    Counts = "\t".join([str(x) for x in GenomeCounts])

    print(ArrayID + "\t" + str(KPlet) + "\t" + str(SpacersQty) + "\t" + Counts)


# import random
# x = {"a":1, "b":2, "c":3, "d":4, "e":5, "f":6}
# list(x.values())
# random.sample(list(x.values()), 5)
#len(x)


# ## get number of shared spacers per K-mer
# import os
#
#
# HitsFolder = "HitsFixed/"
# ArrayHitsInfo = dict()
# ArrayID = "T5Spacers"
# #ArrayID = "LaSpacers"
# ArrayHitsInfo[ArrayID] = dict()
#
# KmerToSpacer = dict()
# for FileName in os.listdir(HitsFolder):
#     if not "KPletsCoverageHitsLinked_T5_uniq_" in FileName:
#     #if not "KPletsCoverageHitsLinked_La_uniq_" in FileName:
#         continue
#
#     KPletSize = int(FileName.split("_")[-1].split(".")[0])
#     if KPletSize > 22:
#         continue
#
#     # if KPletSize != 16:
#     #     continue
#
#     if not KPletSize in KmerToSpacer:
#         KmerToSpacer[KPletSize] = dict()
#
#     if not KPletSize in ArrayHitsInfo[ArrayID]:
#         ArrayHitsInfo[ArrayID][KPletSize] = [set(), set(), set(), set()] # Spacers with hits into genome / viral / rand genome / rand / viral
#
#
#     for Line in open(HitsFolder + FileName):
#         LineValues = Line[:-1].split("\t")
#
#         if LineValues[4] != "":
#             if not "_Random" in LineValues[0]:
#                 if not LineValues[1] in KmerToSpacer[KPletSize]:
#                     KmerToSpacer[KPletSize][LineValues[1]] = set()
#
#                 KmerToSpacer[KPletSize][LineValues[1]].add(LineValues[0])
#
#                 ArrayHitsInfo[ArrayID][KPletSize][0].add(LineValues[0])
#             else:
#                 ArrayHitsInfo[ArrayID][KPletSize][2].add(LineValues[0])
#
#         if LineValues[6] != "":
#             if not "_Random" in LineValues[0]:
#                 # if not LineValues[1] in KmerToSpacer[KPletSize]:
#                 #     KmerToSpacer[KPletSize][LineValues[1]] = set()
#                 #
#                 # KmerToSpacer[KPletSize][LineValues[1]].add(LineValues[0])
#
#                 ArrayHitsInfo[ArrayID][KPletSize][1].add(LineValues[0])
#             else:
#                 ArrayHitsInfo[ArrayID][KPletSize][3].add(LineValues[0])
#
# for KPlet in range(8, 23):
#     for Kmer in KmerToSpacer[KPlet]:
#         print("Org" + "\t" + Kmer + "\t" + str(KPlet) + "\t" + str(len(KmerToSpacer[KPlet][Kmer])))

# ## to get raw numbers for kmers in initial set of spacers
# Spacers = dict()
# for Line in open("Data/T5new-spacers.fasta"):
#     if Line[:-1] in Spacers:
#         Spacers[Line[:-1]] += 1
#     else:
#         Spacers[Line[:-1]] = 1
#
# for Spacer in Spacers:
#     print(Spacer + "\t" + str(Spacers[Spacer]))
# # #for Line in open("KmerPerSpacer_viral_16.tsv"):
# # for Line in open("KmerPerSpacer_viral_22.tsv"):
# #     LineValues = Line[:-1].split("\t")
# #
# #     SpacersCount = 0
# #     for Spacer in Spacers:
# #         if LineValues[1] in Spacer:
# #             SpacersCount += Spacers[Spacer]
# #
# #     print("\t".join(LineValues) + "\t" + str(SpacersCount))

#
# ## to get Hits per spacer
# import os
#
# import sys
# sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
# import Helper1603
#
# def AddSeeds(SeedsDict, HitsText):
#     for ContigHits in HitsText.split("|"):
#         ContigID = ContigHits.split(":")[0]
#         if not ContigID in SeedsDict:
#             SeedsDict[ContigID] = []
#         Hits = ContigHits.split(":")[1]
#
#         for HitCoord in Hits.split(";"):
#             HitStart = int(HitCoord.split(",")[0])
#             if HitCoord.split(",")[1] == -1:
#                 HitStart += pow(10, 10)
#
#             if not Helper1603.IsInSeeds(ContigID, HitStart - 100, HitStart + 200, SeedsDict):
#                 SeedsDict[ContigID].append([HitStart, HitStart + 1])
#
# def PrintSeeds(Type, KmerSize, Seeds):
#     #print(Seeds)
#     for SpacerID in Seeds:
#         #print(Seeds[SpacerID])
#
#         for ContigID in Seeds[SpacerID]:
#             SpacerType = "CRISPR"
#             if "_Random" in SpacerID:
#                 SpacerType = "Random"
#             print(Type + "\t" + SpacerType + "\t" + str(KmerSize) + "\t" + SpacerID + "\t" + str(len(Seeds[SpacerID][ContigID])))
#
#
# HitsFolder = "LinkedResultsFixed/"
#
# Seeds = dict()
# ViralSeeds = dict()
#
# KmerSize = 0
# for I in range(1, 155):
#     for Line in open(HitsFolder + "/" + "KPletsCoverageHitsLinked_" + str(I) + ".tsv"):
#         LineValues = Line[:-1].split("\t")
#
#         if len(LineValues[1]) != KmerSize:
#             if KmerSize != 0:
#                 PrintSeeds("Self", KmerSize, Seeds)
#                 PrintSeeds("Viral", KmerSize, ViralSeeds)
#
#             KmerSize = len(LineValues[1])
#             Seeds = dict()
#             ViralSeeds = dict()
#
#         # if KmerSize != 22:
#         #     continue
#
#         if LineValues[4] != "":
#             if not LineValues[0] in Seeds:
#                 Seeds[LineValues[0]] = dict()
#             AddSeeds(Seeds[LineValues[0]], LineValues[4])
#
#         if LineValues[6] != "":
#             if not LineValues[0] in ViralSeeds:
#                 ViralSeeds[LineValues[0]] = dict()
#             AddSeeds(ViralSeeds[LineValues[0]], LineValues[6])
#
# PrintSeeds("Self", KmerSize, Seeds)
# PrintSeeds("Viral", KmerSize, ViralSeeds)


# ## To get no of kmers per spacer
# import os
#
# import sys
# sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
# import Helper1603
#
# Spacers = dict()
# #for Line in open("Data/T5new-spacers.fasta"):
# for Line in open("Data/lambda-spacers.fasta"):
#     if Line[0] == ">":
#         continue
#
#     if Line[:-1] in Spacers:
#         Spacers[Line[:-1]] += 1
#     else:
#         Spacers[Line[:-1]] = 1
#
# SpacerCounts = dict()
# TotalCount = 0
# #for Line in open("Data/T5-spacers_uniq.fasta"):
# for Line in open("Data/lambda-spacers_uniq.fasta"):
#     if Line[0] == ">":
#         SpacerID = Line[1:-1]
#     else:
#         SpacerCounts[SpacerID] = Spacers[Line[:-1]]
#         TotalCount += Spacers[Line[:-1]]
#
# #SpacersFasta = Helper1603.ReadFasta("Data/T5-spacers_uniq.fasta")
# SpacersFasta = Helper1603.ReadFasta("Data/lambda-spacers_uniq.fasta")
#
# HitsFolder = "HitsFixed/"
# ArrayHitsInfo = dict()
# #ArrayID = "T5Spacers"
# ArrayID = "LaSpacers"
# ArrayHitsInfo[ArrayID] = dict()
#
# for FileName in os.listdir(HitsFolder):
#     #if not "KPletsCoverageHitsLinked_T5_uniq_" in FileName:
#     if not "KPletsCoverageHitsLinked_La_uniq_" in FileName:
#         continue
#
#     KPletSize = int(FileName.split("_")[-1].split(".")[0])
#     if KPletSize > 22:
#         continue
#
#     # if KPletSize != 16:
#     #     continue
#
#     if not KPletSize in ArrayHitsInfo[ArrayID]:
#         ArrayHitsInfo[ArrayID][KPletSize] = [dict(), dict(), dict(), dict()] # Spacers with hits into genome / viral / rand genome / rand / viral
#
#     for Line in open(HitsFolder + FileName):
#         LineValues = Line[:-1].split("\t")
#
#         #!!!
#         if "_Random" in LineValues[0]:
#             continue
#
#         if LineValues[4] != "":
#             if not "_Random" in LineValues[0]:
#                 if not LineValues[0] in ArrayHitsInfo[ArrayID][KPletSize][0]:
#                     ArrayHitsInfo[ArrayID][KPletSize][0][LineValues[0]] = 0
#                 ArrayHitsInfo[ArrayID][KPletSize][0][LineValues[0]] += SpacersFasta[LineValues[0]].count(LineValues[1])
#             else:
#                 if not LineValues[0] in ArrayHitsInfo[ArrayID][KPletSize][2]:
#                     ArrayHitsInfo[ArrayID][KPletSize][2][LineValues[0]] = 0
#                 ArrayHitsInfo[ArrayID][KPletSize][2][LineValues[0]] += 1
#
#         if LineValues[6] != "":
#             if not "_Random" in LineValues[0]:
#                 if not LineValues[0] in ArrayHitsInfo[ArrayID][KPletSize][1]:
#                     ArrayHitsInfo[ArrayID][KPletSize][1][LineValues[0]] = 0
#                 ArrayHitsInfo[ArrayID][KPletSize][1][LineValues[0]] += SpacersFasta[LineValues[0]].count(LineValues[1])
#             else:
#                 if not LineValues[0] in ArrayHitsInfo[ArrayID][KPletSize][3]:
#                     ArrayHitsInfo[ArrayID][KPletSize][3][LineValues[0]] = 0
#                 ArrayHitsInfo[ArrayID][KPletSize][3][LineValues[0]] += 1
#
# for KPlet in range(8, 23):
#     for SpacerID in SpacersFasta:
#         GenomeCounts = [0, 0, 0, 0]
#
#         for I in range(0,4):
#             if SpacerID in ArrayHitsInfo[ArrayID][KPlet][I]:
#                 GenomeCounts[I] = ArrayHitsInfo[ArrayID][KPlet][I][SpacerID]
#
#         Counts = "\t".join([str(x) for x in GenomeCounts])
#
#         print(ArrayID + "\t" + str(KPlet) + "\t" + str(len(SpacersFasta[SpacerID])) + "\t" + Counts + "\t" + str(SpacerCounts[SpacerID] / TotalCount))
#

# ## to find hits for a given spacer
# import os
# SpacerID = "AE006641.1_1305634_1311637_52_spacer_1308869_40_Random"
# HitsPath = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResultsFixed/"
# IDFound = False
# MaxHitsKmer = [0, ""]
# for FileName in os.listdir(HitsPath):
#     for Line in open(HitsPath + FileName):
#         LineValues = Line[:-1].split("\t")
#         if LineValues[0] != SpacerID:
#             continue
#         if len(LineValues[1]) != 22:
#             continue
#         IDFound = True
#         HitsCount = len(LineValues[4].split(";"))
#
#         if HitsCount > MaxHitsKmer[0]:
#             MaxHitsKmer = [HitsCount, Line[:-1]]
#
#     if IDFound:
#         print(MaxHitsKmer[1])
#
#         break


# ## to find hits for a given spacer set
# import os
# SpacerIDs = []
# for Line in open("tmp_spacers.ids"):
#     SpacerIDs.append(Line[:-1])
#
#
# #SpacerID = "AE006641.1_1305634_1311637_52_spacer_1308869_40_Random"
# HitsPath = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResultsFixed/"
# ContigID = SpacerIDs[0].split("_")[0]
# GenomeID = ""
# for FileName in os.listdir(HitsPath):
#     if not "KPletsCoverageLinked" in FileName:
#         continue
#     for Line in open(HitsPath + FileName):
#         if ContigID in Line:
#             GenomeID = FileName.split("_")[1].split(".")[0]
#         break
#     if GenomeID != "":
#         break
#
# HitsFileName = HitsPath + "KPletsCoverageHitsLinked_" + GenomeID + ".tsv"
# for SpacerID in SpacerIDs:
#     IDFound = False
#     MaxHitsKmer = [0, ""]
#     for Line in open(HitsFileName):
#         LineValues = Line[:-1].split("\t")
#         if LineValues[0] != SpacerID:
#             continue
#         if len(LineValues[1]) != 22:
#             continue
#         IDFound = True
#         HitsCount = len(LineValues[4].split(";"))
#
#         if HitsCount > MaxHitsKmer[0]:
#             MaxHitsKmer = [HitsCount, Line[:-1]]
#
#     if IDFound:
#         print(MaxHitsKmer[1])


# ## to get organism names for contig id list
# ContigIDToName = dict()
#
# for Line in open("/panfs/pan1/prokdata/db/Prok1603/Prok1603.pp.txt"):
#     LineValues = Line[:-1].split("\t")
#
#     ContigIDToName[LineValues[1]] = LineValues[0]
#
# for Line in open("TopRandomHitsContigIDs.lst"):
#     print(Line[:-1] + "\t" + ContigIDToName[Line[:-1]])
import sys
sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
import Helper1603

import subprocess
import re
import random
import os
import argparse

Alphabet = {"A", "T", "G", "C"}

def GetKPliets(Sequence, FrameLength):
    KPlets = set()

    for I in range(FrameLength, len(Sequence)):
        if any(x not in Alphabet for x in Sequence[I - FrameLength : I]):
            continue
        KPlets.add(Sequence[I - FrameLength : I])
        #KPlets.add(Helper1603.ReverseComplement(Sequence[I - FrameLength : I]))

    return KPlets

def GetGenomeFNA(GenomeID, CompleteGenomesDict):
    TmpFileName = GenomeID + "_Complete_Tmp.tmp" + str(KPletSize)
    TmpFastaFileName = GenomeID + "_Complete_Fasta.tmp" + str(KPletSize)
    DBFileName = "/panfs/pan1/prokdata/db/Prok1603/Prok1603.nt"

    # find kplets in genomes
    with open(TmpFileName, "w") as TmpFile:
        for Contig in CompleteGenomesDict[GenomeID]:
            TmpFile.write(Contig + "\n")

    subprocess.call("blastdbcmd -db " + DBFileName + " -entry_batch " + TmpFileName + " > " + TmpFastaFileName, shell=True)
    GenomeFNA = Helper1603.ReadFasta(TmpFastaFileName)
    TotalLength = 0
    for ID in GenomeFNA:
        TotalLength += len(GenomeFNA[ID])

    os.remove(TmpFileName)
    os.remove(TmpFastaFileName)

    return GenomeFNA, TotalLength

def MakeRandomKPletsFromGenomeFNA(GenomeFNA, GenomeSpacersFNA, KPletSize):
    RandomKPlets = set()
    Seeds = dict()
    RandomSpacersFNA = dict()

    GenomeFNASequence = ""
    for GenomeID in GenomeFNA:
        GenomeFNASequence += GenomeFNA[GenomeID]

    #RandomSpacers = []
    for SpacerID in GenomeSpacersFNA:
        RandomSpacerPosition = random.randrange(1, len(GenomeFNASequence) - len(GenomeSpacersFNA[SpacerID])) - 1
        #RandomSpacers.append(GenomeFNASequence[RandomSpacerPosition: RandomSpacerPosition + len(GenomeSpacersFNA[SpacerID])])
        RandomSpacersFNA[SpacerID + "_Random"] = GenomeFNASequence[RandomSpacerPosition: RandomSpacerPosition + len(GenomeSpacersFNA[SpacerID])]

    for SpacerID in RandomSpacersFNA:
        RandomKPlets = RandomKPlets.union(GetKPliets(RandomSpacersFNA[SpacerID], KPletSize))

    for SpacerID in RandomSpacersFNA:
        for GenomeID in GenomeFNA:
            Position = GenomeFNA[GenomeID].find(RandomSpacersFNA[SpacerID])
            GenomeIDClear = GenomeID.split(" ")[0]
            if Position != -1:
                if not GenomeIDClear in Seeds:
                    Seeds[GenomeIDClear] = []
                Seeds[GenomeIDClear].append([Position, Position + len(RandomSpacersFNA[SpacerID])])
                break

    return RandomKPlets, Seeds, RandomSpacersFNA

def UpdateSpacerHits(KPlet, SpacersFNA, HitStart, Strand, SpacersHits, Type, ContigID):
    for SpacerID in SpacersFNA:
        if KPlet in SpacersFNA[SpacerID]:
            if not SpacerID in SpacersHits:
                SpacersHits[SpacerID] = dict()
            if not KPlet in SpacersHits[SpacerID]:
                SpacersHits[SpacerID][KPlet] = [[], dict()] ## spacer hits / hits

                for SubstringStart in [m.start() for m in re.finditer(KPlet, SpacersFNA[SpacerID])]: # multiple kplets
                    SpacersHits[SpacerID][KPlet][0].append([SubstringStart])

            ## add genome/viral hit
            if not Type in SpacersHits[SpacerID][KPlet][1]:
                SpacersHits[SpacerID][KPlet][1][Type] = dict()
            if not ContigID in SpacersHits[SpacerID][KPlet][1][Type]:
                SpacersHits[SpacerID][KPlet][1][Type][ContigID] = []

            SpacersHits[SpacerID][KPlet][1][Type][ContigID].append([HitStart, Strand])

    return SpacersHits

def GetKPletCoverage(SpacerKPlets, GenomeFNA, CRISPRSeeds, SpacersFNA):
    # get kplets in spacers
    SpacersKPletsDict = dict()
    SpacersHits = dict()
    for KPlet in SpacerKPlets:
        SpacersKPletsDict[KPlet] = [0, 0]

    for ContigID in GenomeFNA:
        for KPlet in SpacersKPletsDict:
            for SubstringStart in [m.start() for m in re.finditer(KPlet, GenomeFNA[ContigID])]:
                if not Helper1603.IsInSeeds(ContigID.split(" ")[0], SubstringStart - 60, SubstringStart + 60, CRISPRSeeds):
                    SpacersKPletsDict[KPlet][0] += 1
                    SpacersHits = UpdateSpacerHits(KPlet, SpacersFNA, SubstringStart, "+", SpacersHits, "Genome", ContigID)

            for SubstringStart in [m.start() for m in re.finditer(Helper1603.ReverseComplement(KPlet), GenomeFNA[ContigID])]:
                if not Helper1603.IsInSeeds(ContigID.split(" ")[0], SubstringStart - 60, SubstringStart + 60, CRISPRSeeds):
                    SpacersKPletsDict[KPlet][0] += 1
                    SpacersHits = UpdateSpacerHits(KPlet, SpacersFNA, SubstringStart, "-", SpacersHits, "Genome", ContigID)

    # print output
    GenomeHits = 0
    for KPlet in SpacersKPletsDict:
        if SpacersKPletsDict[KPlet][0] > 0:
            GenomeHits += 1

    return GenomeHits, SpacersHits

def WriteSpacerHits(SpacerKPletHitsFile, SpacerHits):
    for SpacerID in SpacerHits:
        for KPlet in SpacerHits[SpacerID]:
            SpacerKPletHitsFile.write(SpacerID + "\t" + KPlet + "\t" +
                                      ";".join([str(SpacerKPletCoord[0]) for SpacerKPletCoord in
                                                SpacerHits[SpacerID][KPlet][0]]) + "\t")
            SpacerKPletHitsFile.write("Genome\t")
            TypeHits = ""
            if "Genome" in SpacerHits[SpacerID][KPlet][1]:
                if len(SpacerHits[SpacerID][KPlet][1]["Genome"]) > 0:
                    for ContigID in SpacerHits[SpacerID][KPlet][1]["Genome"]:
                        if len(SpacerHits[SpacerID][KPlet][1]["Genome"][ContigID]) > 0:
                            TypeHits += ContigID.split(" ")[0] + ":" + ";".join([(str(GenomesKPletCoord[0]) + "," + str(GenomesKPletCoord[1])) for GenomesKPletCoord in SpacerHits[SpacerID][KPlet][1]["Genome"][ContigID]]) + "|"
                    TypeHits = TypeHits[:-1]
            SpacerKPletHitsFile.write(TypeHits + "\n")


CompleteGenomesSetFileName = "/panfs/pan1/prokdata/db/Prok1603/Prok1603.pp.txt"
CompleteGenomesDict = dict()
for Line in open(CompleteGenomesSetFileName, "r"):
    if Line[0] == "#":
        continue
    LineValues = Line[:-1].split("\t")
    if not LineValues[0] in CompleteGenomesDict:
        CompleteGenomesDict[LineValues[0]] = []
    CompleteGenomesDict[LineValues[0]].append(LineValues[1])


ap = argparse.ArgumentParser()
ap.add_argument("-n", help = "LineNo", required = True)
opts = ap.parse_args()

TargetLineNo = int(opts.n)


LinkedInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteLinked.tsv"
SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"

SpacersFNA = Helper1603.ReadFasta(SpacersFNAFileName)

Counter = 0
for LinkedLine in open(LinkedInfoFileName):
    Counter += 1
    if Counter == TargetLineNo:
        break

ResFileName = "KPletsCoverageComplete_" + str(TargetLineNo) + ".tsv"
SpacerKPletHitsFileName = "KPletsCoverageHitsComplete_" + str(TargetLineNo) + ".tsv"

with open(ResFileName, "w") as ResFile:
    with open(SpacerKPletHitsFileName, "w") as SpacerKPletHitsFile:
        for KPletSize in range(8, 23):
        #for KPletSize in range(22, 23):
        #for KPletSize in range(8, 9):
            #Counter = 0
            LinkedLineValues = LinkedLine[:-1].split("\t")
            GenomeID = LinkedLineValues[0]
            CRISPRs = LinkedLineValues[1].split(",")

            CRISPRSeeds = dict()
            for CRISPR in CRISPRs:
                CRISPRValues = CRISPR.split("_")
                if not CRISPRValues[0] in CRISPRSeeds:
                    CRISPRSeeds[CRISPRValues[0]] = []
                CRISPRSeeds[CRISPRValues[0]].append([int(CRISPRValues[1]), int(CRISPRValues[2])])

            CRISPRIDs = []
            for CRISPR in CRISPRs:
                CRISPRIDs.append("_".join(CRISPR.split("_")[0:3]))

            SpacerKPlets = set()
            Spacers = []
            GenomeSpacersFNA = dict()
            for SpacerID in SpacersFNA:
                if any(CRISPRID in SpacerID for CRISPRID in CRISPRIDs):
                    if (len(SpacersFNA[SpacerID]) > 22) and (len(SpacersFNA[SpacerID]) < 50):
                        SpacerKPlets = SpacerKPlets.union(GetKPliets(SpacersFNA[SpacerID], KPletSize))
                        GenomeSpacersFNA[SpacerID] = SpacersFNA[SpacerID]

            GenomeFNA, GenomeLength = GetGenomeFNA(GenomeID, CompleteGenomesDict)

            GenomeHits, SpacerHits = GetKPletCoverage(SpacerKPlets, GenomeFNA, CRISPRSeeds, GenomeSpacersFNA)

            RandomKPlets, RandomSeeds, RandomSpacersFNA = MakeRandomKPletsFromGenomeFNA(GenomeFNA, GenomeSpacersFNA, KPletSize)
            RandomGenomeHits, RandomSpacerHits = GetKPletCoverage(RandomKPlets, GenomeFNA, RandomSeeds, RandomSpacersFNA)

            ResFile.write(str(KPletSize) + "\t" + LinkedLine[:-1] + "\t" +
                          str(GenomeLength) + "\t" +
                          str(len(SpacerKPlets)) + "\t" + str(GenomeHits) + "\t" +
                          str(len(RandomKPlets)) +"\t" + str(RandomGenomeHits) + "\n")


            # SpacersHits[SpacerID][KPlet][1][Type][ContigID]
            WriteSpacerHits(SpacerKPletHitsFile, SpacerHits)
            WriteSpacerHits(SpacerKPletHitsFile, RandomSpacerHits)

#break
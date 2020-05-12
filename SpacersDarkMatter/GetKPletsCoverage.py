import sys
sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
import Helper1603

import subprocess
import re
import random
import os
import argparse

Alphabet = {"A", "T", "G", "C"}

def GetKPlets(Sequence, FrameLength):
    KPlets = set()

    for I in range(FrameLength, len(Sequence) + 1):
        if any(x not in Alphabet for x in Sequence[I - FrameLength : I]):
            continue
        KPlets.add(Sequence[I - FrameLength : I])
        #KPlets.add(Helper1603.ReverseComplement(Sequence[I - FrameLength : I]))

    return KPlets

def GetGenomeFNA(GenomeID, CompleteGenomesDict):
    TmpFileName = GenomeID + "Tmp.tmp"
    TmpFastaFileName = GenomeID + "Fasta.tmp"
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

def GetViralFNA(GenomeID, VirusContigs):
    TmpFileName = GenomeID + "Tmp.tmp"
    TmpFastaFileName = GenomeID + "Fasta.tmp"

    # with open(TmpFileName, "w") as TmpFile:
    #     for Contig in VirusContigs:
    #         TmpFile.write(Contig + "\n")
    #subprocess.call("blastdbcmd -db nt -entry_batch " + TmpFileName + " > " + TmpFastaFileName, shell=True)

    for Contig in VirusContigs:
        if not os.path.exists(TmpFastaFileName):
            subprocess.call("efetch -db nucleotide -id " + Contig + " -format fasta > " + TmpFastaFileName, shell=True)
        else:
            subprocess.call("efetch -db nucleotide -id " + Contig + " -format fasta >> " + TmpFastaFileName, shell=True)

    ViralFNA = Helper1603.ReadFasta(TmpFastaFileName)

    TotalLength = 0
    for ID in ViralFNA:
        TotalLength += len(ViralFNA[ID])

    #os.remove(TmpFileName)
    os.remove(TmpFastaFileName)

    return ViralFNA, TotalLength

def MakeRandomKPletsFromGenomeFNA(GenomeFNA, GenomeSpacersFNA):
    Seeds = dict()
    RandomSpacersFNA = dict()

    GenomeFNASequence = ""
    for GenomeID in GenomeFNA:
        GenomeFNASequence += GenomeFNA[GenomeID]

    for SpacerID in GenomeSpacersFNA:
        #RandomSpacerPosition = random.randrange(1, len(GenomeFNASequence) - len(GenomeSpacersFNA[SpacerID])) - 1
        # # for real sequences from genome
        # RandomSpacersFNA[SpacerID + "_Random"] = GenomeFNASequence[RandomSpacerPosition: RandomSpacerPosition + len(GenomeSpacersFNA[SpacerID])]

        ## for permutated Genome
        RandomSpacersFNA[SpacerID + "_Random"] = ''.join(random.sample(GenomeSpacersFNA[SpacerID], len(GenomeSpacersFNA[SpacerID])))
        # RandomSpacerSequence = GenomeFNASequence[RandomSpacerPosition: RandomSpacerPosition + len(GenomeSpacersFNA[SpacerID])]
        # RandomSpacersFNA[SpacerID + "_Random"] = ''.join(random.sample(RandomSpacerSequence, len(RandomSpacerSequence)))

    for SpacerID in RandomSpacersFNA:
        for GenomeID in GenomeFNA:
            Position = GenomeFNA[GenomeID].find(RandomSpacersFNA[SpacerID])
            GenomeIDClear = GenomeID.split(" ")[0]
            if Position != -1:
                if not GenomeIDClear in Seeds:
                    Seeds[GenomeIDClear] = []
                Seeds[GenomeIDClear].append([Position, Position + len(RandomSpacersFNA[SpacerID])])
                break

    return Seeds, RandomSpacersFNA

def MakeRandomKPletsFromGenomeFNAPermutated(GenomeFNA, GenomeSpacersFNA, KPletSize):
    RandomKPlets = set()
    Seeds = dict()
    RandomSpacersFNA = dict()

    GenomeFNASequence = ""
    for GenomeID in GenomeFNA:
        GenomeFNASequence += GenomeFNA[GenomeID]

    #RandomSpacers = []
    for SpacerID in GenomeSpacersFNA:
        RandomSpacersFNA[SpacerID + "_Random"] = ''.join(random.sample(GenomeSpacersFNA[SpacerID], len(GenomeSpacersFNA[SpacerID])))

    for SpacerID in RandomSpacersFNA:
        RandomKPlets = RandomKPlets.union(GetKPlets(RandomSpacersFNA[SpacerID], KPletSize))

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

def GetKPletCoverage(SpacerKPlets, GenomeFNA, CRISPRSeeds, ViralFNA, SpacersFNA):
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

    for ContigID in ViralFNA:
        for KPlet in SpacersKPletsDict:
            for SubstringStart in [m.start() for m in re.finditer(KPlet, ViralFNA[ContigID])]:
                SpacersKPletsDict[KPlet][1] += 1
                SpacersHits = UpdateSpacerHits(KPlet, SpacersFNA, SubstringStart, "+", SpacersHits, "Viral", ContigID)

            for SubstringStart in [m.start() for m in re.finditer(Helper1603.ReverseComplement(KPlet), ViralFNA[ContigID])]:
                SpacersKPletsDict[KPlet][1] += 1
                SpacersHits = UpdateSpacerHits(KPlet, SpacersFNA, SubstringStart, "-", SpacersHits, "Viral", ContigID)
            # SpacersKPletsDict[KPlet][1] += len([m.start() for m in re.finditer(KPlet, ViralFNA[ContigID])])
            # SpacersKPletsDict[KPlet][1] += len([m.start() for m in re.finditer(Helper1603.ReverseComplement(KPlet), ViralFNA[ContigID])])

    # print output
    GenomeHits = 0
    ViralHits = 0
    for KPlet in SpacersKPletsDict:
        if SpacersKPletsDict[KPlet][0] > 0:
            GenomeHits += 1
        if SpacersKPletsDict[KPlet][1] > 0:
            ViralHits += 1

        #raise Exception("Kplet not found: " + KPlet + " " + "_".join(SpacersFNA[SpacerID].split("_")[0:3]))

    return GenomeHits, ViralHits, SpacersHits

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
            SpacerKPletHitsFile.write(TypeHits + "\t")

            SpacerKPletHitsFile.write("Viral\t")
            TypeHits = ""
            if "Viral" in SpacerHits[SpacerID][KPlet][1]:
                if len(SpacerHits[SpacerID][KPlet][1]["Viral"]) > 0:
                    for ContigID in SpacerHits[SpacerID][KPlet][1]["Viral"]:
                        if len(SpacerHits[SpacerID][KPlet][1]["Viral"][ContigID]) > 0:
                            TypeHits += ContigID.split(" ")[0] + ":" + ";".join([(str(GenomesKPletCoord[0]) + "," + str(GenomesKPletCoord[1])) for GenomesKPletCoord in SpacerHits[SpacerID][KPlet][1]["Viral"][ContigID]]) + "|"
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


LinkedInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/Linked.tsv"
#LinkedInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteLinked.tsv"
#LinkedInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelected.tsv"
SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"

SpacersFNA = Helper1603.ReadFasta(SpacersFNAFileName)

Counter = 0
for LinkedLine in open(LinkedInfoFileName):
    Counter += 1
    if Counter == TargetLineNo:
        break

ResFileName = "KPletsCoverageLinked_" + str(TargetLineNo) + ".tsv"
SpacerKPletHitsFileName = "KPletsCoverageHitsLinked_" + str(TargetLineNo) + ".tsv"
# ResFileName = "PermutatedKPletsCoverageLinked_" + str(TargetLineNo) + ".tsv"
# SpacerKPletHitsFileName = "PermutatedKPletsCoverageHitsLinked_" + str(TargetLineNo) + ".tsv"
# ResFileName = "PermutatedGenomeKPletsCoverageLinked_" + str(TargetLineNo) + ".tsv"
# SpacerKPletHitsFileName = "PermutatedGenomeKPletsCoverageHitsLinked_" + str(TargetLineNo) + ".tsv"


LinkedLineValues = LinkedLine[:-1].split("\t")
GenomeID = LinkedLineValues[0]
CRISPRs = LinkedLineValues[1].split(",")
VirusContigs = LinkedLineValues[2].split(",")

GenomeFNA, GenomeLength = GetGenomeFNA(GenomeID, CompleteGenomesDict)
ViralFNA = dict()
ViralFNA, ViralLength = GetViralFNA(GenomeID, VirusContigs)

CRISPRSeeds = dict()
for CRISPR in CRISPRs:
    CRISPRValues = CRISPR.split("_")
    if not CRISPRValues[0] in CRISPRSeeds:
        CRISPRSeeds[CRISPRValues[0]] = []
    CRISPRSeeds[CRISPRValues[0]].append([int(CRISPRValues[1]), int(CRISPRValues[2])])

CRISPRIDs = []
for CRISPR in CRISPRs:
    CRISPRIDs.append("_".join(CRISPR.split("_")[0:3]))

Spacers = []
GenomeSpacersFNA = dict()
for SpacerID in SpacersFNA:
    if any(CRISPRID in SpacerID for CRISPRID in CRISPRIDs):
        if (len(SpacersFNA[SpacerID]) > 22) and (len(SpacersFNA[SpacerID]) < 50):
            GenomeSpacersFNA[SpacerID] = SpacersFNA[SpacerID]

RandomSeeds, RandomSpacersFNA = MakeRandomKPletsFromGenomeFNA(GenomeFNA, GenomeSpacersFNA)
# RandomKPlets, RandomSeeds, RandomSpacersFNA = MakeRandomKPletsFromGenomeFNAPermutated(GenomeFNA, GenomeSpacersFNA, KPletSize)


with open(ResFileName, "w") as ResFile:
    with open(SpacerKPletHitsFileName, "w") as SpacerKPletHitsFile:
        for KPletSize in range(8, 23):
        #for KPletSize in range(21, 22):
        #for KPletSize in range(8, 9):
            #Counter = 0
            SpacerKPlets = set()
            for SpacerID in GenomeSpacersFNA:
                SpacerKPlets = SpacerKPlets.union(GetKPlets(GenomeSpacersFNA[SpacerID], KPletSize))
            RandomKPlets = set()
            for SpacerID in RandomSpacersFNA:
                RandomKPlets = RandomKPlets.union(GetKPlets(RandomSpacersFNA[SpacerID], KPletSize))

            GenomeHits, ViralHits, SpacerHits = GetKPletCoverage(SpacerKPlets, GenomeFNA, CRISPRSeeds, ViralFNA, GenomeSpacersFNA)
            RandomGenomeHits, RandomViralHits, RandomSpacerHits = GetKPletCoverage(RandomKPlets, GenomeFNA, RandomSeeds, ViralFNA, RandomSpacersFNA)

            ResFile.write(str(KPletSize) + "\t" + LinkedLine[:-1] + "\t" +
                          str(GenomeLength) + "\t" + str(ViralLength) + "\t" +
                          str(len(SpacerKPlets)) + "\t" + str(GenomeHits) + "\t" + str(ViralHits) + "\t" +
                          str(len(RandomKPlets)) +"\t" + str(RandomGenomeHits) + "\t" + str(RandomViralHits) + "\n")
                          # str(GenomeLength) + "\t" + "0" + "\t" +                          #
                          # str(len(SpacerKPlets)) + "\t" + str(GenomeHits) + "\t" + "0" + "\t" +
                          # str(len(RandomKPlets)) +"\t" + str(RandomGenomeHits) + "\t" + "0" + "\n")

            # SpacersHits[SpacerID][KPlet][1][Type][ContigID]
            WriteSpacerHits(SpacerKPletHitsFile, SpacerHits)
            WriteSpacerHits(SpacerKPletHitsFile, RandomSpacerHits)

#break
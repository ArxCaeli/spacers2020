import sys
sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
import Helper1603

# import argparse
import re
import random

Alphabet = {"A", "T", "G", "C"}

def GetKPlets(Sequence, FrameLength):
    KPlets = set()

    for I in range(FrameLength, len(Sequence) + 1):
        if any(x not in Alphabet for x in Sequence[I - FrameLength : I]):
            continue
        KPlets.add(Sequence[I - FrameLength : I])
        #KPlets.add(Helper1603.ReverseComplement(Sequence[I - FrameLength : I]))

    return KPlets


def MakeRandomKPletsFromGenomeFNA(GenomeFNA, GenomeSpacersFNA):
    Seeds = dict()
    RandomSpacersFNA = dict()

    GenomeFNASequence = ""
    for GenomeID in GenomeFNA:
        GenomeFNASequence += GenomeFNA[GenomeID]

    #RandomSpacers = []
    for SpacerID in GenomeSpacersFNA:
        RandomSpacerPosition = random.randrange(1, len(GenomeFNASequence) - len(GenomeSpacersFNA[SpacerID])) - 1
        ## for real sequences from genome
        RandomSpacersFNA[SpacerID + "_Random"] = GenomeFNASequence[RandomSpacerPosition: RandomSpacerPosition + len(GenomeSpacersFNA[SpacerID])]


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
        #print(GenomeFNA[ContigID][0:100])

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

# ap = argparse.ArgumentParser()
# ap.add_argument("-n", help = "KPletSize", required = True)
# opts = ap.parse_args()
#
# KPletSize = int(opts.n)

#SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/Lambda_spacers_overlap8_sample1k.fna"
SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/lambda-spacers_uniq.fasta"
GenomeFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/KD263.fasta"
ViralFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/lambda genome.fasta"

GenomeID = "gb_U00096.3_:c4641652-1"
ViralID = "gi|9626243|ref|NC_001416.1| Enterobacteria phage lambda, complete genome"
CRISPRsID = ""

#SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/T5_spacers_overlap8_sample1k.fna"
# SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/T5-spacers_uniq.fasta"
# GenomeFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/KD263.fasta"
# ViralFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/aclame_mge_1504_dna_seq(T5).fasta"
#
# GenomeID = "gb_U00096.3_:c4641652-1"
# ViralID = "mge:1504.1 Bacteriophage T5 virion, complete genome."
# CRISPRsID = ""


SpacersFNA = Helper1603.ReadFasta(SpacersFNAFileName)
GenomeFNA = Helper1603.ReadFasta(GenomeFNAFileName)
GenomeLength = len(GenomeFNA[GenomeID])
ViralFNA = Helper1603.ReadFasta(ViralFNAFileName)
ViralLength = len(ViralFNA[ViralID])


CRISPRSeeds = dict() # no CRISPR seeds

Spacers = []
GenomeSpacersFNA = dict()
for SpacerID in SpacersFNA:
    if (len(SpacersFNA[SpacerID]) > 29) and (len(SpacersFNA[SpacerID]) < 50):
        GenomeSpacersFNA[SpacerID] = SpacersFNA[SpacerID]

RandomSeeds, RandomSpacersFNA = MakeRandomKPletsFromGenomeFNA(GenomeFNA, GenomeSpacersFNA)

# #Helper1603.WriteFastaList(RandomSpacersFNA, "T5_RandomSpacers.faa")
# Helper1603.WriteFastaList(RandomSpacersFNA, "La_RandomSpacers.faa")
Prefix = "La"
for KPletSize in range(8, 31):
    # ResFileName = "KPletsCoverageLinked_T5_Overlap8_sample1k_" + str(KPletSize) + ".tsv"
    # SpacerKPletHitsFileName = "KPletsCoverageHitsLinked_T5_Overlap8_sample1k_" + str(KPletSize) + ".tsv"

    ResFileName = "KPletsCoverageLinked_" + Prefix + "_uniq_" + str(KPletSize) + ".tsv"
    SpacerKPletHitsFileName = "KPletsCoverageHitsLinked_" + Prefix + "_uniq_" + str(KPletSize) + ".tsv"

    # ResFileName = "KPletsCoverageLinked_lambda_overlap8_sample1k_" + str(KPletSize) + ".tsv"
    # SpacerKPletHitsFileName = "KPletsCoverageHitsLinked_lambda_overlap8_sample1k_" + str(KPletSize) + ".tsv"

    # ResFileName = "KPletsCoverageLinked_lambda_uniq_" + str(KPletSize) + ".tsv"
    # SpacerKPletHitsFileName = "KPletsCoverageHitsLinked_lambda_uniq_" + str(KPletSize) + ".tsv"

    with open(ResFileName, "w") as ResFile:
        with open(SpacerKPletHitsFileName, "w") as SpacerKPletHitsFile:
            #for KPletSize in range(8, 23):
            #for KPletSize in range(21, 22):
            #for KPletSize in range(8, 9):
            #Counter = 0
            print(Prefix + " Generating KPlets " + str(KPletSize))

            SpacerKPlets = set()
            for SpacerID in GenomeSpacersFNA:
                SpacerKPlets = SpacerKPlets.union(GetKPlets(GenomeSpacersFNA[SpacerID], KPletSize))
            RandomKPlets = set()
            for SpacerID in RandomSpacersFNA:
                RandomKPlets = RandomKPlets.union(GetKPlets(RandomSpacersFNA[SpacerID], KPletSize))

            print(Prefix + " KPletCount " + str(len(SpacerKPlets)))

            GenomeHits, ViralHits, SpacerHits = GetKPletCoverage(SpacerKPlets, GenomeFNA, CRISPRSeeds, ViralFNA, GenomeSpacersFNA)

            print(Prefix + " RandomHits")
            print(Prefix + " KPletCount " + str(len(RandomKPlets)))

            #RandomKPlets, RandomSeeds, RandomSpacersFNA = MakeRandomKPletsFromGenomeFNAPermutated(GenomeFNA, GenomeSpacersFNA, KPletSize)
            RandomGenomeHits, RandomViralHits, RandomSpacerHits = GetKPletCoverage(RandomKPlets, GenomeFNA, RandomSeeds, ViralFNA, RandomSpacersFNA)

            ResFile.write(str(KPletSize) + "\t" + GenomeID + "\t" + CRISPRsID + "\t" + ViralID + "\t" +
                          str(GenomeLength) + "\t" + str(ViralLength) + "\t" +
                          str(len(SpacerKPlets)) + "\t" + str(GenomeHits) + "\t" + str(ViralHits) + "\t" +
                          str(len(RandomKPlets)) +"\t" + str(RandomGenomeHits) + "\t" + str(RandomViralHits) + "\n")

            # SpacersHits[SpacerID][KPlet][1][Type][ContigID]
            WriteSpacerHits(SpacerKPletHitsFile, SpacerHits)
            WriteSpacerHits(SpacerKPletHitsFile, RandomSpacerHits)


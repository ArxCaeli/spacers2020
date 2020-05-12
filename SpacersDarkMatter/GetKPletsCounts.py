import sys
sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
import Helper1603

import subprocess
import re

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
    TmpFileName = "Tmp.tmp" + str(KPletSize)
    TmpFastaFileName = "Fasta.tmp" + str(KPletSize)
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

    return GenomeFNA, TotalLength

def GetViralFNA(VirusContigs):
    TmpFileName = "Tmp.tmp" + str(KPletSize)
    TmpFastaFileName = "Fasta.tmp" + str(KPletSize)

    with open(TmpFileName, "w") as TmpFile:
        for Contig in VirusContigs:
            TmpFile.write(Contig + "\n")

    subprocess.call("blastdbcmd -db nt -entry_batch " + TmpFileName + " > " + TmpFastaFileName, shell=True)

    ViralFNA = Helper1603.ReadFasta(TmpFastaFileName)

    TotalLength = 0
    for ID in ViralFNA:
        TotalLength += len(ViralFNA[ID])

    return ViralFNA, TotalLength



CompleteGenomesSetFileName = "/panfs/pan1/prokdata/db/Prok1603/Prok1603.pp.txt"
CompleteGenomesDict = dict()
for Line in open(CompleteGenomesSetFileName, "r"):
    if Line[0] == "#":
        continue
    LineValues = Line[:-1].split("\t")
    if not LineValues[0] in CompleteGenomesDict:
        CompleteGenomesDict[LineValues[0]] = []
    CompleteGenomesDict[LineValues[0]].append(LineValues[1])


LinkedInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/Linked.tsv"
SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"

SpacersFNA = Helper1603.ReadFasta(SpacersFNAFileName)

KPletSize = 6

ResFileName = "KPlets" + str(KPletSize) + ".tsv"
Counter = 0
with open(ResFileName, "w") as ResFile:
    for LinkedLine in open(LinkedInfoFileName):
        Counter += 1
        print(str(Counter) + "\t" + LinkedLine[:-1])
        LinkedLineValues = LinkedLine[:-1].split("\t")
        GenomeID = LinkedLineValues[0]
        CRISPRs = LinkedLineValues[1].split(",")
        VirusContigs = LinkedLineValues[2].split(",")

        CRISPRSeeds = dict()
        for CRISPR in CRISPRs:
            CRISPRValues = CRISPR.split("_")
            if not CRISPRValues[0] in CRISPRSeeds:
                CRISPRSeeds[CRISPRValues[0]] = []
            CRISPRSeeds[CRISPRValues[0]].append([int(CRISPRValues[1]), int(CRISPRValues[2])])

        # get kplets in spacers
        SpacerKPlets = set()
        CRISPRIDs = []
        for CRISPR in CRISPRs:
            CRISPRIDs.append("_".join(CRISPR.split("_")[0:3]))

        for SpacerID in SpacersFNA:
            if any(CRISPRID in SpacerID for CRISPRID in CRISPRIDs):
                SpacerKPlets = SpacerKPlets.union(GetKPliets(SpacersFNA[SpacerID], KPletSize))

        SpacersKPletsDict = dict()
        for KPlet in SpacerKPlets:
            SpacersKPletsDict[KPlet] = [0, 0]

        GenomeFNA, GenomeLength = GetGenomeFNA(GenomeID, CompleteGenomesDict)

        for ContigID in GenomeFNA:
            for KPlet in SpacersKPletsDict:
                for SubstringStart in [m.start() for m in re.finditer(KPlet, GenomeFNA[ContigID])]:
                    if not Helper1603.IsInSeeds(ContigID.split(" ")[0], SubstringStart - 100, SubstringStart + 100, CRISPRSeeds):
                        SpacersKPletsDict[KPlet][0] += 1

                for SubstringStart in [m.start() for m in re.finditer(Helper1603.ReverseComplement(KPlet), GenomeFNA[ContigID])]:
                    if not Helper1603.IsInSeeds(ContigID.split(" ")[0], SubstringStart - 100, SubstringStart + 100, CRISPRSeeds):
                        SpacersKPletsDict[KPlet][0] += 1

        # find kplets in viruses
        ViralFNA, ViralLength = GetViralFNA(VirusContigs)

        for ContigID in ViralFNA:
            for KPlet in SpacersKPletsDict:
                SpacersKPletsDict[KPlet][1] += len([m.start() for m in re.finditer(KPlet, ViralFNA[ContigID])])
                SpacersKPletsDict[KPlet][1] += len([m.start() for m in re.finditer(Helper1603.ReverseComplement(KPlet), ViralFNA[ContigID])])

        # print output
        GenomeHits = 0
        ViralHits = 0
        for KPlet in SpacersKPletsDict:
            GenomeHits += SpacersKPletsDict[KPlet][0]
            ViralHits += SpacersKPletsDict[KPlet][1]

        ResFile.write(LinkedLine[:-1] + "\t" + str(GenomeLength) + "\t" + str(ViralLength) + "\t" + str(len(SpacerKPlets)) +
                      "\t" + str(GenomeHits) + "\t" + str(ViralHits) + "\n")
        # for KPlet in SpacersKPletsDict:
        #     ResFile.write(KPlet + "\t" + str(SpacersKPletsDict[KPlet][0]) + "\t" + str(SpacersKPletsDict[KPlet][1]) + "\n")

        #break




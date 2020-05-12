import os
import subprocess

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


def GetCoverage(ContigRegionFileName, SpacerFNAFileName, SpacerSequence, ContigPatch, FixedLength):
    SpacerLength = len(SpacerSequence)

    HitsFileName = SpacerID + "_hits.tsv"
    subprocess.call(
        "blastn -outfmt \"7 qseqid sseqid slen sstart send evalue qseq sseq qstart qend bitscore score\" -word_size 5 -dust no -task=\"blastn-short\" -subject " +
        ContigRegionFileName + " -query " + SpacerFNAFileName + " > " + HitsFileName, shell=True)

    SpacerCoverage = [0] * SpacerLength
    for Line in open(HitsFileName):
        if Line[0] == "#":
            continue

        LineValues = Line[:-1].split("\t")
        StartPosition = int(LineValues[3]) - int(LineValues[8])

        for I in range(0, len(SpacerSequence)):
            if (StartPosition + I > 0) and (StartPosition + I < len(ContigPatch)):
                if SpacerSequence[I] == ContigPatch[StartPosition + I]:
                    SpacerCoverage[I] = 1
        break # skip other hits - take first one as a seed

    AdjustedCoverage = [0] * FixedLength

    for I in range(0, FixedLength):
        AdjustedCoverage[I] = SpacerCoverage[round(I * SpacerLength / FixedLength)]

    os.remove(HitsFileName)

    return AdjustedCoverage


def GetKPletSpacerCoverage(KPletHits, LookAroundDistance, SpacerID, FixedLength):
    Coverages = []

    ContigsFolder = "ViralContigs/"
    SpacersFolder = "Spacers/"

    SpacerFNAFileName = SpacerID + ".fna"
    if not os.path.exists(SpacersFolder + SpacerFNAFileName):
        subprocess.call("grep " + SpacerID + " /panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna -A 1 > " +
                        SpacersFolder + SpacerFNAFileName, shell = True)
    for Line in open(SpacersFolder + SpacerFNAFileName):
        if Line[0] == ">":
            continue
        SpacerSequence = Line[:-1]

    for Contig in KPletHits:
        if not os.path.exists(ContigsFolder + Contig + ".fna"):
            subprocess.call("blastdbcmd -db nt -entry " + Contig + " > " + ContigsFolder + Contig + ".fna", shell=True)

        ContigRegionFileName = Contig + "_region.fna"
        ContigFNA = ""
        for Line in open(ContigsFolder + Contig + ".fna"):
            if Line[0] != ">":
                ContigFNA += Line[:-1]


        Coordinates = sorted(KPletHits[Contig][0])
        LastCheckedCoord = 0
        for I in range(0, len(Coordinates)):
            Coord = Coordinates[I]
            if Coord <= LastCheckedCoord:
                continue

            ContigPatch = ContigFNA[Coord - LookAroundDistance - 1 : Coord + 2*LookAroundDistance - 1]

            with open(ContigRegionFileName, "w") as ContigRegionFile:
                ContigRegionFile.write(">" + Contig + "\n")
                ContigRegionFile.write(ContigPatch + "\n")


            Coverages.append(GetCoverage(ContigRegionFileName, SpacersFolder + SpacerFNAFileName, SpacerSequence, ContigPatch, FixedLength))
            LastCheckedCoord = Coord + 2*LookAroundDistance

            os.remove(ContigRegionFileName)

        #os.remove(Contig + ".fna")

    #os.remove(SpacerFNAFileName)

    return Coverages

FolderName = "../LinkedResults/"

LookAroundDistance = 32
FixedLength = 32
#Threshold = 22


#for Threshold in range(8, FixedLength + 1):
for Threshold in range(22, FixedLength + 1):
    KPletToCoverage = dict()

    ResFileName = "Coverage_" + str(Threshold) + ".tsv"
    with open(ResFileName, "w") as ResFile:
        for FileName in os.listdir(FolderName):
            if not "KPletsCoverageHitsLinked_" in FileName:
                continue

            Hits = dict()
            #ResFileName = FileName.split(".")[0] + "_densities.tsv"
            #with open(ResFileName, "w") as ResFile:
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
            #for KPletSize in [14,18,22]:
                #TotalCoverage = [0] * FixedLength
                if not KPletSize in KPletToCoverage:
                    KPletToCoverage[KPletSize] = [0] * FixedLength

                for SpacerID in Hits[KPletSize]:
                    HitTypePrefix = ""
                    if "_Random" in SpacerID:
                        continue
                    # if "_Random" in SpacerID:
                    #     HitTypePrefix = "Random"

                    Coverages = GetKPletSpacerCoverage(Hits[KPletSize][SpacerID][1], LookAroundDistance, SpacerID, FixedLength)
                    for Cov in Coverages:
                        if sum(Cov) >= Threshold:
                            for I in range(0, len(Cov)):
                                KPletToCoverage[KPletSize][I] += Cov[I]

                                # if Cov[I] == 1:
                                #     print(str(Threshold) + "\t" + HitTypePrefix + "ViralHit" + "\t" + str(KPletSize) + "\t" + str(I + 1))



        for KPletSize in KPletToCoverage:
            for I in range(0, len(KPletToCoverage[KPletSize])):
                ResFile.write(str(Threshold) + "\t" + HitTypePrefix + "ViralHit" + "\t" + str(KPletSize) + "\t" + str(I + 1) + "\t" + str(KPletToCoverage[KPletSize][I]) + "\n")

        # for Hit in HitDensities:
                #     ResFile.write(HitTypePrefix + "Genome" + "\t" + str(KPletSize) + "\t" + str(LookAroundDistance) + "\t" + str(Hit) + "\n")
                # HitDensities = GetLocalisationMetric(Hits[KPletSize][SpacerID][1], LookAroundDistance, KPletSize)
                # for Hit in HitDensities:
                #     ResFile.write(HitTypePrefix + "Viral" + "\t" + str(KPletSize) + "\t" + str(LookAroundDistance) + "\t" + str(Hit) + "\n")
                #break
            #break

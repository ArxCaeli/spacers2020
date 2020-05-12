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

def ReverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


# def GetOverlap(Seq1, Seq2):
#     for I in reversed(range(0, len(Seq1))): # equal lengths
#
# def Overlapped(Seq1, Seq2):


def GetLocalisationMetric(KPletHits, LookAheadDistance, HitsArrayID):
    SpacerDistalHits = 0
    SpacerDistalHits += GetLocalisationMetricForStrand(KPletHits, LookAheadDistance, HitsArrayID, 0)
    SpacerDistalHits += GetLocalisationMetricForStrand(KPletHits, LookAheadDistance, HitsArrayID, 1)

    return SpacerDistalHits

def GetLocalisationMetricForStrand(KPletHits, LookAheadDistance, HitsArrayID, Strand):
    DistalHits = 0
    for KPlet in KPletHits:
        KPletHitsCount = 0
        for Contig in KPletHits[KPlet][HitsArrayID]:
            Coordinates = sorted(KPletHits[KPlet][HitsArrayID][Contig][Strand])

            for OtherKPlet in KPletHits:
                if Contig in KPletHits[OtherKPlet][HitsArrayID]:
                    if (KPlet == OtherKPlet) or (KPlet == ReverseComplement(OtherKPlet)):
                        continue
                    OtherCoordinates = sorted(KPletHits[OtherKPlet][HitsArrayID][Contig][Strand])

                    for Coordinate in Coordinates:
                        for OtherCoordinate in OtherCoordinates:
                            if (Coordinate > OtherCoordinate + LookAheadDistance) or (Coordinate < OtherCoordinate - LookAheadDistance):
                                #DistalHits += 1
                                KPletHitsCount += 1
        if KPletHitsCount > 0:
            DistalHits += 1

    return DistalHits


FolderName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/LinkedResults/"

LookAheadDistance = 200

for FileName in os.listdir(FolderName):
    if not "KPletsCoverageHitsLinked_" in FileName:
        continue

    Hits = dict()
    ResFileName = FileName.split(".")[0] + "_distal_hits.tsv"
    print(ResFileName)
    with open(ResFileName, "w") as ResFile:
        for Line in open(FolderName + FileName):
            #CP001891.1_998980_1000411_21_spacer_1000229_32  GCCAGTCCGCATCTCGCCAAAA  1       Genome          Viral   NC_004775.2:34952,-
            LineValues = Line[:-1].split("\t")
            ArrayID = "_".join(LineValues[0].split("_")[0:3])

            KPletSize = len(LineValues[1])
            KPlet = LineValues[1]

            if not KPletSize in Hits:
                Hits[KPletSize] = dict() # for spacers

            if not LineValues[0] in Hits[KPletSize]:
                Hits[KPletSize][LineValues[0]] = dict()#[dict(), dict()] # Spacers with hits into genome / viral

            if not KPlet in Hits[KPletSize][LineValues[0]]:
                Hits[KPletSize][LineValues[0]][KPlet] = [dict(), dict()]


            if LineValues[4] != "":
                Hits[KPletSize][LineValues[0]][KPlet][0] = AddKPletHits(Hits[KPletSize][LineValues[0]][KPlet][0], LineValues[4])

            if LineValues[6] != "":
                Hits[KPletSize][LineValues[0]][KPlet][1] = AddKPletHits(Hits[KPletSize][LineValues[0]][KPlet][1], LineValues[6])

            #break
        for KPletSize in Hits:
            # if KPletSize != 22:
            #     continue

            for SpacerID in Hits[KPletSize]:
                HitTypePrefix = ""
                if "_Random" in SpacerID:
                    HitTypePrefix = "Random"

                SpacerDistalHits = GetLocalisationMetric(Hits[KPletSize][SpacerID], LookAheadDistance, 0)
                ResFile.write(SpacerID + "\t" + HitTypePrefix + "Genome" + "\t" + str(KPletSize) + "\t" + str(SpacerDistalHits) + "\n")
                SpacerDistalHits = GetLocalisationMetric(Hits[KPletSize][SpacerID], LookAheadDistance, 1)
                ResFile.write(SpacerID + "\t" + HitTypePrefix + "Viral" + "\t" + str(KPletSize) + "\t" + str(SpacerDistalHits) + "\n")

    #break

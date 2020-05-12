# to get CRISPRs-Genome-Virus tables

# CompleteGenomesContigsFileName = "/panfs/pan1/prokdata/db_tmp/Prok1603/Prok1603.pp.txt"
# CompleteGenomesContigsDict = dict()
# GenomeToContigsDict = dict()
#
# for Line in open(CompleteGenomesContigsFileName):
#     if Line[0] == "#":
#         continue
#
#     LineValues = Line[:-1].split("\t")
#
#     if LineValues[0] not in GenomeToContigsDict:
#         GenomeToContigsDict[LineValues[0]] = []
#     GenomeToContigsDict[LineValues[0]].append(LineValues[1])
#
#     CompleteGenomesContigsDict[LineValues[1]] = [LineValues[0], LineValues[5]]


GenomeToCRISPRs = dict()
CRISPRInfoFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_info_known_021718.tsv"

for Line in open(CRISPRInfoFileName):
    LineValues = Line[:-1].split("\t")

    if not LineValues[10] in GenomeToCRISPRs:
        GenomeToCRISPRs[LineValues[10]] = []

    GenomeToCRISPRs[LineValues[10]].append([LineValues[1], LineValues[2], LineValues[3], LineValues[9]])


## genome to virus
# SpaceHitsFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/SpacersSummary_clust_new.tsv"
# GenomeToCRISPRToVirus = dict()
#
# for Line in open(SpaceHitsFileName):
#     LineValues = Line[:-1].split("\t")
#
#     if LineValues[1] != "Phage":
#         continue
#
#     CRISPRID = LineValues[0].split("_")[0:3]
#     ContigID = CRISPRID[0]
#
#     if not ContigID in CompleteGenomesContigsDict:
#         continue
#
#     OrgID = CompleteGenomesContigsDict[ContigID][0]
#     if not OrgID in GenomeToCRISPRToVirus:
#         GenomeToCRISPRToVirus[OrgID] = [GenomeToCRISPRs[OrgID],[]]
#
#     GenomeToCRISPRToVirus[OrgID][1].append(LineValues[2])
#
# with open("Linked.tsv", "w") as ResFile:
#     for Genome in GenomeToCRISPRToVirus:
#         ResFile.write(Genome + "\t" + ",".join(["_".join(x) for x in GenomeToCRISPRToVirus[Genome][0]]) + "\t" + ",".join(GenomeToCRISPRToVirus[Genome][1]) + "\n")SpaceHitsFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/SpacersSummary_clust_new.tsv"

# get complete genomes to CRISPR
CompleteGenomesFileName = "/panfs/pan1/prokdata/db_tmp/Prok1603/Prok1603.gg.txt"
CompleteGenomes = set()

for Line in open(CompleteGenomesFileName):
    if Line[0] == "#":
        continue

    LineValues = Line[:-1].split("\t")

    CompleteGenomes.add(LineValues[0])


CompleteGenomeToCRISPR = dict()

for Genome in GenomeToCRISPRs:
    if not Genome in CompleteGenomes:
        continue

    for CRISPR in GenomeToCRISPRs[Genome]:
        if not Genome in CompleteGenomeToCRISPR:
            CompleteGenomeToCRISPR[Genome] = []

        CompleteGenomeToCRISPR[Genome].append(CRISPR)

with open("CompleteLinked.tsv", "w") as ResFile:
    for Genome in CompleteGenomeToCRISPR:
        ResFile.write(Genome + "\t" + ",".join(["_".join(x) for x in CompleteGenomeToCRISPR[Genome]]) + "\n")
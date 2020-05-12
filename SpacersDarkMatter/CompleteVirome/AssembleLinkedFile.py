import os


AccessionFolder = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/"

OrganismToVirus = dict()
for FileName in os.listdir(AccessionFolder):
    #if ".acc" in FileName:
    if "Mycobacterium.acc" in FileName:
        OrgName = FileName.split(".")[0]
        Accessions = []
        for Line in open(FileName):
            Accessions.append(Line[:-1])
        OrganismToVirus[OrgName] = Accessions

for Line in open("CompleteLinked.tsv"):
    LineValues = Line[:-1].split("\t")

    for OrgName in OrganismToVirus:
        if OrgName in LineValues[0]:
            print(Line[:-1] + "\t" + ",".join(OrganismToVirus[OrgName]))
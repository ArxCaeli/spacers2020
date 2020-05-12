import os

Path = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/"
for FileName in os.listdir(Path):
    Prefix = FileName.split(".")[0]
    if not "spacers" in Prefix:
        continue
    print(Prefix)

    ResFileName = Path + Prefix + "_named.fasta"
    with open(ResFileName, "w") as ResFile:
        Counter = 0
        for Line in open(Path + FileName):
            if Line[0] == ">":
                Counter += 1
                ResFile.write(">" + Prefix + "_" + str(Counter) + "\n")
            else:
                ResFile.write(Line)
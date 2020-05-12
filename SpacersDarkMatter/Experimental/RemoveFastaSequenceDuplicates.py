#InputFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/lambda-spacers_named.fasta"
InputFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/T5new-spacers_named.fasta"

Unique = False
LineCount = 0
SpacersSet = set()
with open(InputFileName, "r") as InputFile:
    for Line in InputFile:
        if Line[0] == ">":
            ID = Line[:-1]
        else:
            if not Line[:-1] in SpacersSet:
                print(ID)
                print(Line[:-1])
                SpacersSet.add(Line[:-1])
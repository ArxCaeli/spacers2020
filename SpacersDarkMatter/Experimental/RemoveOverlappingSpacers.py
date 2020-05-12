import sys
sys.path.insert(0, "/home/shmakovs/Fishing/Pipeline/Git/Fishing/")
import Helper1603

def IsOverlapped(Spacer1, Spacer2):
    if Spacer1 in Spacer2:
        return True
    if Spacer2 in Spacer1:
        return True

    OverlapThreshold = 8

    for I in range(0, len(Spacer1) - OverlapThreshold + 1):
        Substring = Spacer1[I:len(Spacer1)]
        if Substring == Spacer2[0:len(Substring)]:
            return True

    for I in range(0, len(Spacer1) - OverlapThreshold + 1):
        Substring = Spacer1[0:len(Spacer1) - I]
        if Substring == Spacer2[-len(Substring):]:
            return True


    return False

# IsOverlapped("abc", "abcd")
# IsOverlapped("1234567890abc", "0abcdzxcvbnm,")
# IsOverlapped("1234567890abc", "qwerty12")


#SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/T5-spacers_uniq.fasta"
SpacersFNAFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/ExperimentalResults/Data/lambda-spacers_uniq.fasta"
SpacersFasta = Helper1603.ReadFasta(SpacersFNAFileName)

GoodSpacerSet = set()
for SpacerID in SpacersFasta:
    if not any(IsOverlapped(SpacersFasta[SpacerID], x) for x in GoodSpacerSet):
        print(">" + SpacerID)
        print(SpacersFasta[SpacerID])
        GoodSpacerSet.add(SpacersFasta[SpacerID])

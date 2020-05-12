def FillKPletInSpacer(SpacerCoverage, StartPositions, KPletSize, SpacerLength, FixedLength, SpacerDirection):
    # for Position in StartPositions.split(";"):
    #     Start = int(Position)
    #     AdjustedStart = round((Start / SpacerLength) * (FixedLength - 1))
    #     AdjustedEnd = round(((Start + KPletSize) / SpacerLength) * (FixedLength - 1))
    #     #print(Position + "\t" + str(Start) + "\t" + str(KPletSize) + "\t" + str(AdjustedStart) + "\t" + str(AdjustedEnd))
    #
    #     if SpacerDirection == "-":
    #         tmp = AdjustedStart
    #         AdjustedStart = FixedLength - AdjustedEnd
    #         AdjustedEnd = FixedLength - tmp
    #     for I in range(AdjustedStart, AdjustedEnd):
    #         SpacerCoverage[I] = 1
    # return SpacerCoverage

    #circular variant
    for Position in StartPositions.split(";"):
        Start = int(Position)
        End = int(Start + KPletSize)
        InitialArray = [0] * SpacerLength

        for I in range(Start, End):
            InitialArray[I] = 1
            if I < KPletSize:
                InitialArray[SpacerLength - I - 1] = 1
            if I >= SpacerLength - KPletSize:
                InitialArray[SpacerLength - Start - KPletSize + I - (SpacerLength - KPletSize)] = 1
                InitialArray[I - Start] = 1

        print(InitialArray)
        ArrayLen = len(InitialArray)
        for I in range(0, FixedLength):
            AdjustedPosition = min(round(I * ArrayLen / FixedLength), len(InitialArray) - 1)
            if InitialArray[AdjustedPosition] == 1:
                SpacerCoverage[I] = 1
    return SpacerCoverage

FixedLength = 32
SpacerCoverage = [0] * FixedLength
SpacerLength = 26
StartPositions = "8"
KPletSize = 16

# FixedLength = 16
# SpacerCoverage = [0] * FixedLength
# SpacerLength = 8
# StartPositions = "0"
# KPletSize = 2

FillKPletInSpacer(SpacerCoverage, StartPositions, KPletSize, SpacerLength, FixedLength, "+")

print(SpacerCoverage)

#print(round(15 * 0.5))


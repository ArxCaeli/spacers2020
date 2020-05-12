import subprocess
import time

import os

qsub_parameters = "-v SGE_FACILITIES -v SGE_NOMAIL -v SGE_SUMMARY=\"stdout\""
qsub_parameters += " -P unified -l h_rt=72000,h_vmem=10G,mem_free=10G,m_mem_free=10G,ul2=5 -b y"

#GenomesToCRISPRsFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteVirome/LinkedSelected.tsv"
#GenomesToCRISPRsFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/CompleteLinked.tsv"
GenomesToCRISPRsFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/DarkMatter2018/Linked.tsv"

Counter = 0
for Line in open(GenomesToCRISPRsFileName):
    Counter += 1

    # if Counter < 326:
    #     continue

    FileName = "KPletsCoverageLinked_" + str(Counter) + ".tsv"
    # if os.path.exists(FileName) and os.stat(FileName).st_size != 0:
    #     continue

    #     if os.stat(FileName).st_size != 0:
    #         for Line in open(FileName):
    #             LineValues = Line[:-1].split("\t")
    #         if LineValues[0] != str(22):
    #             print(FileName)
    #
    #         continue

    #print(FileName)

    command = "sh RunCoverageCountJob.sh " + str(Counter)
    #print(command)
    subprocess.call(
        "qsub " + qsub_parameters + " -N Spacers_" + str(Counter) + " -e std_err.txt -o std_out.txt " + command, shell=True)  # , stdout = devnull, stderr = devnull)
    subprocess.call("$qsub", shell=True)  # , stdout = devnull, stderr = devnull)

    time.sleep(0.1)

    #break
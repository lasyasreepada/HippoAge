#!/usr/bin/env python
'''
Run this script on a personal computer to download dicoms from ida.loni.usc.edu. 
The script will not work on the bscsub cluster.
'''

import os

#####################################
''' UPDATE THESE VARIABLES BEFORE RUNNING '''
## date for download file name
date="aibl_dl_Mar2025"

## Advanced download filepaths from ida.loni.usc.edu
downloads = {
    "ALL_MPRAGE":["https://ida.loni.usc.edu/download/files/ida1/828f1358-3b31-4c0f-a288-d2edf04c58b1/AIBL.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/42ef5b33-6e27-440c-8d6a-a94db5141f36/AIBL1.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/43f05e86-3583-4a45-8414-a770606846ec/AIBL2.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/346ed3e3-6b2e-476f-9da0-658d13a4a13a/AIBL3.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/eab6905e-6074-4d35-9773-835999a60de5/AIBL4.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/55d86462-0dd9-4771-923a-9b1a1427035a/AIBL5.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/7af2a4c1-cd7a-46fe-89c8-c40034a97cab/AIBL6.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/5fb4c647-ea56-4d1d-ab81-e875a34869b8/AIBL7.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/3d5d6f69-156e-4ffd-b4d9-f5ffcdc59272/AIBL8.zip",
                  "https://ida.loni.usc.edu/download/files/ida1/cf260cbd-91d8-44a6-93b4-52e6abad460d/AIBL9.zip"],
}
######################################

## directory on personal computer
basedir="/Users/lasyasreepada/Projects/hippoage/data/AIBL/MRI"
targetdir=f"{basedir}/{date}" 
filestounzip=set(()) #only once for each modality

os.system(f"mkdir {targetdir}")

for key, value in downloads.items():
    for i in range(0,len(value)):
        if len(value) > 1:
            output=f"{targetdir}/{i}_{key}.zip"
            filestounzip.add(key) 
        else:
            output = f"{targetdir}/{key}.zip"
        print(f"ready to download {key} {i}")
        os.system(f"curl {value[i]} -o {output}")

print()
# print("Command to copy all files & containing folder to cluster:")
# print(f"rsync -avh --progress --ignore-existing {targetdir} sreepada@bscsub.pmacs.upenn.edu:/project/wolk/aibl")
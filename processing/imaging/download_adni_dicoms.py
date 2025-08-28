#!/usr/bin/env python
'''
Run this script on a personal computer to download dicoms from ida.loni.usc.edu.
The script will not work on the bscsub cluster.
'''

import os

#####################################
''' UPDATE THESE VARIABLES BEFORE RUNNING '''
## directory on personal computer
basedir="/Users/lasyasreepada/Projects/hippoage/data/ADNI/MRI"
## cluster username for scp to bscsub cluster
bsc_username="sreepada"
## date for download file name
date="adni_dl_Feb2025"
## Advanced download filepaths from ida.loni.usc.edu
downloads = {
    "ALL_MPRAGE":['https://ida.loni.usc.edu/download/files/ida1/a9928837-754e-4368-8eaf-72567a9e77a5/ADNI.zip',
                  'https://ida.loni.usc.edu/download/files/ida1/132ef37e-255d-4f44-a66f-b2d558f4ec0b/ADNI1.zip',
                  'https://ida.loni.usc.edu/download/files/ida1/3a955538-3ca2-4d95-b670-346ea002c78c/ADNI2.zip',
                  'https://ida.loni.usc.edu/download/files/ida1/50f7f82e-8017-4721-b7e7-b27df9eb15d6/ADNI3.zip',
                  'https://ida.loni.usc.edu/download/files/ida1/9eb1b4fd-9aeb-4af5-8d11-4bf55563c698/ADNI4.zip'],
}
######################################

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

# print()
# print("Command to copy all files & containing folder to cluster:")
# print(f"rsync -avh --progress --ignore-existing {targetdir} {bsc_username}@bscsub.pmacs.upenn.edu:/project/wolk/lps")

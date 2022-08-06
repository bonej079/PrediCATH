#!/usr/bin/env bash

# More info at https://unix.stackexchange.com/questions/189878/parallelise-rsync-using-gnu-parallel

sudo mkdir -p /tmp/ramdisk
sudo mount -t tmpfs -o size=50G tmpfs /tmp/ramdisk

# rsync -ahrP ~/apps/cath-tools-genomescan /mnt/ramdisk
rsync -ahrP ./cath-tools-genomescan /mnt/ramdisk

# CATH4_1
# rsync -ahrP /mnt/DataDrive/Data/UCL/CATHv4_1/hmms/ /mnt/ramdisk/cath-tools-genomescan/cath4_1/

# CATH4_2
# rsync -ahrP /mnt/DataDrive/Data/UCL/CATHv4_2/hmms/ /mnt/ramdisk/cath-tools-genomescan/cath4_2/

# CATH4_3
rsync -ahrP ./cath-tools-genomescan/CATHv4_3/hmms/ /mnt/ramdisk/cath-tools-genomescan/cath4_3/

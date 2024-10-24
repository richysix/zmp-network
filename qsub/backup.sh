#!/usr/bin/env bash
# backup.sh - Script to backup results directory to Sharepoint
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=1G
#$ -o copy.o
#$ -e copy.e

module load rclone/1.65.2

dir=zmp-network/nf/results
echo "Starting $dir backup 1" 1>&2
rclone copy /data/scratch/bty114/$dir/ sharepoint-qmul-buschlab:Projects/$dir/ \
 --include "*.{doc,docx,xls,xlsx,xlsm,ppt,pptx,html,png}" \
 --ignore-size --ignore-checksum --drive-acknowledge-abuse --update --copy-links

echo "Starting $dir backup 2" 1>&2
rclone copy /data/scratch/bty114/$dir/ sharepoint-qmul-buschlab:Projects/$dir/ \
 --exclude "*.{doc,docx,xls,xlsx,xlsm,ppt,pptx,html,png}" \
 --drive-acknowledge-abuse --update --copy-links

module unload rclone/1.65.2

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2024. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007

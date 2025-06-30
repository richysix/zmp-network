#!/usr/bin/env bash
# aggregate_networks.sh - Script to test the aggregate_networks.R script
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=64G
#$ -o aggregate.o
#$ -e aggregate.e

module load R/4.4.0

# test script
# Rscript ~/checkouts/zmp-network/bin/aggregate_networks.R \
# --annotation /data/scratch/bty114/zmp-network/nf/reference/danio_rerio-e99-annotation.test.txt \
# --go_annotation /data/scratch/bty114/zmp-network/nf/reference/danio_rerio_e109_go.test.txt \
# --samples_file /data/scratch/bty114/zmp-network/nf/test/expt-sample-condition-tfap2-plus.tsv \
# --expts_file /data/scratch/bty114/zmp-network/nf/test/expts.txt \
# --num_orderings 10 --debug \
# /data/scratch/bty114/zmp-network/nf/test/zmp_ph*-tpm-filtered-orig.test.mat

Rscript ~/checkouts/zmp-network/bin/aggregate_networks.R \
--annotation /data/scratch/bty114/zmp-network/nf/reference/danio_rerio-e99-annotation.txt \
--go_annotation /data/scratch/bty114/zmp-network/nf/reference/danio_rerio_e109_go.txt \
--samples_file /data/scratch/bty114/zmp-network/nf/test/expt-sample-condition-tfap2-plus.tsv \
--expts_file /data/scratch/bty114/zmp-network/nf/test/expts.txt \
--num_orderings 10 --debug \
/data/scratch/bty114/zmp-network/nf/test/zmp_ph*-tpm-filtered-orig.mat

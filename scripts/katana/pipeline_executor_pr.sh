#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-12 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-24 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 25-30 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 31-36 pipeline_executor_pr.pbs
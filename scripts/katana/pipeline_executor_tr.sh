#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-12 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-18 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 19-22 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 23-26 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 27-30 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 31-34 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 35-38 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 39-42 pipeline_executor_tr.pbs

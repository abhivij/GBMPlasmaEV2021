#!/bin/bash

qsub -l select=1:ncpus=8:mem=32gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-9 RFE_execution.pbs
qsub -l select=1:ncpus=8:mem=64gb,walltime=100:00:00 -M a.vijayan@unsw.edu.au -m ae -J 10-18 RFE_execution.pbs
qsub -l select=1:ncpus=8:mem=64gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 19-27 RFE_execution.pbs
#!/bin/bash

qsub -l select=1:ncpus=8:mem=32gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-3 RFE_execution.pbs
qsub -l select=1:ncpus=8:mem=32gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 4-6 RFE_execution.pbs
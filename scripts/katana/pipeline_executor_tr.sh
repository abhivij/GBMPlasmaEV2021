#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-12 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-18 pipeline_executor_tr.pbs
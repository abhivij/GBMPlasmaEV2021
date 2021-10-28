#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 1-6 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 7-12 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 13-18 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 19-24 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 25-30 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 31-40 pipeline_executor.pbs
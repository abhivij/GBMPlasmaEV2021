#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -J 1-6 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -J 7-18 pipeline_executor.pbs
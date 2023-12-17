#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-12 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-24 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 25-30 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 31-36 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 37-48 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 49-60 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 61-64 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 65-68 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 69-72 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 73-76 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 77-80 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 81-84 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 85-88 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 89-92 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 93-96 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 97-100 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 101-104 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 105-108 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 109-112 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 113-116 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 117-120 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 121-124 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 125-128 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 129-132 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 133-136 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 137-140 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 141-144 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 145-147 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 148-150 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 151-153 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 154-156 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 157-161 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 162-164 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 165-168 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 169-172 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 173-176 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 177-181 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 182-186 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 187-190 pipeline_executor_pr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 191-194 pipeline_executor_pr.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 195-198 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 199-202 pipeline_executor_pr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 203-206 pipeline_executor_pr.pbs
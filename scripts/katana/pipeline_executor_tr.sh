#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-12 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-18 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 19-22 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 23-26 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 27-30 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 31-34 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 35-38 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 39-42 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 43-46 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 47-50 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 51-54 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 55-58 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 59-62 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 63-66 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 67-74 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 75-84 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 85-88 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 89-92 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 93-96 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 97-100 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 101-104 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 105-108 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 109-112 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 113-116 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=48:00:00 -M a.vijayan@unsw.edu.au -m ae -J 117-120 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 121-126 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 127-132 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 133-137 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 138-140 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 141-147 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=64gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 148-152 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 153-156 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 157-160 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 161-164 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 165-169 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 170-173 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 174-181 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 182-185 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 186-189 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 190-193 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 194-197 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 198-203 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 204-209 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 210-213 pipeline_executor_tr.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 214-217 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 218-221 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 222-225 pipeline_executor_tr.pbs


qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 226-227 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 228-229 pipeline_executor_tr.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 230-234 pipeline_executor_tr.pbs
#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 1-6 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 7-12 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 13-18 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 19-24 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 25-30 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 31-40 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 41-46 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 47-52 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 53-58 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 59-72 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 73-74 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 75-84 pipeline_executor.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 85-94 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 95-104 pipeline_executor.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 105-114 pipeline_executor.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=10:00:00 -J 115-124 pipeline_executor.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 125-136 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 137-138 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 139-158 pipeline_executor.pbs


qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -J 159-162 pipeline_executor.pbs


qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -J 163-189 pipeline_executor.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -J 190-218 pipeline_executor.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -J 219-222 pipeline_executor.pbs




#transcriptomics new comparisons
qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 223-232 pipeline_executor.pbs

#proteomics new comparisons
qsub -l select=1:ncpus=16:mem=124gb,walltime=24:00:00 -J 233-242 pipeline_executor.pbs



#new comparisons subset
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -J 243-247 pipeline_executor.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -J 248-254 pipeline_executor.pbs
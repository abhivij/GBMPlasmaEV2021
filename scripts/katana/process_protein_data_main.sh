#!/bin/bash

qsub -l select=1:ncpus=8:mem=124gb,walltime=24:00:00 -J 1-23 process_protein_data_main.pbs
#!/bin/bash

qsub -l select=1:ncpus=8:mem=128gb,walltime=24:00:00 -J 1-16 process_protein_data_main.pbs
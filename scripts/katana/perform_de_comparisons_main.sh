#!/bin/bash

qsub -l select=1:ncpus=8:mem=64gb,walltime=12:00:00 -J 2-7 perform_de_comparisons_main.pbs
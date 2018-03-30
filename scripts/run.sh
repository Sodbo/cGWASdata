#!/bin/bash

# Sodbo Sharapov

echo "Running cGAS analysis using biochemical distances"
echo "Running main_BN.R"
Rscript main_BN.R

echo "Running cGAS analysis using GGM distances"
echo "Running main_GGM.R"
Rscript main_GGM.R

echo "Creating Figure 1"
echo "Running Figure_1.R"
Rscript Figure_1.R

echo "Creating Figure 2"
echo "Running Figure_2.R"
Rscript Figure_2.R



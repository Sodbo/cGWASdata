#!/bin/bash

# Sodbo Sharapov

echo "Running cGAS analysis using biochemical distances"
echo "Running main_BN.R"
Rscript code/main_BN.R

echo "Running cGAS analysis using GGM distances"
echo "Running main_GGM.R"
Rscript code/main_GGM.R

echo "Creating Supplementary Table 1A"
echo "Running Suppl.table.1A.R"
Rscript code/Suppl.table.1A.R

echo "Creating Supplementary Table 1B"
echo "Running Suppl.table.1B.R"
Rscript code/Suppl.table.1B.R

echo "Creating Supplementary Table 2"
echo "Running Suppl.table.2.R"
Rscript code/Suppl.table.2.R

echo "Creating Figure 1"
echo "Running Figure_1.R"
Rscript code/Figure_1.R

echo "Creating Figure 2"
echo "Running Figure_2.R"
Rscript code/Figure_2.R

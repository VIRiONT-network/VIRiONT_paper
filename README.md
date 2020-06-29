# MinION_HBV

# Quick presentation

This pipeline is designed to analyse viral data from nanopore sequencing. It was primary developated to analyse HBV data, but was extended also to HDV. However, this pipeline can be adapted for other virus aswell.

# Workflow

The pipeline take as input unmerged fastq from demultiplexed nanopore data.
For each barcode, here is the global workflow:
Step1 => merging these fastq
Step2 => filtering merged fastq
Step3 => blast research on each reference given
Step4 => build consensus sequence on the majoritary reference

# Requirements & Tools

This pipeline use several tools.
We strongly recommend you using conda for a better tool management. an environment.yml file is provided for a quick and easy install of the necessary tools.


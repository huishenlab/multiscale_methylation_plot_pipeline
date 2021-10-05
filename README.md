# Pipeline to create the multiscale methylation plot

## Components of the workflow
  1. Calculate x-values for plot
  2. Add chromosome lengths to x-values file (needed for next step)
  3. Create windows around x-values for calculating average methylation value
  4. Calculate average methylation values for each window
  5. Calculate cross-sample mean and standard deviation for each window size

## Running the pipeline
The pipeline can be run with on the command line `snakemake --cores 1 --use-envmodules` or submitted to a job scheduler
on a cluster (a PBS script is provided: `qsub bin/run_snakemake_workflow.sh`)

## Creating the multiscale plot
An example script to generate the multiscale plot is provided in `scripts/create_multiscale_plot.R`.

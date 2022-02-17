# Pipeline to create the multiscale methylation plot

## Components of the workflow
  1. Calculate x-values for plot
  2. Add chromosome lengths to x-values file (needed for next step)
  3. Create windows around x-values for calculating average methylation value
  4. Calculate average methylation values for each window
  5. Calculate cross-sample mean and standard deviation for each window size

## Dependencies
The following dependencies are downloaded with `--use-conda`, otherwise you must have them in your PATH.
  - `snakemake` (version 6.0+)
  - `bedtools`
  - `python` (version 3.7+)
    - `pandas`

## Running the pipeline
  - Clone or download the repo.
    - `git clone git@github.com:huishenlab/multiscale_methylation_plot_pipeline.git`, OR
    - Download from the [releases page](https://github.com/huishenlab/multiscale_methylation_plot_pipeline/releases).
  - Place *gzipped* BED files in the `raw_data` directory.
    - Alternatively, you can specify the path to your *gzipped* BED files in `config/config.yaml`.
  - Replace the temporary sample names in `config/samples.tsv` with your sample names (everything before `.bed.gz` in
  your input files.
    - Alternatively, you can specify the path to your sample sheet in `config/config.yaml`. If you use your own file, be
    sure to include "sample" as the header line.
  - Modify the `config/config.yaml` file with your chosen inputs.
    - At minimum, you will need to specify the FAI index file location for your genome. All inputs have defaults set.
  - The pipeline can then be run on the command line (`snakemake --cores 1 --use-conda`) or submitted to a job scheduler
  on a cluster (a PBS script is provided: `qsub bin/run_snakemake_workflow.sh`).

## After the pipeline finishes running
  - Three directories will be created in your specified output directory.
    - `analysis/`: results from the pipeline
      - `means/`: mean methylation values in bins of varying size for each input sample
      - `stats/`: average and standard deviation for values in each bin across all samples
      - `x_values/`: bins used for calculating mean values
    - `benchmarks/`: benchmarking data compiled when each rule runs
    - `logs/`: log files generated as rules finish

## Creating the multiscale plot
The multiscale plot can be created using [bisplotti](https://github.com/huishenlab/bisplotti), specifically the
`multiscaleMethylationPlot()` function. Input is a specific sample in `analysis/means/` or the average for all samples
from `analysis/stats/`. If you have two (or more) sample groups that you processed together, you can also find per-bin
average values for a specified set of samples using `multiscaleGroupAverage()` in `bisplotti`.

## Example dataset
An example dataset has been provided as a `.zip` file on the releases page. Instructions for using can be found in a
README once the file has been unzipped.

import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version

# Set minimum version of Snakemake
min_version('6.0')

# Set configuration file location
configfile:
    'config/config.yaml'

# Constrain possible values of power
wildcard_constraints:
    power = '\d+\.\d*',

# Load sample base names from sample file
samples = pd.read_table(config['samples'], dtype=str).set_index(['sample'], drop=False)

# Generate list of requested power values
# Will create bins of size 10^X (for X in powers_list) that are centered around each value in create_x_values
powers_list = [
    round(i, 2) for i in np.arange(
        config['bounds']['min'],
        config['bounds']['max'] + config['bounds']['step'],
        config['bounds']['step']
    )
]
powers = pd.DataFrame({'pows': powers_list}).set_index(['pows'], drop=False)

# Create a tag to use in output files based on the x value step size.
def create_tag(step):
    if step < 1000: # Just bp
        return f'{step}bp'
    elif step < 1000000: # kilobases
        return f'{int(step/1000)}kb'
    else: # Anything bigger will be defined as megabases
        return f'{int(step/1000000)}Mb'

x_tag = create_tag(config['x_step'])

# Run all the things
rule all:
    input:
        f"analysis/x_values/x_values_no_lengths_{x_tag}.bed",
        f"analysis/x_values/x_values_{x_tag}.bed",
        expand(f"analysis/x_values/x_values_{x_tag}_{{powers.pows}}.bed", powers=powers.itertuples()),
        expand(f"analysis/means/{{samples.sample}}_{{powers.pows}}.bed.gz", samples=samples.itertuples(), powers=powers.itertuples()),
        expand(f"analysis/stats/cross_sample_mean_{{powers.pows}}.bed.gz", powers=powers.itertuples()),
        expand(f"analysis/stats/cross_sample_stddev_{{powers.pows}}.bed.gz", powers=powers.itertuples()),

# Rule to create x-axis values that will be used when plotting the multiscale plot
rule create_x_values:
    input:
        config['reference_fai'],
    output:
        bed = f"analysis/x_values/x_values_no_lengths_{x_tag}.bed",
    envmodules:
        config['envmodules']['bedtools'],
    params:
        outdir = directory("analysis/x_values/"),
        step = config['x_step'],
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['low_memory_gb'],
        walltime = config['hpc_parameters']['low_walltime'],
    shell:
        """
        set +o pipefail
        mkdir -p {params.outdir}
        bedtools makewindows \
            -w 1 -s {params.step} \
            -g {input} | \
        sort -k1,1 -k2,2n > {output.bed}
        """

# Rule to add the chromosome length as a fourth column in the x-axis values BED file, as well as restricting to
# canonical chromsomes
# Fourth column will be needed when creating the windows centered around the x-values
# Don't want the windows running past the end of the chromosomes!
rule add_chr_length:
    input:
        bed = f"analysis/x_values/x_values_no_lengths_{x_tag}.bed",
        ref = config['reference_fai'],
    output:
        bed = f"analysis/x_values/x_values_{x_tag}.bed",
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['low_memory_gb'],
        walltime = config['hpc_parameters']['low_walltime'],
    script:
        "scripts/add_chrom_length.py"

# Rule to create the windows that will be used to calculate the average methylation values that go into the multiscale
# plot
# Windows are centered around the x-values calculated previously and have a width of 10^(power)
# Beginning and end are adjusted as needed for the start and end of the chromosome
rule create_x_value_windows:
    input:
        bed = f"analysis/x_values/x_values_{x_tag}.bed",
    output:
        bed = f"analysis/x_values/x_values_{x_tag}_{{power}}.bed",
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['low_memory_gb'],
        walltime = config['hpc_parameters']['low_walltime'],
    script:
        "scripts/create_x_val_windows.py"

# Rule to calculate the average methylation value in the windows calculated previously
# If config['target'] is defined, then the input BED file will be subset to data that overlaps config['target']
rule calculate_means:
    input:
        windows = f"analysis/x_values/x_values_{x_tag}_{{power}}.bed",
        data    = f"raw_data/{{sample}}.bed.gz",
    output:
        bed_gz = f"analysis/means/{{sample}}_{{power}}.bed.gz",
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['mid_memory_gb'],
        walltime = config['hpc_parameters']['mid_walltime'],
    params:
        outdir = directory("analysis/means/"),
        target = config['targets']
    envmodules:
        config['envmodules']['bedtools'],
    shell:
        """
        set +o pipefail
        mkdir -p {params.outdir}

        if [[ -f {params.target} ]]; then
            bedtools intersect \
                -a {input.data} \
                -b {params.target} \
                -sorted | \
            bedtools map \
                -prec 2 \
                -a {input.windows} \
                -b - \
                -c 4 -o mean | \
            cut -f 1,4-6 | \
            sort -k1,1 -k2,2n | \
            gzip > {output.bed_gz}
        else
            bedtools map \
                -prec 2 \
                -a {input.windows} \
                -b {input.data} \
                -c 4 -o mean | \
            cut -f 1,4-6 | \
            sort -k1,1 -k2,2n | \
            gzip > {output.bed_gz}
        fi
        """

# Collect all files that have the same power to calculate a cross-sample mean and standard deviation
def get_all_samples(wildcards):
    files = []
    for s in list(samples['sample']):
        files.append(f'analysis/means/{s}_{wildcards.power}.bed.gz')
    files.sort()
    return files

# Rule to calculate the cross-sample mean and standard deviation for each power value
rule sample_mean_std_dev:
    input:
        files = get_all_samples,
    output:
        mean_gz = f"analysis/stats/cross_sample_mean_{{power}}.bed.gz",
        stdv_gz = f"analysis/stats/cross_sample_stddev_{{power}}.bed.gz",
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['mid_memory_gb'],
        walltime = config['hpc_parameters']['low_walltime'],
    params:
        outdir = directory("analysis/stats/"),
        sample = list(samples['sample']),
    script:
        "scripts/calculate_stats.py"

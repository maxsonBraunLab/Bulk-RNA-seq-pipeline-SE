__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""
"""Edits by Garth Kong at OHSU"""


import datetime
import sys
import os
import pandas as pd
import json


timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']


with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))


for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


def format_plot_columns():
    factors = config['meta_columns_to_plot'].keys()
    reformat_factors = '"' + '","'.join(factors) + '"'
    return 'c({})'.format(reformat_factors)


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")


rule all:
    input:
        # read alignment and QC -------------------------------------------------------------------
        expand("samples/fastp/{sample}_{dir}.fastq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("samples/fastqscreen/{sample}/{sample}_{dir}_screen.{fastqscreen_ext}", sample=SAMPLES, dir = ["R1", "R2"], fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        # expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        # expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        expand("results/diffexp/pairwise/{contrast}.pca_plot.pdf", contrast = config["diffexp"]["contrasts"]),
        "results/multiqc_report/multiqc_report.html",
        # DESeq2 group ----------------------------------------------------------------------------
        "data/{project_id}_norm.txt".format(project_id = config["project_id"]),
        "data/{project_id}_log2norm.txt".format(project_id = config["project_id"]),
        "results/diffexp/group/LRT_pca.pdf",
        "results/diffexp/group/MDS_table.txt",
        "results/diffexp/group/LRT_density_plot.pdf",
        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot.html"],contrast = config["diffexp"]["contrasts"]),
        "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
        # DESeq2 pairwise -------------------------------------------------------------------------
        expand(["results/diffexp/pairwise/{contrast}.qplot.pdf","results/diffexp/pairwise/{contrast}.qhist.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv"],contrast=config["diffexp"]["contrasts"]),
        expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
        expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
        # expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf", contrast = config["diffexp"]["contrasts"]),
        expand("samples/bigwig/{sample}.bw", sample = SAMPLES)

include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"

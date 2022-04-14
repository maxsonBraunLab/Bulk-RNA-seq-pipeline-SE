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

SAMPLES, = glob_wildcards("samples/raw/{sample}.fastq.gz")

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','pdf','xls']
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
        expand("samples/fastp/{sample}.fastq.gz", sample = SAMPLES),
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_fastqc.zip", sample = SAMPLES),
        expand("samples/fastqscreen/{sample}/{sample}_screen.{fastqscreen_ext}", sample=SAMPLES, fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
        expand("samples/bigwig/{sample}.bw", sample = SAMPLES),
        "results/multiqc_report/multiqc_report.html",
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"]),
        "results/multiqc_report/multiqc_report.html",

        # DESeq2 ----------------------------------------------------------------------------------
        "data/{project_id}_norm.txt".format(project_id=config["project_id"]),
        "results/diffexp/group/MDS_plot.pdf",

        # group analysis
        "results/diffexp/group/LRT_pca.pdf",

        # pairwise analysis
        expand("results/diffexp/pairwise/{contrast}.pca_plot.pdf", contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/pairwise/{contrast}.qplot.pdf", contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/pairwise/{{contrast}}.diffexp.{adjp}.VolcanoPlot.pdf".format(adjp=config["adjp"]), contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"]), contrast = config["diffexp"]["contrasts"])

include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"

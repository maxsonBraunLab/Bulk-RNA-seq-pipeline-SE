rule fastp:
    input:
        fwd = "samples/raw/{sample}_R1.fastq.gz",
        rev = "samples/raw/{sample}_R2.fastq.gz"
    output:
        fwd = "samples/fastp/{sample}_R1.fastq.gz",
        rev = "samples/fastp/{sample}_R2.fastq.gz"
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{sample}.fastp.json"
    threads: 8
    shell:
        "fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev} "
        "--detect_adapter_for_pe --thread {threads} -j {log} -h /dev/null"

rule fastqscreen:
    input:
        fwd = "samples/fastp/{sample}_R1.fastq.gz",
        rev = "samples/fastp/{sample}_R2.fastq.gz"
    output:
        "samples/fastqscreen/{sample}/{sample}_R1_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R1_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R1_screen.txt",
        "samples/fastqscreen/{sample}/{sample}_R2_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R2_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R2_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input.fwd} {input.rev}"""


rule fastqc:
    input:
        fwd = "samples/fastp/{sample}_R1.fastq.gz",
        rev = "samples/fastp/{sample}_R2.fastq.gz"
    output:
        fwd = "samples/fastqc/{sample}/{sample}_R1_fastqc.zip",
        rev = "samples/fastqc/{sample}/{sample}_R2_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """fastqc --outdir samples/fastqc/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""

rule STAR:
    input:
        fwd = "samples/fastp/{sample}_R1.fastq.gz",
        rev = "samples/fastp/{sample}_R2.fastq.gz"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out"
    threads: 12
    params:
        gtf=config["gtf_file"],
        genome_index=config["star_index"]
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.genome_index} \
        --readFilesIn {input.fwd} {input.rev} \
        --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
        --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
        --sjdbGTFtagExonParentGene gene_name \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --twopassMode Basic"
# install star via conda rather than use CEDAR version. 

rule index:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/samtools_env.yaml"
    shell:
        """samtools index {input} {output}"""

rule bigwig:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    output:
        "samples/bigwig/{sample}.bw"
    conda:
        "../envs/deeptools.yaml"
    threads: 12
    shell:
        "bamCoverage -b {input[0]} -o {output} -p {threads} --normalizeUsing CPM --binSize 10"

rule star_statistics:
    input:
        expand("samples/star/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"


rule compile_star_counts:
    input:
        expand("samples/star/{sample}_bam/ReadsPerGene.out.tab",sample=SAMPLES)
    params:
        samples=SAMPLES
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_star_counts.py"

rule filter_counts:
    input:
        countsFile="data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    params:
        anno=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito']
    script:
        "../scripts/RNAseq_filterCounts.R"

rule multiqc:
    input:
        expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES)
    output:
        "results/multiqc_report/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc logs/ samples/ -f -o results/multiqc_report"
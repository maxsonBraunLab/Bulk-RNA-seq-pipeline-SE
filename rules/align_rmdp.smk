rule fastp:
    input:
        "samples/raw/{sample}.fastq.gz"
    output:
        "samples/fastp/{sample}.fastq.gz"
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{sample}.fastp.json"
    threads: 8
    shell:
        "fastp -i {input} -o {output} --thread {threads} -j {log} -h /dev/null"

rule fastqscreen:
    input:
        "samples/fastp/{sample}.fastq.gz"
    output:
        "samples/fastqscreen/{sample}/{sample}_screen.html",
        "samples/fastqscreen/{sample}/{sample}_screen.png",
        "samples/fastqscreen/{sample}/{sample}_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        "fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input}"

rule fastqc:
    input:
        "samples/fastp/{sample}.fastq.gz"
    output:
        "samples/fastqc/{sample}/{sample}_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    message:
        "--- Quality check of raw data with Fastqc."
    shell:
        "fastqc --outdir samples/fastqc/{wildcards.sample} --extract -f fastq {input}"

rule STAR:
    input:
        "samples/fastp/{sample}.fastq.gz"
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
        --readFilesIn {input} \
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
        "samtools index {input} {output}"

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
        "bamCoverage -b {input[0]} -o {output} -p {threads} --normalizeUsing CPM --binSize 10 --smoothLength 50"

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

rule insertion_profile:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/rseqc/insertion_profile/{sample}.insertion_profile.r",
        "samples/rseqc/insertion_profile/{sample}.insertion_profile.pdf",
        "samples/rseqc/insertion_profile/{sample}.insertion_profile.xls",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "insertion_profile.py -i {input} -o samples/rseqc/insertion_profile/{wildcards.sample} -s 'SE'"

rule clipping_profile:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/rseqc/clipping_profile/{sample}.clipping_profile.r",
        "samples/rseqc/clipping_profile/{sample}.clipping_profile.pdf",
        "samples/rseqc/clipping_profile/{sample}.clipping_profile.xls",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "clipping_profile.py -i {input} -s 'SE' -o samples/rseqc/clipping_profile/{wildcards.sample}"

rule read_distribution:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    params:
        bed=config['bed_file']
    output:
        "samples/rseqc/read_distribution/{sample}.read_distribution.txt",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -i {input} -r {params.bed} > {output}"

rule read_GC:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/rseqc/read_GC/{sample}.GC.xls",
        "samples/rseqc/read_GC/{sample}.GC_plot.r",
        "samples/rseqc/read_GC/{sample}.GC_plot.pdf",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o samples/rseqc/read_GC/{wildcards.sample}"

rule multiqc:
    input:
        expand("samples/fastp/{sample}.fastq.gz", sample = SAMPLES),
        expand("samples/fastqscreen/{sample}/{sample}_screen.txt", sample = SAMPLES),
        expand("samples/fastqc/{sample}/{sample}_fastqc.zip", sample = SAMPLES),
        expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        expand("samples/rseqc/insertion_profile/{sample}.insertion_profile.xls", sample = SAMPLES),
        expand("samples/rseqc/clipping_profile/{sample}.clipping_profile.xls", sample = SAMPLES),
        expand("samples/rseqc/read_distribution/{sample}.read_distribution.txt", sample = SAMPLES),
        expand("samples/rseqc/read_GC/{sample}.GC.xls", sample = SAMPLES)
    output:
        "results/multiqc_report/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc samples/ -f -o results/multiqc_report"

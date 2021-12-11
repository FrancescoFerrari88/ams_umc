import os

rule SambambaExonCoverage:
    input:
        bam = "Bowtie2/{sample}.sorted.bam",
        bai = "Bowtie2/{sample}.sorted.bam.bai"
    output:
        "ExonGeneCoverage/{sample}_exonCov.txt"
    log: "ExonGeneCoverage/logs/{sample}_exonCov.log"
    params:
        extended_exon_bed=extended_exon_bed
    threads: lambda wildcards: 8 if 8 < maxJobs else maxJobs
    shell:"""
        sambamba depth region \
        -L {params.extended_exon_bed} -t {threads} {input.bam} > {output} 2> {log}
    """

rule SambambaPerBaseExonCoverage:
    input:
        bam = "Bowtie2/{sample}.sorted.bam",
        bai = "Bowtie2/{sample}.sorted.bam.bai"
    output:
        temp("ExonGeneCoverage/{sample}_perBaseExonCov.txt")
    log: "ExonGeneCoverage/logs/{sample}_perBaseExonCov.log"
    params:
        extended_exon_bed=extended_exon_bed
    threads: lambda wildcards: 8 if 8 < maxJobs else maxJobs
    shell:"""
        sambamba depth base -L {params.extended_exon_bed} -t {threads} \
        --min-coverage=1 {input.bam} \
        | awk -v OFS="\t" '{{print $1,$2,$2,$10,$3}}' \
        | tail -n +2 \
        | bedtools intersect -a - -b {params.extended_exon_bed} -wao \
        | awk -v OFS="\t" '{{print $1,$2,$3,$9,$5}}' > {output} 2> {log}
    """

rule AggGeneCoverage:
    input:
        "ExonGeneCoverage/{sample}_perBaseExonCov.txt"
    output:
        "ExonGeneCoverage/{sample}_aggGeneCov.txt"
    log: "ExonGeneCoverage/logs/{sample}_aggGeneCov.log"
    params:
        script_path=os.path.join(maindir,'scripts','AggGeneCoverage.py'),
        extended_exon_bed=extended_exon_bed
    shell:"""
        python {params.script_path} {input} {params.extended_exon_bed} 2> {log}
    """

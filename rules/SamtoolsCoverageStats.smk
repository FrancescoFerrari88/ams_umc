rule SamtoolsDepth:
    input:
        "CRAM/{sample}.cram"
    output:
        temp("CoverageStats/{sample}.depth")
    log: "CoverageStats/log/{sample}.depth.log"
    threads: maxJobs
    shell:"""
        samtools depth {input} > {output} 2> {log}
    """

rule SamtoolsCoverageStats:
    input:
        "CoverageStats/{sample}.depth"
    output:
        "CoverageStats/{sample}.stats"
    log: 
        mean_cov = "CoverageStats/log/{sample}.meanCov.log",
        median_cov = "CoverageStats/log/{sample}.medianCov.log",
        pct30 = "CoverageStats/log/{sample}.pct30.log"
    params:
        median_cov_script=os.path.join(maindir,"scripts","medianCoverage.py"),
        genome_length=genome_length,
        min_coverage_depth=30
    threads: maxJobs
    shell:"""
        cat {input} | awk '{{sum+=$3}} END {{ print "Mean Coverage = ",sum/{params.genome_length} }}' >> {output} 2> {log.mean_cov};
        python {params.median_cov_script} {input} {params.genome_length} >> {output} 2> {log.median_cov};
        cat {input} | awk '{{if($3>{params.min_coverage_depth}){{sum+=1}}}} END {{ print "PCT 30X = ",sum/{params.genome_length} }}' >> {output} 2> {log.pct30};
    """
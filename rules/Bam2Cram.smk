rule Bam2Cram:
    input:
        aligner+"/{sample}.sorted.bam"
    output:
        "CRAM/{sample}.cram"
    log: "CRAM/logs/{sample}.cram.log"
    params:
        fasta_genome=fasta_genome
    shell: """
            samtools view -C -T {params.fasta_genome} -o {output} {input}
        """
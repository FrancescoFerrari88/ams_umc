import os

rule origFASTQ1:
      input:
          os.path.join(indir,"{sample}"+reads[0]+ext)
      output:
          "originalFASTQ/{sample}"+reads[0]+".fastq.gz"
      params:
            cmd = lambda wildcards, input, output: "ln -s {} {}".format(input[0],output[0])
      shell: """
               {params.cmd}
          """

if pairedEnd:
    rule origFASTQ2:
        input:
            os.path.join(indir,"{sample}"+reads[1]+ext)
        output:
            "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
        params:
            cmd = lambda wildcards, input, output: "ln -s {} {}".format(input[0],output[0])
        shell: """
               {params.cmd}
            """
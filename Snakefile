import os
import common_functions as cf

### snakemake_workflows initialization ########################################
maindir = workflow.basedir

# load config file
globals().update(cf.load_configfile(workflow.overwrite_configfiles[0], config["verbose"]))

# load organism-specific data, i.e. genome indices, annotation, etc.
globals().update(cf.load_organism_data(genome, maindir, config["verbose"]))

# do workflow specific stuff
include: os.path.join(workflow.basedir, "internals.smk")

# import rules

### execute before workflow starts #############################################
################################################################################
onstart:
    if "verbose" in config and config["verbose"]:
        print("--- Workflow parameters --------------------------------------------------------")
        print("samples:", samples)
        print("paired:", pairedEnd)
        print("read extension:", reads)
        print("maximum insert size (Bowtie2 -X):", insertSizeMax)
        print("-" * 80, "\n")

        print("--- Environment ----------------------------------------------------------------")
        print("$TMPDIR: ",os.getenv('TMPDIR', ""))
        print("-" * 80, "\n")

### main rule ##################################################################
################################################################################


### execute after workflow finished ############################################
################################################################################
onsuccess:
    if "verbose" in config and config["verbose"]:
        print("\n--- DNA mapping workflow finished successfully! --------------------------------\n")

onerror:
    print("\n !!! ERROR in DNA mapping workflow! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")


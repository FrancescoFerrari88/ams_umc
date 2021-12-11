import glob
import os

### Initialization #############################################################

### Check if input files exist, get sample names and idxRange (single/paired end)
infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*'+ext)))

if infiles == []:
    sys.exit("Error! Samples extension in {} are not {}. "
             "Please change the extensions to it or update the config.yaml file "
             "with your desired extension.".format(indir,ext))
samples = cf.get_sample_names(infiles,ext,reads)

pairedEnd = cf.is_paired(infiles,ext,reads)

del infiles

if not samples:
    sys.exit("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")

idxRange = 1
if pairedEnd:
    idxRange = 2

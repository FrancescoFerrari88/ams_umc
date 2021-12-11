#!/usr/bin/env python3

import sys
import os
import pandas as pd

def get_sample_outdir_info(cov_infile_path):
    outdir = os.path.dirname(cov_infile_path)
    sample_ext = os.path.basename(cov_infile_path)
    sample = "_".join(sample_ext.split("_")[:-1])
    return sample, outdir

def main():

    per_base_exon_coverage_file = sys.argv[1]
    extended_exons_file = sys.argv[2]

    # import files into pandas df
    df_cov = pd.read_csv(per_base_exon_coverage_file, header=None, sep="\t")
    df_cov.columns = ["chr","start","end","gene","cov"]

    df_exons = pd.read_csv(extended_exons_file, header=None, sep="\t")
    df_exons.columns = ["chr","start","end","gene"]

    # compute sum of extended exons per gene
    df_exons["exon_len"] = df_exons["end"] - df_exons["start"]
    agg_exon_len = df_exons[["gene","exon_len"]].groupby("gene").sum()["exon_len"]
    print(agg_exon_len.shape[0])

    # compute sum of coverage per gene
    sum_cov = df_cov[["gene","cov"]].groupby("gene").sum()["cov"]
    print(sum_cov.shape[0])

    # merge on gene and compute mean coverage per gene based on extended exons
    merge = sum_cov.to_frame().merge(agg_exon_len.to_frame(), 
    how="right", right_index=True, left_index=True)
    merge.fillna(0, axis=1, inplace=True)
    
    merge["mean_per_base_gene_cov"] = merge["cov"] / merge["exon_len"]

    # write to file
    sample, dirname = get_sample_outdir_info(per_base_exon_coverage_file)

    outfile_path_exonCov = os.path.join(dirname,f"{sample}_aggGeneCov.txt")
    merge[["mean_per_base_gene_cov"]].to_csv(outfile_path_exonCov, sep="\t", index=True)

    
if __name__ == '__main__':
    main()
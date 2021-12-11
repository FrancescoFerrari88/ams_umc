#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def main():

    # get command line arguments
    covFile = sys.argv[1]
    genome_len = int(sys.argv[2])

    # read in depth file
    df = pd.read_csv(covFile, sep="\t", header=None)
    df.columns = ["chr","base","cov"]

    # calculate uncovered bases
    n_zeros = int(genome_len - df.shape[0])

    # get full coverage vector
    full_vector_cov = np.concatenate([df["cov"].values, np.zeros(n_zeros,dtype=int)])

    # print median coverage
    print("Median Coverage =", np.median(full_vector_cov))


if __name__ == "__main__":
    main()
#!/usr/bin/python3

import csv
from collections import Counter
from Bio import SeqIO
import argparse


def run_wh_on_fasta(in_fa, out_csv):
    """
    :param in_fa: input fasta file
    :param out_csv: name of output csv file
    """

    with open(out_csv, "w") as f_out:
        fa = SeqIO.parse(in_fa, "fasta")
        wh_results = csv.writer(f_out)
        wh_results.writerow(['sid', 'WH.solubility'])
        for record in fa:
            seq = str(record.seq)
            wh_results.writerow([record.id, wh_model(seq)])


def wh_model(seq):
    """
    modified Wilkinson-Harrison model for prediction of protein solubility
    :param seq: amino acid sequence
    :return: solubility
    """
    lambda_1 = 15.43
    lambda_2 = -29.56
    discr_CV = 1.71

    # aminoacid/letter frequencies
    amino = Counter(seq.upper())
    # number of letter
    n = len(seq)
    # canonical variable
    CV = \
        lambda_1*((amino['N']+amino['G']+amino['P']+amino['S'])/n)\
        + \
        lambda_2*abs(((amino['R']+amino['K']) - (amino['D']+amino['E']))/n - 0.03)
    # Davis: "If CV-discr_CV is positive, the protein is predicted to be
    # insoluble, while if CV-discr_CV is negative, the protein is predicted to be
    # soluble."
    diff_CV = CV - discr_CV
    # Probability of result
    prob = 0.4934 + 0.276 * abs(diff_CV) - 0.0392 * (diff_CV*diff_CV)
    # probability is kind of "solubility" in this context
    if diff_CV > 0:
        solubility = 1 - prob
    else:
        solubility = prob
        # following corrections are based on script that is used in Loschmidt
        # laboratories
        if solubility < 0:
            solubility = 0.98
        elif diff_CV < -3.99:
            solubility = 0.98
    return solubility


def main():
    args = arguments()
    run_wh_on_fasta(args.i_fa, args.o_csv)


def arguments():

    parser = argparse.ArgumentParser(description="Modified WH model for "
                                                 "prediction of solubility")
    parser.add_argument('--i_fa', help="Input file in fasta format",
                        required=True)
    parser.add_argument('--o_csv', help="Name of output file", required=True)
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()












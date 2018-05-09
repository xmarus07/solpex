#!/usr/bin/python3

import csv
import json
import argparse
from tqdm import tqdm
import sys
import numpy as np

class EspritzProcessor:

    def __init__(self, input, type, round_to=4, progress=False):
        self.input = input
        self.type = type
        self.round_to = round_to
        self.progress = progress

    def _run_on_each(self, fun, out_csv):
        # increasing limit (some items are quite large)
        csv.field_size_limit(sys.maxsize)
        with open(self.input, "r") as f:
            # get number of lines for time estimation
            if self.progress:
                num_lines = sum(1 for line in f) - 1
                f.seek(0)
            fess_csv = csv.reader(f)
            # skipping & checking header
            assert next(fess_csv) == ['sid', 'espritz.json']
            if self.type == "mean":
                SID, JSON = 0, 1
                fess_iter = fess_csv
                if self.progress:
                    fess_iter = tqdm(fess_iter, total=num_lines)
                for row in fess_iter:
                    out_csv.writerow([row[SID]] + fun(row[JSON]))

    def _mean_energy(self, energy_json):
        energy = json.loads(energy_json)
        if len(energy) == 1:
            if len(energy[0]) == 0:
                return [np.nan]
        s = sum([i[0] for i in energy])
        return [round(s/len(energy), self.round_to)]

    def mean_energy(self, output):
        with open(output, "w") as f:
            out_csv = csv.writer(f)
            out_csv.writerow(['sid', 'mean'])
            self._run_on_each(self._mean_energy, out_csv)


def main():
    args = arguments()
    esp = EspritzProcessor(args.i_raw_csv, args.type, progress=True)
    if args.type == 'mean':
        esp.mean_energy(args.o_csv)


def arguments():

    parser = argparse.ArgumentParser(
        description=("Process results of espritz\n"
                     "Run example:\n"
                     "./process_espritz.py "
                     "--i_raw_csv ../../features/espritz_results.csv "
                     "--o_csv espritz_mean.csv --type mean"
                     ))
    parser.add_argument('--i_raw_csv', help="Input file - results of "
                                      "./run_espritz.py ", required=True)
    parser.add_argument('--o_csv', help="Name of output file", required=True)
    parser.add_argument('--type', help="Type of processing",
                        choices=['mean'], default='mean')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()

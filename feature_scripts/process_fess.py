#!/usr/bin/python3
import csv
import json
import argparse
import pandas
from tqdm import tqdm
import sys


class FessProcessor:

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
            assert next(fess_csv) == ['sid', 'fess.json.results']
            SID, JSON = 0, 1
            fess_iter = fess_csv
            if self.progress:
                fess_iter = tqdm(fess_iter, total=num_lines)
            if self.type == "ratio":
                for row in fess_iter:
                    out_csv.writerow([row[SID]] + fun(row[JSON]))

    @staticmethod
    def load_fess_json(fess_json):
        return pandas.DataFrame(json.loads(fess_json))

    def _count_ratio(self, fess_json):
        results = json.loads(fess_json)
        if results == [[]]:
            ratio = ["NA"]*3
        else:
            counts = {'H': 0, 'E': 0, 'C': 0}
            for aa in results:
                counts[aa[0]] += 1
            h, e, c = counts['H'], counts['E'], counts['C']
            ratio = [round(i/len(results), self.round_to) for i in (h, e, c)]
        return ratio

    def count_ratio(self, output):
        with open(output, "w") as f:
            out_csv = csv.writer(f)
            out_csv.writerow(['sid', 'H', 'E', 'C'])
            self._run_on_each(self._count_ratio, out_csv)


def main():
    args = arguments()
    fess = FessProcessor(args.i_raw_csv, args.type, progress=True)
    if args.type == 'ratio':
        fess.count_ratio(args.o_csv)


def arguments():

    parser = argparse.ArgumentParser(
        description=("Process results of fess.\n Run example: \n"
                     "./fess_process_results.py "
                     "--i_raw_csv ../../features/fess_results.csv "
                     "--o_csv fess_HEC_ratio.csv --type ratio"
                     ))
    parser.add_argument('--i_raw_csv', help="Input file - results of "
                                             "./run_fess.py script",
                                             required=True)
    parser.add_argument('--o_csv', help="Name of output file", required=True)
    parser.add_argument('--type', help="Type of processing", required=True,
                        choices=['ratio'])
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()




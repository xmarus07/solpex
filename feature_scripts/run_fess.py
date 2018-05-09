#!/usr/bin/python3

import multiprocessing as mp
import argparse
import subprocess
import tempfile
import json
import pandas as pd
import csv
from tqdm import tqdm
import shutil
# import from upper directory common
import os
from Bio import SeqIO
from common.seq_count_fa import seq_count_fa
from common.clear_dir import clear_dir

class FessExecutionFailed(Exception):
    pass

class FessRunner:

    __queue = None
    __exists = False

    def __init__(self, pools, fess_path, tmp_dir, round_to=4,
                 progress=False):
        if FessRunner.__exists:
            raise ValueError("Only one instance of runner is allowed.")
        self.pools = pools
        self.fess_path = fess_path
        self.fess_inst_dir = os.path.dirname(fess_path)
        self.round_to = round_to
        self.tmp_dir = tmp_dir
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        else:
            clear_dir(self.tmp_dir)
        self.progress = progress

    def __init_que(self):

        FessRunner.__queue = mp.Queue()
        for i in range(self.pools*2):
            FessRunner.__queue.put(i)

    def process_results(self, f_name):
        with open(f_name, 'r') as f:
            results = pd.read_csv(f, delimiter=r'\s+')
            results[['H', 'E', 'C']] = \
                results[['H', 'E', 'C']].apply(lambda x: round(x, self.round_to))
            return results.values.tolist()

    def _worker(self, seq_record):
        file_id = FessRunner.__queue.get()
        f_name_in = os.path.join(self.tmp_dir, str(file_id)+'.fasta')
        f_name_out = os.path.join(self.tmp_dir, str(file_id)+'.fess')

        SeqIO.write(seq_record, f_name_in, format='fasta')

        result = [[]]
        try:
            subprocess.run([self.fess_path, f_name_in, f_name_out], check=True,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT,
                           cwd=self.fess_inst_dir)
            result = self.process_results(f_name_out)
        except subprocess.CalledProcessError:
            raise FessExecutionFailed()
        finally:
            FessRunner.__queue.put(file_id)

        return [seq_record.id, json.dumps(result)]

    def run_fess(self, fasta, output):
        fa = SeqIO.parse(fasta, format='fasta')
        with open(output, "w") as f, \
                mp.Pool(self.pools, initializer=self.__init_que()) as pools:
            pool_iter = pools.imap(self._worker, fa)
            csv_out = csv.writer(f)
            csv_out.writerow(['sid', 'fess.json.results'])
            job_iter = pool_iter
            if self.progress:
                total_iter = seq_count_fa(fasta)
                job_iter = tqdm(job_iter, total=total_iter)
            for result in job_iter:
                csv_out.writerow(result)
        shutil.rmtree(self.tmp_dir)


def main():
    args = arguments()
    fa = os.path.abspath(args.i_fa)
    o = os.path.abspath(args.o_raw_csv)
    tmp_dir = os.path.join(os.path.dirname(__file__), 'tmp_fess_dir')
    tmp_dir = os.path.abspath(tmp_dir)
    esp_path = os.path.abspath('../additional_software/fess/fess')
    fess_runner = FessRunner(args.p, esp_path, tmp_dir, progress=True)
    fess_runner.run_fess(fa, o)


def arguments():

    parser = argparse.ArgumentParser(description="Run fess on fasta file that "
                                                 "contains more than one "
                                                 "sequence and"
                                                 "store results as csv + "
                                                 "json(value)")
    parser.add_argument('--i_fa', help="Input file in fasta format",
                        required=True,)
    parser.add_argument('--o_raw_csv', help="Name of output file",
                        required=True)
    parser.add_argument('-p', help="Number of pools", required=True,
                        type=int)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()


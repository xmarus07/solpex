#!/usr/bin/python3

import multiprocessing as mp
import argparse
import subprocess
import json
import pandas as pd
import csv
from tqdm import tqdm
import shutil
# import from upper directory common
import os
import math
from Bio import SeqIO
from common.FastaChunk import FastaChunk
from common.seq_count_fa import seq_count_fa
from common.clear_dir import clear_dir

class EspritzExecutionFailed(Exception):
    pass

class EspritzRunner:

    __queue = None
    __exists = False

    def __init__(self, pools, esp_path, dataset, p_decision,
                 tmp_dir, skip=9, round_to=4, progress=False):
        if EspritzRunner.__exists:
            raise ValueError("Only one instance of runner is allowed.")
        EspritzRunner.__exists = True
        self.pools = pools
        self.esp_path = esp_path
        self.esp_inst_dir = os.path.dirname(self.esp_path)
        self.p_decision = p_decision
        self.dataset = dataset
        self.round_to = round_to
        self.tmp_dir = tmp_dir
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        else:
            clear_dir(self.tmp_dir)
        self.skip = skip
        self.progress = progress

    def __init_que(self):
        EspritzRunner.__queue = mp.Queue()
        for i in range(self.pools*2):
            EspritzRunner.__queue.put(i)

    def process_results(self, f_names):
        chunk = list()
        for f_name in f_names:
            with open(f_name, 'r') as f:
                results = pd.read_csv(f, delimiter=r'\s+',
                                      skiprows=self.skip,
                                      names=["Order", "value"])
                results[['value']] = \
                    results[['value']].apply(lambda x: round(x, self.round_to))
                chunk.append(json.dumps(results[['value']].values.tolist()))
        return chunk

    def _worker(self, seq_records):
        dir_id = EspritzRunner.__queue.get()
        dir_name = os.path.join(self.tmp_dir, "tmp_" + str(dir_id))
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        out_fs = list()
        result = [[]]
        results = list()
        for record in seq_records:
            f_name_in = os.path.join(dir_name, record.id+'.fasta')
            SeqIO.write(record, f_name_in, format='fasta')
            f_name_out = os.path.join(dir_name, record.id+'.espritz')
            out_fs.append(f_name_out)
            results.append([record.id, result])
        try:
            subprocess.run([self.esp_path, dir_name, self.dataset,
                            self.p_decision],
                           check=True, cwd=self.esp_inst_dir,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT)
            values = self.process_results(out_fs)
            assert len(values) == len(results)
            for i in range(len(results)):
                results[i][1] = values[i]
        except subprocess.CalledProcessError:
            raise EspritzExecutionFailed()
        finally:
            EspritzRunner.__queue.put(dir_id)
            shutil.rmtree(dir_name)
        return results

    def run_espritz(self, fasta, output, chunk_size=100):
        fa = FastaChunk(chunk_size, fasta)
        with open(output, "w") as f, \
                mp.Pool(self.pools, initializer=self.__init_que()) as pools:
            pool_iter = pools.imap(self._worker, fa)
            csv_out = csv.writer(f)
            csv_out.writerow(['sid', 'espritz.json'])
            job_iter = pool_iter
            if self.progress:
                total_iter = math.ceil(seq_count_fa(fasta) / chunk_size)
                job_iter = tqdm(job_iter, total=total_iter)
            for result in job_iter:
                csv_out.writerows(result)
        shutil.rmtree(self.tmp_dir)


def main():
    args = arguments()
    fa = os.path.abspath(args.i_fa)
    o = os.path.abspath(args.o_csv)
    tmp_dir = os.path.join(os.path.dirname(__file__), 'tmp_esp_dir')
    tmp_dir = os.path.abspath(tmp_dir)
    esp_path = os.path.abspath('../additional_software/espritz/espritz.pl')
    esp_runner = EspritzRunner(args.p, esp_path, 'X', '1', tmp_dir, progress=True)
    esp_runner.run_espritz(fa, o)


def arguments():

    parser = argparse.ArgumentParser(description="Run espritz on fasta file "
                                                 "that "
                                                 "contains more than one "
                                                 "sequence and"
                                                 "store results as csv + "
                                                 "json(value)")
    parser.add_argument('--i_fa', help="Input file in fasta format",
                        required=True)
    parser.add_argument('--o_csv', help="Name of output file", required=True)
    parser.add_argument('-p', help="Number of pools", required=True,
                        type=int)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()

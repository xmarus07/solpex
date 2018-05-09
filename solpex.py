#!/usr/bin/python3
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import ProtParam
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from sklearn.externals import joblib
import argparse
import subprocess as sb
import pandas as pd
import numpy as np
import sys
import os

from feature_scripts.common.get_abs_path import get_abs_path

# for loading other modules inside directory feature_scripts
sys.path.insert(1, get_abs_path(__file__, './feature_scripts'))
sys.path.insert(1, get_abs_path(__file__, './data'))

from feature_scripts.dimers_comb import DimerComb
from feature_scripts.common.prefix import set_df_prefix
from feature_scripts.physico_chemical import fracnumcharge, kr_ratio
from feature_scripts.blast6_to_max_id_csv import process_blast6
from feature_scripts.run_espritz import EspritzRunner
from feature_scripts.run_espritz import EspritzExecutionFailed
from feature_scripts.process_espritz import EspritzProcessor
from feature_scripts.run_fess import FessRunner
from feature_scripts.run_fess import FessExecutionFailed
from feature_scripts.process_fess import FessProcessor
from data.RandomForestModel import RandomForestModel

AA = IUPAC.protein.letters

class Paths:
    """
    Relative paths to solpex scripts
    """
    _SCRIPT = __file__
    _MODEL = './data/rf_model.pkl'
    _USEARCH = './additional_software/usearch/usearch10.0.240_i86linux32'
    _PDB_ECOLI_FA = './data/Ecoli_xray_nmr_pdb.fa'
    _TMHMM = './additional_software/tmhmm-2.0c.Linux/tmhmm-2.0c/bin/tmhmm'
    _ESPRITZ = './additional_software/espritz/espritz.pl'
    _FESS = './additional_software/fess/fess'

    @staticmethod
    def _get_abs_path(rel_path):
        return get_abs_path(Paths._SCRIPT, rel_path)

    @staticmethod
    def get_and_check_abs_file_path(rel_path, exp):
        abs_path = get_abs_path(Paths._SCRIPT, rel_path)
        if not os.path.exists(abs_path) or os.path.isdir(abs_path):
            raise exp
        return abs_path

    @staticmethod
    def model():
        return Paths.get_and_check_abs_file_path(Paths._MODEL,
                                                 ModelInvalidPath())

    @staticmethod
    def usearch():
        return Paths.get_and_check_abs_file_path(Paths._USEARCH,
                                                 UsearchInvalidPath())

    @staticmethod
    def pdb_db():
        return Paths.get_and_check_abs_file_path(Paths._PDB_ECOLI_FA,
                                                 PdbDatabaseNotFound())

    @staticmethod
    def tmhmm():
        return Paths.get_and_check_abs_file_path(Paths._TMHMM,
                                                 TmhmmInvalidPath())

    @staticmethod
    def espritz():
        return Paths.get_and_check_abs_file_path(Paths._ESPRITZ,
                                                 EspritzInvalidPath())
    @staticmethod
    def fess():
        return Paths.get_and_check_abs_file_path(Paths._FESS,
                                                 FessInvalidPath())

class InvalidAlphabet(Exception):
    pass

class ShortSequence(Exception):
    pass

class DuplicatedSid(Exception):
    pass

class UsearchInvalidPath(Exception):
    pass

class UsearchExcecutionFailed(Exception):
    pass

class PdbDatabaseNotFound(Exception):
    pass

class TmhmmInvalidPath(Exception):
    pass

class TmhmmExecutionFailed(Exception):
    pass

class TmhmmParsingError(Exception):
    pass

class EspritzInvalidPath(Exception):
    pass

class EspritzProcessingFailed(Exception):
    pass

class FessInvalidPath(Exception):
    pass

class FessProcessingError(Exception):
    pass

class ModelInvalidPath(Exception):
    pass

class ModelIsNotCompatible(Exception):
    pass

class MissingModelFeatures(Exception):
    pass

class Predictor:

    _PRE_MONOMERS = "monomers"
    _PRE_DIMERS = "dimers_comb"
    _PRE_PHYSICO_CHEM = "physico_chemical"
    _PRE_IDENTITY = "ecoli_usearch_identity"
    _PRE_TMHMM = "tmhmm"
    _PRE_ESPRITZ = "espritz_mean"
    _PRE_FESS = "fess_HEC_ratio"

    _THREADS = 3
    _MIN_SEQ_LENGTH = 20

    def __init__(self, fasta_file, tmp_dir, model_path, usearch, pdb_db,
                 tmhmm, fess, espritz, usearch_threads=_THREADS,
                 espritz_threads=_THREADS, fess_threads=_THREADS
                 ):

        model_path = model_path
        try:
            self.model = joblib.load(model_path)
        except AttributeError:
            raise ModelIsNotCompatible()
        if type(self.model) is not RandomForestModel:
            raise ModelIsNotCompatible()
        # loading fasta to dataframe
        fasta_parser = SeqIO.parse(fasta_file, "fasta")
        fa_id = []
        sequences = []
        for record in fasta_parser:
            fa_id.append(record.id)
            seq = str(record.seq)
            if len(seq) < Predictor._MIN_SEQ_LENGTH:
                raise ShortSequence()
            for aa in seq:
                if aa not in AA:
                    raise InvalidAlphabet()
            sequences.append(seq)

        self.seq = pd.DataFrame({"sequence": sequences, "fa_id": fa_id},
                                index=range(len(sequences)))
        self.seq.index = self.seq.index.astype(str)
        self.seq.index.name = "sid"
        if self.seq.index.nunique() != self.seq.index.shape[0]:
            raise DuplicatedSid()
        self.features = pd.DataFrame(index=self.seq.index)

        # other arguments
        self.tmp_dir = os.path.abspath(tmp_dir)
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        self.fasta_path = None
        self.usearch_threads = usearch_threads
        self.usearch = usearch
        self.pdb_db = pdb_db
        self.tmhmm = tmhmm
        self.espritz_threads = espritz_threads
        self.espritz = espritz
        self.fess_threads = fess_threads
        self.fess = fess

    def file_path(self, file_name):
        return os.path.join(self.tmp_dir, file_name)

    def compute_features(self):
        self.create_fasta("query.fa")
        self._add_monomers()
        self._add_dimers()
        self._add_physico_chemical()
        self._add_tmhmm()
        self._add_fess()
        self._add_espritz()
        self._add_usearch_identity()

    def predict(self, round_to=4):
        if len(self.features.columns) != len(self.model.order):
            raise MissingModelFeatures()
        is_null = self.features.isnull().any(axis=1)
        null_features = self.features[is_null]
        for index, row in null_features.iterrows():
            for col in null_features.columns:
                if pd.isnull(row[col]):
                    print("Warning feature {f} can not be calculated for "
                          "sequence with id {id}, mean of training set will be "
                          "used.".format(f=col, id=index), file=sys.stderr)
                    self.features.at[index, col] = self.model.features_mean[col]

        pred = self.model.predict(self.features)
        pred = np.round(pred, round_to)
        results = pd.DataFrame({"soluble":pred}, index=self.features.index)
        results = results.join(self.seq["fa_id"])
        results.index.name = "runtime_id"
        return results[["fa_id", "soluble"]]

    def _join(self, feature:pd.DataFrame, prefix):
        feature.index = feature.index.astype(str)
        set_df_prefix(feature, prefix)
        cols = feature.columns[feature.columns.isin(self.model.order)]
        feature = feature[cols]
        self.features = self.features.join(feature, how="left")

    def _add_monomers(self):
        print("Computing monomers")
        monomers = pd.DataFrame(index=self.features.index, columns=list(AA))
        for index, row in self.seq.iterrows():
            analysis = ProtParam.ProteinAnalysis(row["sequence"])
            aa_freq = analysis.get_amino_acids_percent()
            monomers.loc[index] = aa_freq
        self._join(monomers, Predictor._PRE_MONOMERS)

    def _add_dimers(self):
        print("Computing dimers")
        dim = DimerComb()
        dimers_comb = pd.DataFrame(index=self.features.index, columns=dim.combs)
        for index, row in self.seq.iterrows():
            dimers_comb.loc[index] = dim.get_comb_ratio(row["sequence"])
        self._join(dimers_comb, Predictor._PRE_DIMERS)

    def _add_physico_chemical(self):
        print("Computing physico-chemical features")
        cols = ["fracnumcharge", "kr_ratio", "aa_helix", "aa_sheet",
                "aa_turn", "molecular_weight", "aromaticity",
                "avg_molecular_weight", "flexibility", "gravy",
                "isoelectric_point"]
        physico_chem = pd.DataFrame(index=self.features.index, columns=cols)
        for index, row in self.seq.iterrows():
            pc_row = dict.fromkeys(cols, np.nan)
            analysis = ProtParam.ProteinAnalysis(row["sequence"])
            aa_freq = analysis.get_amino_acids_percent()
            pc_row["fracnumcharge"] = fracnumcharge(aa_freq)
            pc_row["kr_ratio"] = kr_ratio(aa_freq)
            h, s, t = analysis.secondary_structure_fraction()
            pc_row["aa_helix"] = h
            pc_row["aa_sheet"] = s
            pc_row["aa_turn"] = t
            pc_row['molecular_weight'] = analysis.molecular_weight()
            pc_row['length'] = analysis.length
            pc_row['avg_molecular_weight'] = pc_row['molecular_weight']/pc_row['length']
            pc_row['aromaticity'] = analysis.aromaticity()
            pc_row['flexibility'] = np.mean(analysis.flexibility())
            pc_row['gravy'] = analysis.gravy()
            pc_row['isoelectric_point'] = analysis.isoelectric_point()
            physico_chem.loc[index] = pc_row
        self._join(physico_chem, Predictor._PRE_PHYSICO_CHEM)

    def _add_usearch_identity(self, b6="identity.b6"):
        print("Computing identity")
        b6_path = self.file_path(b6)
        check_remove_file(b6_path)
        usearch_arguments = ['-search_global', self.fasta_path,
                             '-db', self.pdb_db, '-id', '0.0', '-blast6out',
                             b6_path, '-threads', str(self.usearch_threads),
                             '-top_hits_only']
        try:
            sb.run([self.usearch] + usearch_arguments, check=True,
                   stdout=sb.DEVNULL)
        except sb.CalledProcessError:
            raise UsearchExcecutionFailed()
        if not os.path.exists(b6_path):
            raise UsearchExcecutionFailed()
        identity = process_blast6(b6_path, "sid", "identity")
        identity.set_index('sid', inplace=True)
        self._join(identity, Predictor._PRE_IDENTITY)

    def _add_tmhmm(self, tmhmm="tmhmm.tmhmm"):
        print("Computing transmembrane regions")
        tmhmm_path = self.file_path(tmhmm)
        check_remove_file(tmhmm_path)
        try:
            with open(tmhmm_path, "w") as tm_f:
                sb.run([self.tmhmm, self.fasta_path], check=True, stdout=tm_f)
        except sb.CalledProcessError:
            raise TmhmmExecutionFailed()
        if not os.path.exists(tmhmm_path):
            raise TmhmmExecutionFailed()
        tm_df = tmhmm_to_df(tmhmm_path, "sid")
        tm_df.set_index('sid', inplace=True)
        self._join(tm_df, Predictor._PRE_TMHMM)

    def _add_espritz(self, esp_raw_out="espritz_raw.csv",
                     esp_mean_path="espritz_mean.csv", tmp_esp="tmp_espritz"):
        print("Computing disorder")
        esp_raw_out = self.file_path(esp_raw_out)
        check_remove_file(esp_raw_out)
        tmp_esp = self.file_path(tmp_esp)
        # running espritz
        esp_runner = EspritzRunner(self.espritz_threads, self.espritz,
                                   'X', '1', tmp_esp)
        esp_runner.run_espritz(self.fasta_path, esp_raw_out)
        if not os.path.exists(esp_raw_out):
            raise EspritzExecutionFailed()
        # processing espritz
        esp_proc = EspritzProcessor(esp_raw_out, 'mean')
        esp_mean_path = self.file_path(esp_mean_path)
        check_remove_file(esp_mean_path)
        esp_proc.mean_energy(esp_mean_path)
        if not os.path.exists(esp_mean_path):
            raise EspritzProcessingFailed()
        esp_mean = pd.read_csv(esp_mean_path)
        esp_mean.set_index('sid', inplace=True)
        self._join(esp_mean, Predictor._PRE_ESPRITZ)

    def _add_fess(self, fess_raw_out="fess_raw.csv",
                  fess_hec_ratio_path="fess_ratio.csv", tmp_fess="tmp_fess"):
        print("Computing secondary structure")
        fess_raw_out = self.file_path(fess_raw_out)
        check_remove_file(fess_raw_out)
        tmp_fess = self.file_path(tmp_fess)
        # running fess
        fess_runner = FessRunner(self.fess_threads, self.fess, tmp_fess)
        fess_runner.run_fess(self.fasta_path, fess_raw_out)
        if not os.path.exists(fess_raw_out):
            raise FessExecutionFailed()
        fess_hec_ratio_path = self.file_path(fess_hec_ratio_path)
        check_remove_file(fess_hec_ratio_path)
        fess_proc = FessProcessor(fess_raw_out, 'ratio')
        fess_proc.count_ratio(fess_hec_ratio_path)
        if not os.path.exists(fess_hec_ratio_path):
            raise FessProcessingError()
        fess_hec_ratio = pd.read_csv(fess_hec_ratio_path)
        fess_hec_ratio.set_index('sid', inplace=True)
        self._join(fess_hec_ratio, Predictor._PRE_FESS)

    def create_fasta(self, file_name):
        file_path = self.file_path(file_name)
        with open(file_path, "w") as f_fasta:
            fasta_wr = FastaIO.FastaWriter(f_fasta, wrap=None)
            fasta_wr.write_header()  # does nothing, but is required
            for index, row in self.seq.iterrows():
                record = SeqIO.SeqRecord(
                    seq=Seq(row["sequence"], IUPAC.protein),
                    id=str(index), description='')
                fasta_wr.write_record(record)
            fasta_wr.write_footer()  # does nothing, but is required
        self.fasta_path = file_path


def check_remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)


def tmhmm_to_df(tmhmm_path, id_col):
        F_COUNT = 5
        df_dict = {id_col:[], "len":[], "ExpAA":[], "First60":[],
                   "PredHel":[], "Topology":[]}
        with open(tmhmm_path, "r") as tm_f:
            for line in tm_f:
                line = line.rstrip()
                values = line.split('\t')
                seq_id, features = values[0], values[1:]
                df_dict[id_col].append(seq_id)
                for f in features:
                    key, val = f.split('=')
                    df_dict[key].append(val)
                if len(features) != F_COUNT:
                    raise TmhmmParsingError()
        df = pd.DataFrame(df_dict)
        df.rename(columns={"ExpAA": "exp_aa", "First60": "first_60",
                           "PredHel": "pred_hel", "Topology": "topology"},
                  inplace=True)
        return df


def main():
    args = arguments()
    if not os.path.isfile(args.i_fa):
        print("Invalid path to input FASTA file.", file=sys.stderr)
        return 1
    Paths._MODEL = args.model
    Paths._USEARCH = args.usearch
    Paths._ESPRITZ = args.espritz
    Paths._FESS = args.fess
    Paths._TMHMM = args.tmhmm
    Paths._PDB = args.pdb
    try:
        pred = Predictor(args.i_fa, args.tmp_dir, model_path=Paths.model(),
                       usearch=Paths.usearch(), espritz=Paths.espritz(),
                       fess=Paths.fess(), tmhmm=Paths.tmhmm(),
                       pdb_db=Paths.pdb_db())
        pred.compute_features()
        res = pred.predict()
        res.to_csv(args.o_csv)
    except UsearchInvalidPath:
        print("Path to USEARCH is invalid:", Paths._USEARCH, file=sys.stderr)
    except EspritzInvalidPath:
        print("Path to Espritz is invalid:", Paths._ESPRITZ, file=sys.stderr)
    except FessInvalidPath:
        print("Path to Fess is invalid:", Paths._FESS, file=sys.stderr)
    except TmhmmInvalidPath:
        print("Path to TMHMM is invalid:",Paths._TMHMM ,file=sys.stderr)
    except PdbDatabaseNotFound:
        print("PDB database was not found on given path:", Paths._PDB,
              file=sys.stderr)
    except ModelInvalidPath:
        print("Model does not exists, path:", Paths._MODEL, file=sys.stderr)
    except InvalidAlphabet:
        print("Invalid amino acid alphabet, sequences can contain only "
              "standard amino acids.", file=sys.stderr)
    except ShortSequence:
        print("Some sequences are too short, minimum length of sequence is 20 "
              "amino acids.", file=sys.stderr)
    except DuplicatedSid:
        print("Duplicated identifier in FASTA file, each sequence must "
              "contain unique identifier, duplicated sequences with same "
              "identifier are not allowed.", file=sys.stderr)
    except UsearchExcecutionFailed:
        print("Execution of USEARCH failed.", file=sys.stderr)
    except TmhmmExecutionFailed:
        print("Execution of TMHMM failed.", file=sys.stderr)
    except EspritzExecutionFailed:
        print("Execution of Espritz failed.", file=sys.stderr)
    except FessExecutionFailed:
        print("Execution of Fess failed.", file=sys.stderr)
    except TmhmmParsingError:
        print("Processing of TMHMM results failed.", file=sys.stderr)
    except EspritzProcessingFailed:
        print("Processing of Espritz results failed.", file=sys.stderr)
    except FessProcessingError:
        print("Processing of Fess results failed.", file=sys.stderr)
    except ModelIsNotCompatible:
        print("Loaded model is not compatible with predictor.", file=sys.stderr)
    except MissingModelFeatures:
        print("Model requires features that are not computed by predictor.",
              file=sys.stderr)


def arguments():
    parser = argparse.ArgumentParser(
        description="Protein solubility predictor Solpex. "
                    "NOTE: All paths of external tools "
                    "should be either relative to position of solpex.py "
                    "script or absolute.")
    parser.add_argument('--i_fa', help="Input sequences in FASTA format.",
                        required=True)
    parser.add_argument('--o_csv', help="Prediction results in csv format.",
                        required=True)
    parser.add_argument('--tmp_dir', required=True,
                        help=("Directory for temporary results"
                              "and computations. It is recommended to create "
                              "new directory for this "
                              "purpose as files in this directory may be "
                              "rewritten."))
    parser.add_argument('--model', default=Paths._MODEL,
                        help="Relative or absolute path to the model")
    parser.add_argument('--usearch', default=Paths._USEARCH,
                        help="Relative or absolute path to USEARCH executable")
    parser.add_argument('--espritz', default=Paths._ESPRITZ,
                        help="Relative or absolute path to Espritz executable")
    parser.add_argument('--fess', default=Paths._FESS,
                        help="Relative or absolute path to Fess executable")
    parser.add_argument('--tmhmm', default=Paths._TMHMM,
                        help="Relative or absolute path to TMHMM executable")
    parser.add_argument('--pdb', default=Paths._PDB_ECOLI_FA,
                        help="Relative or absolute path to PDB FASTA file.")

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()

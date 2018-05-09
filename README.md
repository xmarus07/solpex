# Solpex
Protein solubility predictor.

# Installation

  - Install Anaconda, see: https://conda.io/docs/user-guide/install/linux.html
  - Add anaconda bin directory to PATH if it  is not already there.
  - Clone repository: git clone https://github.com/xmarus07/Solpex.git
  - Create environment for Solpex:
	  1. cd Solpex
	  2. conda env create -f environment.yml

- Download and install additional tools:
	1. USEARCH: https://www.drive5.com/usearch/
	2. TMHMM: http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm
	3. Fess: http://protein.bio.unipd.it/download/
	4. Espritz: http://protein.bio.unipd.it/download/
- Change paths in class Paths located in file solpex.py or set them by arguments see --help defualt paths:
   
    _USEARCH = './additional_software/usearch/usearch10.0.240_i86linux32'  
	_PDB_ECOLI_FA = './data/Ecoli_xray_nmr_pdb_no_nesg.fa'  
    _TMHMM = './additional_software/tmhmm-2.0c.Linux/tmhmm-2.0c/bin/tmhmm'  
	_ESPRITZ = './additional_software/espritz/espritz.pl'  
	_FESS = './additional_software/fess/fess'

# Run

Before running Solpex, set envrionment to solpex: source activate solpex
	
	Exmaple for ./data/test.fa:
		python3 solpex.py --i_fa ./data/test.fa --o_csv test.csv --tmp_dir tmp
		
		--i_fa - Input FASTA file
		--o_csv - Output CSV file
		--tmp_dir - Directory for temporary results of additional tools
		
	Help:
		python3 solpex.py --help


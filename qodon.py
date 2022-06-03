'''

    mRNA Codon Optimization with Quantum Computers 
    Copyright (C) 2021  Dillion M. Fox, Ross C. Walker

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Example Usage:
> (qiskit) -bash$ python qodon.py 
IBM: 0.6752671392317706 0.6752671392317706
GA: 0.6752671392317706 0.6752671392317706 0.0
IBM: 0.7994909558630391 0.7370755249557397
GA: 0.7370755249557397 0.7370755249557397 0.0

'''

from src.classical_ga import CodonOptimization
from src.codon_bqm import DWaveBQM, QiskitBQM, ibm_score
from src.constants import *
from Bio import SeqIO
import numpy as np

# In order to use D-Wave or IBM, you must have access to appropriate
# libraries/devices.
use_dwave = False
use_ibm = True
use_ga = True
use_paper_seqs = False

if use_paper_seqs:
    # Path to fasta
    fasta_path = 'sample_seqs/covid_sequences.fasta'
    seq_list = [str(bio_seq.seq) for bio_seq in SeqIO.parse(fasta_path, 'fasta')]
    seq_name = [str(bio_seq.id) for bio_seq in SeqIO.parse(fasta_path, 'fasta')]
else:
    seq_list = ['AGM', 'GSK']
    seq_name = ['test_1', 'test_2']

for p_seq, seq_name in zip(seq_list, seq_name):

    if use_dwave:
        # Must have D-Wave libraries installed. If using hybrid or QPU methods,
        # must run code on Leap interface with appropriate resource access.
        # Quantum code is programmed to run n_execs times
        dwave_results = DWaveBQM(p_seq, exact=True, hybrid=False)
        print('D-Wave:',dwave_results.score, dwave_results.score_mean, dwave_results.score_std)

    if use_ibm:
        # Must have qiskit installed or run code on IBM Experience.
        # QiskitBQM will run 1 time.
        ibm_results = QiskitBQM(p_seq, exact=False, noise=False)
        print('IBM:',ibm_results.score, ibm_results.exact_score)

    if use_ga:
        # Run classical GA n_execs times
        c_scores = [CodonOptimization(p_seq).score for _ in range(n_execs)]
        print('GA:',min(c_scores),np.mean(c_scores), np.std(c_scores))


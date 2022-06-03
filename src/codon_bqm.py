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

'''
from src.penalties import OneBodyPenaltyMatrix, TwoBodyPenaltyMatrix
from collections import Counter, OrderedDict
from itertools import groupby, combinations
from src.constants import *
from src.scoring import SeqScorer
from typing import Optional
import numpy as np
import time, os


def valid_seq(p_seq,tb,i,raw):
    cond1 = list(raw.record[i][0]).count(1) == len(p_seq)
    cond2 = set([
        tb._codon_pos[_][0]
        for _ in np.where(raw.record[i][0] == 1)[0]
    ]) == set(range(len(p_seq)))
    if cond1 and cond2:
        return True
    return False


class DWaveBQM(object):
    def __init__(self,
                 p_seq: str,
                 hybrid: Optional[bool] = False,
                 exact: Optional[bool] = False):
        self.p_seq = p_seq
        self.offset = gc_constant * target_score**2
        self.exact = exact
        self.hybrid = hybrid
        self.execute()

    def execute(self):

        self._compute_penalties()
        self.score = 0
        self.n_seq = ''
        solutions = []

        for _ in range(n_execs):
            self._run_sampler()
            self._interpret_results()
            solutions.append(self.score)

        self.solutions = solutions
        self.score_mean = np.mean(self.solutions)
        self.score_std = np.std(self.solutions)
        self.score = min(solutions)

    def _compute_penalties(self):
        # Compute one-body penalties
        ob = OneBodyPenaltyMatrix(self.p_seq)
        self.h = ob.dwave_format
        self.ob = ob

        # Compute two-body penalties
        tb = TwoBodyPenaltyMatrix(self.p_seq)
        self.tb = tb
        self.J = tb.dwave_format

    def _interpret_results(self):

        # Find lowest energy nucleotide sequence
        self.result = min(self.filtered_samples, key=lambda x: x[1])
        self.n_seq = ''.join(
            list(np.array(self.tb._codon_code)[np.where(self.result[0] != 0)]))
        self.score = SeqScorer(self.n_seq).score

    def _run_sampler(self):

        # Import D-Wave libraries
        import dimod 

        # Construct custon binary quadratic model
        bqm = dimod.BinaryQuadraticModel(self.h,
                                         self.J,
                                         self.offset,
                                         dimod.BINARY,
                                         auto_scale=True)

        # Run the problem on the sampler and print the results
        if self.exact:
            sampler = dimod.ExactSolver()
            sampleset = sampler.sample(bqm)
        elif self.hybrid:
            from dwave.system import LeapHybridSampler
            sampler = LeapHybridSampler(minimum_time_limit=500.0)
            sampleset = sampler.sample(bqm)
        else:
            from dwave.system import EmbeddingComposite, DWaveSampler
            sampler = EmbeddingComposite(DWaveSampler())
            sampleset = sampler.sample(bqm,
                                       num_reads=5000,
                                       annealing_time=20,
                                       return_embedding=False)

        # Select samples with correct codon --> amino acid mapping
        self.filtered_samples = [[_[0], _[1]]
                                 for it, _ in enumerate(sampleset.record)
                                 if valid_seq(self.p_seq,self.tb,it,sampleset)]

        if len(self.filtered_samples) == 0:
            raise ValueError(
                '''No sequences with correct number of codons were found''')


class QiskitBQM(object):

    def __init__(self,
                 p_seq,
                 rseed = 1045,
                 noise = True,
                 exact = False):
        self.p_seq = p_seq
        self.offset = gc_constant * target_score**2
        self.use_noise = noise
        self.exact = exact
        self.score = -1
        self.exact_score = -1
        self.n_seq = ''
        self.exact_n_seq = ''
        self.rseed = rseed
        self.execute()

    def execute(self):
        self._compute_penalties()
        self._run_sampler()
        self._interpret_results()

    def _compute_penalties(self):
        # Compute one-body penalties
        ob = OneBodyPenaltyMatrix(self.p_seq)
        self.h = ob.dwave_format
        self.ob = ob

        # Compute two-body penalties
        tb = TwoBodyPenaltyMatrix(self.p_seq)
        self.tb = tb
        self.J = tb.dwave_format

    def _interpret_results(self):
        self.n_seq = ''.join(
            list(np.array(self.tb._codon_code)[np.where(self.result == 1.0)[0]]))
        self.best_score = SeqScorer(self.n_seq).score

        self.score = self.filter_score()

        self.exact_n_seq = ''.join(
            list(np.array(self.tb._codon_code)[np.where(self.exact_result == 1.0)[0]]))
        self.exact_score = SeqScorer(self.exact_n_seq).score

    def filter_score(self):
        record = list(self.raw.min_eigen_solver_result['eigenstate'].keys())
        filtered_samples = list(filter(lambda x: list(x).count('1') == len(self.p_seq), [[char for char in _] for _ in record]))
        valid_results = np.array([filtered_samples[i] for i in range(len(filtered_samples)) \
         if set([self.tb._codon_pos[_][0] for _ in np.where(np.array(filtered_samples[i])=='1')[0]]) == set(range(len(self.p_seq)))])
        n_seq_list = [''.join(_) for _ in [np.array(self.tb._codon_code)[np.where(valid_results[i] == '1')[0]] \
                                           for i in range(len(valid_results))]]
        scores = [SeqScorer(n).score for n in n_seq_list]
        if len(scores) > 0:
            best_score = min(scores)
            return best_score
        else:
            return np.nan
    
    def _run_sampler(self):

        # Import qiskit libraries
        from qiskit import BasicAer, Aer
        from qiskit.aqua import aqua_globals, QuantumInstance
        from qiskit.aqua.algorithms import QAOA, NumPyMinimumEigensolver
        from qiskit.optimization.algorithms import MinimumEigenOptimizer
        from qiskit.optimization import QuadraticProgram
        aqua_globals.massive=True
        os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'

        # Construct custon binary quadratic model
        mod = QuadraticProgram('codon_opt')
        for i in self.h.keys():
            mod.binary_var(name='q{}'.format(i))

        # Convert qubit labels to qiskit-friendly format
        qh = { 'q{}'.format(ind) : val for ind,val in self.h.items() }
        qJ = { ('q{}'.format(t[0]),'q{}'.format(t[1])) : val for t,val in self.J.items() }

        # Compose objective function in qiskit format
        mod.minimize(constant=self.offset, linear=qh, quadratic=qJ)
        
        # Convert problem to QAOA operator
        op, offset = mod.to_ising()

        if self.use_noise:
            from qiskit.providers.aer import QasmSimulator
            from qiskit.providers.aer.noise import NoiseModel
            from qiskit.test.mock import FakeVigo

            backend = Aer.get_backend('qasm_simulator')
            device = QasmSimulator.from_backend(FakeVigo())
            noise_model = NoiseModel.from_backend(device)
            #coupling_map = device.configuration().coupling_map

            quantum_instance = QuantumInstance(backend=backend, seed_simulator=self.rseed, seed_transpiler=self.rseed,
                                 noise_model=noise_model,)

        else:
            quantum_instance = QuantumInstance(BasicAer.get_backend('qasm_simulator'),
                                               seed_simulator=aqua_globals.random_seed,
                                               seed_transpiler=aqua_globals.random_seed,
                                               optimization_level=3)
        qaoa_mes = QAOA(quantum_instance=quantum_instance)
            
        qaoa = MinimumEigenOptimizer(qaoa_mes)   # using QAOA
        qaoa_result = qaoa.solve(mod)
        self.raw = qaoa_result
        self.result = np.array(list(qaoa_result))
        
        exact_mes = NumPyMinimumEigensolver()    # exact
        exact = MinimumEigenOptimizer(exact_mes)
        exact_result = exact.solve(mod)
        self.exact_result = np.array(list(exact_result))

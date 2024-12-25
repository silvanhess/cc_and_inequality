from pomegranate import *
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

CLC_data = np.load(file('BN_sample_data'))
CLC_data_valid = np.load(file('BN_sample_data_validation'))

overall_results = np.zeros([1,10])

VLAU_bins = pd.qcut(CLC_data[:, 2], 9)
VLAU = DiscreteDistribution.from_samples(VLAU_bins.codes)
VLAU_data = np.expand_dims(VLAU_bins.codes, axis=1)
ylHr0_whe_bins = pd.qcut(CLC_data[:, 4], 5)
ylHr0_whe = DiscreteDistribution.from_samples(ylHr0_whe_bins.codes)
ylHr0_whe_data = np.expand_dims(ylHr0_whe_bins.codes, axis=1)
gras200a_yld_bins = pd.qcut(CLC_data[:, 3], 10)
gras200a_yld = DiscreteDistribution.from_samples(gras200a_yld_bins.codes)
gras200a_yld_data = np.expand_dims(gras200a_yld_bins.codes, axis=1)

s_pop = Node(VLAU, name="Pop100m")
s_whe = Node(ylHr0_whe, name="ylHr0_whe")
s_add = Node(gras200a_yld, name="ExtraVar")

BN_data = np.concatenate((VLAU_data, ylHr0_whe_data, gras200a_yld_data , CLC_data[:, [0]]), axis=1)
CLC_before = ConditionalProbabilityTable.from_samples(BN_data, [VLAU, ylHr0_whe, gras200a_yld])
s_before = Node(CLC_before, name="CLC_before")

BN_data = np.concatenate((CLC_data[:, [0]], VLAU_data, ylHr0_whe_data, gras200a_yld_data , CLC_data[:, [1]]), axis=1)
cpt = ConditionalProbabilityTable.from_samples(BN_data, [CLC_before, VLAU, ylHr0_whe, gras200a_yld])
s_final = Node(cpt, name="CLC_after")

BN_data = np.concatenate((VLAU_data, ylHr0_whe_data, gras200a_yld_data , CLC_data[:, [1]]), axis=1)
CLC_after = ConditionalProbabilityTable.from_samples(BN_data, [VLAU, ylHr0_whe, gras200a_yld])
BN_data = np.concatenate((CLC_data[:, [1]], VLAU_data, ylHr0_whe_data, gras200a_yld_data , CLC_data[:, [0]]), axis=1)
cpt2 = ConditionalProbabilityTable.from_samples(BN_data, [CLC_after, VLAU, ylHr0_whe, gras200a_yld])

print(pd.Timestamp.now())

# validation for land-use
for k, L in enumerate([1,3,4,3,4]):
    if k<3:
        CLC_data_valid_u = CLC_data_valid[~(CLC_data_valid[:,0]==L),:]
        n_transitions = sum((CLC_data_valid[:,1]==L) & (~(CLC_data_valid[:,0]==L)))
    else:
        CLC_data_valid_u = CLC_data_valid[~(CLC_data_valid[:, 1] == L), :]
        n_transitions = sum((CLC_data_valid[:, 0] == L) & (~(CLC_data_valid[:, 1] == L)))
    data_points = len(CLC_data_valid_u)

    sims = range(0,data_points)
    results = np.zeros([data_points,2])
    if k < 3:
        results[:, 0] = CLC_data_valid_u[:, 1]
    else:
        results[:, 0] = CLC_data_valid_u[:, 0]
    for s in sims:
        s_before_v = CLC_data_valid_u[s, 0]
        s_pop_v = pd.cut(CLC_data_valid_u[s, [9]], VLAU_bins.categories).codes[0]
        s_whe_v = pd.cut(CLC_data_valid_u[s, [23]], ylHr0_whe_bins.categories).codes[0]
        s_add_v = pd.cut(CLC_data_valid_u[s, [13]], gras200a_yld_bins.categories).codes[0]
        s_final_v = CLC_data_valid_u[s, 1]
        if k < 3:
            class_probability = cpt.parameters[0][cpt.keymap[(s_before_v, s_pop_v, s_whe_v, s_add_v, L)]][5]
        else:
            class_probability = cpt2.parameters[0][cpt2.keymap[(s_final_v, s_pop_v, s_whe_v, s_add_v, L)]][5]
        results[s, 1] = class_probability

    results_rank = data_points - rankdata(results[:,1], method='max')
    ix = results_rank <= n_transitions
    results2 = results[ix,:]
    hits = sum(results2[:,0] == L)
    success_index = hits / sum(ix)
    random_result = n_transitions / data_points
    overall_results[0,k] = success_index
    overall_results[0,k+5] = random_result
    print([L,pd.Timestamp.now()])

np.savetxt(file('BN_validation'), overall_results, fmt='%10.4f', delimiter=',', )
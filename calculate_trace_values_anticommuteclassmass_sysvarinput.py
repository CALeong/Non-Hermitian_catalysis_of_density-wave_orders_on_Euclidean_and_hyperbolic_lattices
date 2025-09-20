import scipy.linalg
import numpy as np
from Fundamental.NonHermitian_Hamiltonian import NonHermitian_Honeycomb_periodic_boundary
from Fundamental.NonHermitian_Hamiltonian import NonHermitian_Hamiltonian
import sys

###
savedir = '/home/cal422/bir218_proj/cal422/NonHermitian_Matsubara_Summation_Calculations/Trace_Value_Computation_Results/HoneycombPBC_nl16_Data_AntiCommuting'
temperature = float(sys.argv[1])
num_levels = 16
alpha = float(sys.argv[2])
nindexbound = 3000
###

total_num_sites = np.size(NonHermitian_Honeycomb_periodic_boundary(num_levels, 1, 0), 0)

pauli_zero = np.array([[1, 0], [0, 1]])
pauli_one = np.array([[0, 1], [1, 0]])
pauli_three = np.array([[1, 0], [0, -1]])

def omegaval_from_index(nindex, tempval):
    return((2*nindex + 1)*np.pi*tempval)

tbham = np.kron(pauli_zero, NonHermitian_Honeycomb_periodic_boundary(num_levels, 1, 0))
massterm = np.diag(np.concatenate((np.repeat(1, int(total_num_sites / 2)), np.repeat(-1, int(total_num_sites / 2)))))
hmass = np.kron(pauli_one, massterm)

nonhermham = tbham + alpha*np.matmul(hmass, tbham)
massterm = np.kron(pauli_three, massterm)
del tbham
del hmass

final_trace_results = np.array([], dtype=np.complex128)
progress_counter = 0
for n in range(-nindexbound, nindexbound+1):
    if n % 100 == 0:
        print(n)

    omegaval = omegaval_from_index(n, temperature)

    nhterm = np.diag(np.repeat(1j*omegaval, 2*total_num_sites)) - nonhermham
    invterm = scipy.linalg.inv(nhterm)
    del nhterm

    massinvproduct = np.matmul(massterm, invterm)
    finalproduct = np.matmul(massinvproduct, massinvproduct)
    del massinvproduct

    traceval = np.trace(finalproduct)
    final_trace_results = np.append(final_trace_results, traceval)
    
    progress_counter = progress_counter + 1
    if progress_counter % 500 == 0:
        np.savetxt(savedir + '/honeycombPBC_nl{}a{}_temperature{}_upto{}_AntiCommuteMassClass_traceresults_INTERMEDIARY_BACKUP.txt'.format(num_levels, alpha, temperature, -nindexbound+progress_counter-1), final_trace_results)

np.savetxt(savedir + '/honeycombPBC_nl{}a{}_temperature{}_summationindexbound{}_AntiCommuteMassClass_traceresults.txt'.format(num_levels, alpha, temperature, nindexbound),
           final_trace_results)


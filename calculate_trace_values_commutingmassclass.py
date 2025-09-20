import scipy.linalg
import numpy as np
from Fundamental.NonHermitian_Hamiltonian import NonHermitian_Honeycomb_periodic_boundary
from Fundamental.NonHermitian_Hamiltonian import NonHermitian_Hamiltonian
import sys

###
savedir = '/home/cal422/bir218_proj/cal422/NonHermitian_Matsubara_Summation_Calculations/Trace_Value_Computation_Results/HoneycombPBC_nl20_Data'
temperature = float(sys.argv[1])
num_levels = 20
alpha = float(sys.argv[2])
nonhermham = NonHermitian_Honeycomb_periodic_boundary(num_levels, 1, alpha)
nindexbound = 4000
###

def omegaval_from_index(nindex, tempval):
    return((2*nindex + 1)*np.pi*tempval)

final_trace_results = np.zeros(2*nindexbound + 1, dtype=np.complex128)
for n in range(-nindexbound, nindexbound+1):
    if n%100 == 0:
        print(n)

    omegaval = omegaval_from_index(n, temperature)

    nhterm = np.diag(np.repeat(1j*omegaval, np.size(nonhermham, 0))) - nonhermham
    invterm = scipy.linalg.inv(nhterm)
    del nhterm

    massterm = np.diag(np.concatenate((np.repeat(1, int(np.size(nonhermham, 0) / 2)),
                                       np.repeat(-1, int(np.size(nonhermham, 0) / 2)))))

    massinvproduct = np.matmul(massterm, invterm)
    finalproduct = np.matmul(massinvproduct, massinvproduct)
    del massinvproduct

    traceval = np.trace(finalproduct)
    final_trace_results[n+nindexbound] = traceval

np.savetxt(savedir + '/honeycombPBC_nl{}a{}_temperature{}_summationindexbound{}_traceresults.txt'.format(num_levels, alpha, temperature, nindexbound),
           final_trace_results)


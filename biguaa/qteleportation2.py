

from qutip import *
import numpy as np
import math
import matplotlib.pyplot as plt
import qutip.qip
import scipy.stats
from qutip.qip.operations import snot, cnot, rx, ry, rz
from qutip.qobjevo import proj


# ############ FUNCTIONS FOR ANY REPRESENTATIONS OF QUANTUM STATES ################

def get_dict_form_qutip(wave_vector):
    dict_qubits = dict()

    n_qubits = str(int(math.log(len(wave_vector.full()), 2)))
    i = 0
    for x in wave_vector.full():
        if np.abs(x) > 0.00000000001:
            dict_qubits[format(i, '0'+n_qubits+'b')] = x
        i += 1
    return dict_qubits


def get_qutip_form_dict(dict_qubits):
    list_qubits = list(dict_qubits.items())
    n_qubits = len(list_qubits[0][0])
    dim_space = pow(2, n_qubits)
    wave_vector = np.zeros(dim_space, dtype=complex)
    for x in list_qubits:
        wave_vector[int(x[0], 2)] += x[1]
    return Qobj(wave_vector)


def print_dict(dictionary, precision=3):
    np.set_printoptions(precision=3)
    list_data = list(dictionary.items())
    pr_ch = str(precision)
    for data in list_data:
        print(('({0.real:.'+pr_ch+'f} + {0.imag:.'+pr_ch+'f}j)').format(data[1][0])+"|{}>".format(data[0]))

# ########################### Main structure ####################

# Initialize a random 3-qubit quantum state
def initialize(state):
    bell = bell_state("00")
    return tensor(state, bell).unit()

# Make transformations until the projective measurements
def evolute(state):
    state = cnot(3, 0, 1)*state
    state = snot(3, 0)*state

    return state.unit()

# Make projective measurements on 0, 1 qubits
def measure(state):
    # Trace over 2 qubit
    rho_trace = ket2dm(state).ptrace((0, 1))
    # Get amplitudes of Bell's state
    pr_set = rho_trace.diag()
    # Random projective measurement
    mres = np.random.choice(range(0, 4), 1, False, pr_set)[0]
    return mres


# Make teleportation of an input state
def teleport(state, mres):
    X = rx(-math.pi, 3, 2)
    Z = rz(-math.pi, 3, 2)
    print(mres)

    if mres == 0:
        bra0 = bra("000")
        bra1 = bra("001")
        state =  state
    if mres == 1:
        bra0 = bra("010")
        bra1 = bra("011")
        state = Z*state
    if mres == 2:
        bra0 = bra("100")
        bra1 = bra("101")
        state = X*state
    if mres == 3:
        bra0 = bra("110")
        bra1 = bra("111")
        state = Z*X*state

    pr0 = (bra0*state).tr()
    pr1 = (bra1*state).tr()

    return (pr0*basis(2, 0) + pr1*basis(2, 1)).unit()

if __name__ == '__main__':
    state = rand_ket(N=2, density=1, dims=[[2], [1]], seed=None)
    state0 = state.copy()
    print("Begin state of 1th qubit:")
    print_dict(get_dict_form_qutip(state))
    state = initialize(state=state)
    state = evolute(state)
    mres = measure(state)
    state = teleport(state, mres)
    print("\nState of 3th qubit after teleportation: ")
    print_dict(get_dict_form_qutip(state))
    print("Fidelity between the input and the output: ", fidelity(state, state0))


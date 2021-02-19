#!/usr/bin/env python
# coding: utf-8

# In[10]:


from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from qutip.qip.circuit import QubitCircuit, Gate
from qutip.qip.operations import gate_sequence_product
import math


def initialize(state):
    zero = basis(2, 0)
    one = basis(2, 1)
    phi_plus = (tensor(one, one) + tensor(zero, zero)).unit()
    return tensor(state, phi_plus)


def evolute(state):
    qc = QubitCircuit(N=3)
    qc.add_gate("CNOT", controls=0, targets=1)
    qc.add_gate("SNOT", targets=0)
    U_list = qc.propagators()
    return gate_sequence_product(U_list) * state


def measure(state):
    pho12 = state.ptrace([0,1])
    basisstate = [basis(2,0), basis(2,1)]
    proj = [[tensor(basisstate[j],basisstate[i]) 
                      * tensor(basisstate[j],basisstate[i]).dag() for i in range(2)] for j in range(2)]
    prob = np.array([expect(proj[(i >> 1)][i & 1], pho12) for i in range(4)])
    return np.random.choice([0,1,2,3], p=prob)


def teleport(state, mres):
    qc = QubitCircuit(N=1)
    if mres ==0:
        state = state[0][0][0] * basis(2,0) + state[1][0][0] * basis(2,1)
    elif mres == 1:
        state = state[2][0][0] * basis(2,0) + state[3][0][0] * basis(2,1)
        state = -1j * state
        qc.add_gate("RX", targets=0, arg_value=-math.pi)
    elif mres == 2:
        state = state[4][0][0] * basis(2,0) + state[5][0][0] * basis(2,1)
        state = -1j * state
        qc.add_gate("RZ", targets=0, arg_value=-math.pi)
    else:
        state = state[6][0][0] * basis(2,0) + state[7][0][0] * basis(2,1)
        state = -1 * state
        qc.add_gate("RX", targets=0, arg_value=-math.pi)
        qc.add_gate("RZ", targets=0, arg_value=-math.pi)
    U_list = qc.propagators()
    state = gate_sequence_product(U_list) * state
    return state.unit()
if __name__ == '__main__':
    st = evolute(initialize(rand_ket(N=2)))
    print(teleport(st, measure(st)))


# In[ ]:





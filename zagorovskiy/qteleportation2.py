from qutip import *
from numpy import *
import random

# телепортируемое состояние
#q1=rand_ket(2,1)
#print(q1)
H=Qobj(array([[1,1],[1,-1]])/2**0.5)
CX=Qobj(array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]),[[2, 2], [2, 2]])
I=Qobj(eye(2))
# вид Белловского состояния
s2=1; s3=1

def initialize(state):
    q23=ket([s2,s3],2)
    q23=CX*tensor(H,I)*q23
    phi=tensor(state,q23)
    return phi
def evolute(state):
    phi=tensor(tensor(H,I)*CX,I)*state
    return phi
def measure(state):
    p=0; r=random.random()
    for i in range(4):
        v=ket([i>>1,i%2],2)
        P=v*v.trans()
        Ps=tensor(P,I)
        phi=Ps*state
        p=p+real(trace(phi.trans().conj()*phi))
        if(p>r):
            break
    return i
def teleport(state, mres):
    v = ket([mres>>1,mres%2], 2)
    P = v * v.trans()
    Ps = tensor(P, I)
    phi = Ps * state
    phi=phi/(real(trace(phi.trans().conj()*phi)))**0.5
    psi=Qobj(array([phi[mres*2][0],phi[mres*2+1][0]]))
    psi=sigmax()**(mod(s3+mres%2,2))*psi
    psi = sigmaz() ** (mod(s2 + (mres >> 1), 2)) * psi
    return psi

# схема телепортации
"""phi1=initialize(q1)
phi2=evolute(phi1)
m=measure(phi2)
#print(m)
psi=teleport(phi2,m)
#print(psi)"""

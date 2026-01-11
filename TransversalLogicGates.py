from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info import Pauli
from qiskit.circuit.library import UnitaryGate
import numpy as np

#Add Logical gates X_L and Z_L and a transversal implementation of Logical CNOT between two logical qubits of the 5 qubit code.


#Logical X and Z Operators

X_L = Pauli("XXXXX").to_instruction()
Z_L = Pauli("ZZZZZ").to_instruction()
Y_L = Pauli("YYYYY").to_instruction()

#Transversal K_{s_x,s_y,s_z} operators.

def K_gate(sx,sy,sz):
    label_string = ''
    for value in [sx,sy,sz]:
        if value == 1:
            label_string.append('+')
        elif value == -1:
            label_string.append('-')
    
    X = np.array([[0.0,1.0],[1.0,0.0]])
    Y = np.array([[0.0,-1.0j],[1.0j,0.0]])
    Z = np.array([[1.0,0.0],[0.0,-1.0]])

    K = np.exp((1j*np.pi)/(3*np.sqrt(3))(sx*X + sy*Y + sz*Z))
    K = UnitaryGate(K, label = 'K_('+label_string+')')

    return K

#Y_L = Pauli("iXXXXX").dot(Pauli("ZZZZZ")).to_instruction()
#Y_L_ = np.sqrt(Pauli("iXXXXX").dot(Pauli("ZZZZZ")).to_matrix()) #Y_L = i*(X_L)(Z_L)


#Transversal Logical CX and CZ Operators

def transversal_CX():
    qc = QuantumCircuit(10)
    for i in range(int(qc.num_qubits/2)):
        qc.cx(i,i+int(qc.num_qubits/2))

    return qc

def transversal_CZ():
    qc = QuantumCircuit(10)
    for i in range(int(qc.num_qubits/2)):
        qc.cz(i,i+int(qc.num_qubits/2))

    return qc

def transversal_CY():
    qc = QuantumCircuit(10)
    for i in range(int(qc.num_qubits/2)):
        qc.cy(i,i+int(qc.num_qubits/2))

    return qc

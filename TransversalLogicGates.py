from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info import Pauli

#Add Logical gates X_L and Z_L and a transversal implementation of Logical CNOT between two logical qubits of the 5 qubit code.


#Logical X and Z Operators

X_L = Pauli("XXXXX").to_instruction()
Z_L = Pauli("ZZZZZ").to_instruction()

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
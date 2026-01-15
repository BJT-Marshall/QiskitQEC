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
    """
    Constructs a unitary gate QisKit instruction :math:`K_{s_x,s_y,s_z} = \exp{\\frac{i\pi(s_xX+s_yY+s_zZ)}{3\sqrt{3}}}`. All 8 possible :math:`K` gates are transversal in the 5-qubit-code
    and can therefore be implemented across a logical qubit as :math:`K := K_{s_x,s_y,s_z}^{\otimes 5}`.
    
    :param sx: Coefficient for the Pauli :math:`X` gate in the construction of the :math:`K_{s_x,s_y,s_z}` gate. Can take the value :math:`+1` or :math:`-1`.
    :type sx: int
    :param sy: Coefficient for the Pauli :math:`Y` gate in the construction of the :math:`K_{s_x,s_y,s_z}` gate. Can take the value :math:`+1` or :math:`-1`.
    :type sy: int
    :param sz: Coefficient for the Pauli :math:`Z` gate in the construction of the :math:`K_{s_x,s_y,s_z}` gate. Can take the value :math:`+1` or :math:`-1`.
    :type sz: int

    :returns K: The single qubit unitary gate :math:`K_{s_x,s_y,s_z}`.
    :rtype: UnitaryGate 
    """
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
    """
    Generates a transversal implementation of the CNOT operation in the 5-qubit-code to be applied to a pair of encoded logical qubits.

    :returns: Transversal implementation of the CNOT operation in the 5-qubit-code.
    :rtype: QuantumCircuit
    """
    qc = QuantumCircuit(10)
    for i in range(int(qc.num_qubits/2)):
        qc.cx(i,i+int(qc.num_qubits/2))

    return qc

def transversal_CZ():
    """
    Generates a transversal implementation of the CZ operation in the 5-qubit-code to be applied to a pair of encoded logical qubits.

    :returns: Transversal implementation of the CZ operation in the 5-qubit-code.
    :rtype: QuantumCircuit
    """
    qc = QuantumCircuit(10)
    for i in range(int(qc.num_qubits/2)):
        qc.cz(i,i+int(qc.num_qubits/2))

    return qc

def transversal_CY():
    """
    Generates a transversal implementation of the CY operation in the 5-qubit-code to be applied to a pair of encoded logical qubits.

    :returns: Transversal implementation of the CY operation in the 5-qubit-code.
    :rtype: QuantumCircuit
    """
    qc = QuantumCircuit(10)
    for i in range(int(qc.num_qubits/2)):
        qc.cy(i,i+int(qc.num_qubits/2))

    return qc

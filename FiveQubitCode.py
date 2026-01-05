from qiskit import QuantumCircuit
from qiskit.quantum_info import Pauli

#Syndrome Operators of the five-qubit code:

M1 = Pauli("XZZXI").to_instruction()
M2 = Pauli("IXZZX").to_instruction()
M3 = Pauli("XIXZZ").to_instruction()
M4 = Pauli("ZXIXZ").to_instruction()

#Logical Qubits 
# 
#|0>_L = 1/4 (I + M1)(I + M2)(I + M3)(I+ M4)|00000>
#|1>_L = 1/4 (I + M1)(I + M2)(I + M3)(I+ M4)|11111>
#State vector given in computational basis in increasing binary order.

logi_0_vector = [0.25,0,0,0.25,0,-0.25,0.25,0,0,-0.25,-0.25,0,0.25,0,0,-0.25,0,0.25,-0.25,0,-0.25,0,0,-0.25,0.25,0,0,-0.25,0,-0.25,-0.25,0]
logi_1_vector = [0,-0.25,-0.25,0,-0.25,0,0,0.25,-0.25,0,0,-0.25,0,-0.25,0.25,0,-0.25,0,0,0.25,0,-0.25,-0.25,0,0,0.25,-0.25,0,0.25,0,0,0.25]

def logi_0():
    qc = QuantumCircuit(5)
    qc.initialize(logi_0_vector,range(5))

    return qc

def logi_1():
    qc = QuantumCircuit(5)
    qc.initialize(logi_1_vector,range(5))

    return qc

#Logical X and Z Operators

X_L = Pauli("XXXXX").to_instruction()
Z_L = Pauli("ZZZZZ").to_instruction()

#Syndrome Detection on a single logical qubit

def syndrome_detection(logical_qubit):
    syndromes = [M1,M2,M3,M4]
    qc = QuantumCircuit(logical_qubit.num_qubits+4,4)
    qc.append(logical_qubit,range(4,4+logical_qubit.num_qubits))
    for i in range(qc.num_qubits - logical_qubit.num_qubits):
        qc.h(i)
        #controled syndrome measurements
        CS = syndromes[i].control(1)
        qubit_inds = [j + qc.num_qubits - logical_qubit.num_qubits for j in range(logical_qubit.num_qubits)]
        qubit_inds.insert(0,i)
        qc.append(CS,qubit_inds) #[qubit_inds.insert(0,i)]
        qc.h(i)
    
    qc.measure(range(qc.num_qubits - logical_qubit.num_qubits),range(qc.num_qubits - logical_qubit.num_qubits))

    return qc

#Test:

qc = syndrome_detection(logi_0())
qc.draw(output = 'mpl', filename = 'ErrorFreeSyndromeCirc')

qc1 = QuantumCircuit(9,4)
qc1.x(4) #X1 Error
qc1.append(syndrome_detection(logi_0()),range(9),range(4))
qc1.draw(output = 'mpl', filename = 'X1ErrorSyndromeCirc')


#Error lookup table. Each single qubit error will result in a unique key measured from the ancillary qubits in the computational basis, represented in binary.

error_binaries = [format(i,'b') for i in range(16)]
error_paulis = [None,"X1","Z3","X5","Z5","Z2","X4","Y5","X2","Z4","Z1","Y1","X3","Y2","Y3","Y4"]






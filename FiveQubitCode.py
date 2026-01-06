from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
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

logi_0_vector = [0.25,0,0,-0.25,0,0.25,-0.25,0,0,0.25,0.25,0,-0.25,0,0,-0.25,0,-0.25,0.25,0,0.25,0,0,-0.25,-0.25,0,0,-0.25,0,-0.25,-0.25,0]
logi_1_vector = [0,-0.25,-0.25,0,-0.25,0,0,-0.25,-0.25,0,0,0.25,0,0.25,-0.25,0,-0.25,0,0,-0.25,0,0.25,0.25,0,0,-0.25,0.25,0,-0.25,0,0,0.25]

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
    clr = ClassicalRegister(4)
    qr = QuantumRegister(logical_qubit.num_qubits+clr.size)
    qc = QuantumCircuit(qr,clr)
    qc.append(logical_qubit,range(4,4+logical_qubit.num_qubits))
    for i in range(clr.size):
        qc.h(i)
        #controled syndrome measurements
        CS = syndromes[i].control(1)
        qubit_inds = [j + clr.size for j in range(logical_qubit.num_qubits)]
        qubit_inds.insert(0,i)
        qc.append(CS,qubit_inds) #[qubit_inds.insert(0,i)]
        qc.h(i)
    
    qc.measure(range(clr.size),range(clr.size))

    return qc, clr

#Test:

"""qc,clr = syndrome_detection(logi_0())
qc.draw(output = 'mpl', filename = 'ErrorFreeSyndromeCirc')

qc1 = QuantumCircuit(9,4)
qc1.x(4) #X1 Error
temp,clr_temp = syndrome_detection(logi_0())
qc1.append(temp,range(9),range(4))
qc1.draw(output = 'mpl', filename = 'X1ErrorSyndromeCirc')"""


#Error lookup table. Each single qubit error will result in a unique key measured from the ancillary qubits in the computational basis, represented in binary.

#error_binaries = [format(i,'b') for i in range(16)]
error_binaries = [0b0000,0b0001,0b0010,0b0011,0b0100,0b0101,0b0110,0b0111,0b1000,0b1001,0b1010,0b1011,0b1100,0b1101,0b1110,0b1111]
error_paulis = [None,"X1","Z3","X5","Z5","Z2","X4","Y5","X2","Z4","Z1","Y1","X3","Y2","Y3","Y4"]


#Need to think of a clever way to integrate the dynamic qc.if_test statements that allow 4 if tests to correctly identify the single qubit error.

#step 1) measure ALL 
#step 2) run error search
#step 3) convert error code to circuit instruction
#step 4) apply error
# 

def single_qb_EC(qc,e_code):
    """Step 3,4"""
    ind = error_binaries.index(e_code)
    EC_str = error_paulis[ind]
    offset = qc.num_clbits-1
    if EC_str == None:
        print("No Error")
    elif EC_str[0] == 'X':
        qc.x(offset+int(EC_str[1]))
        #print(EC_str + " error corrected.")
    elif EC_str[0] == 'Y':
        qc.y(offset+int(EC_str[1]))
        #print(EC_str + " error corrected.")
    elif EC_str[0] == 'Z':
        qc.z(offset+int(EC_str[1]))
        #print(EC_str + " error corrected.")

    return qc


def dynamic_error_correction(qc,clr):
    """syndrome detection qc inputted, steps 1,2"""
    for e_code in error_binaries:
        with qc.if_test((clr,e_code)):
            qc = single_qb_EC(qc,e_code)
    
    qc.draw(output = 'mpl', filename = '0')

    return qc



def QEC(logical_qubit):
    """Compile the above functions to perform single qubit QEC on a logical qubit in the 5 qubit code"""
    qc, clr = syndrome_detection(logical_qubit)
    qc = dynamic_error_correction(qc,clr)

    return qc


#Want to run the QEC returned circuit

zero_L = logi_0()
#manual error
#zero_L.x(0)
#
qc = QEC(zero_L)

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram, plot_state_city

simulator = AerSimulator()
circ = transpile(qc, simulator)

# Run and get counts
result = simulator.run(circ, shots = 10000).result()

counts = result.get_counts(circ)
plot_histogram(counts, title='Final Logical Qubit Counts', filename = '1')









from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info import Pauli

#Syndrome Operators of the five-qubit code:
#Pauli strings need to be 'reversed' for correct implementation as the Pauli class takes the string as P_(n-1) \otimes ... NOT P_0 \otimes ...

M1 = Pauli("IXZZX").to_instruction() #XZZXI syndrome operator
M2 = Pauli("XZZXI").to_instruction() #IXZZX syndrome operator
M3 = Pauli("ZZXIX").to_instruction() #XIXZZ syndrome operator
M4 = Pauli("ZXIXZ").to_instruction() #ZXIXZ syndrome operator

#Logical Qubits 
# 
#|0>_L = 1/4 (I + M1)(I + M2)(I + M3)(I+ M4)|00000>
#|1>_L = 1/4 (I + M1)(I + M2)(I + M3)(I+ M4)|11111>
#State vector given in computational basis in increasing binary order.

logi_0_vector = [0.25,0,0,-0.25,0,0.25,-0.25,0,0,0.25,0.25,0,-0.25,0,0,-0.25,0,-0.25,0.25,0,0.25,0,0,-0.25,-0.25,0,0,-0.25,0,-0.25,-0.25,0]
logi_1_vector = [0,-0.25,-0.25,0,-0.25,0,0,-0.25,-0.25,0,0,0.25,0,0.25,-0.25,0,-0.25,0,0,-0.25,0,0.25,0.25,0,0,-0.25,0.25,0,-0.25,0,0,0.25]

def logi_0():
    """
    Returns a 5 qubit QuantumCircuit object initialised in the 5-qubit-code logical zero state, |0>_L. 
    """
    qc = QuantumCircuit(5)
    qc.initialize(logi_0_vector,range(5))

    return qc

def logi_1():
    """
    Returns a 5 qubit QuantumCircuit object initialised in the 5-qubit-code logical one state, |1>_L. 
    """
    qc = QuantumCircuit(5)
    qc.initialize(logi_1_vector,range(5))

    return qc

#Logical X and Z Operators

X_L = Pauli("XXXXX").to_instruction()
Z_L = Pauli("ZZZZZ").to_instruction()

#Syndrome Detection on a single logical qubit

def syndrome_detection(logical_qubit):
    """
    Takes in a single logical qubit as a QuantumCircuit object, applies all syndrome operators controlled via seperate ancillary qubits and measures
    the resulting syndrome measurments to a classical register in the order (Syndrome Result, cr) = ((M1,3),(M2,2),(M3,1),(M4,0)). This reverse ordering
    is done to produce the correct binary keys for error identification.
    
    :param logical_qubit: 5 qubit QuantumCircuit object containing the current state of a single logical qubit.
    :return qc: The complete syndrome detection circuit including the intial logical qubit state and syndrome measurments.
    :return clr: The classical register to which the syndrome measurement results are stored.
    """
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
    
    qc.measure(range(clr.size),range(clr.size-1,-1,-1))

    return qc, clr

#Error lookup table. Each single qubit error will result in a unique key measured from the ancillary qubits in the computational basis, represented in binary.

#error_binaries = [format(i,'b') for i in range(16)] #Needed the leading zeros to compare to counts so manually inputted binaries here.
error_binaries = [0b0000,0b0001,0b0010,0b0011,0b0100,0b0101,0b0110,0b0111,0b1000,0b1001,0b1010,0b1011,0b1100,0b1101,0b1110,0b1111]
error_paulis = [None,"X1","Z3","X5","Z5","Z2","X4","Y5","X2","Z4","Z1","Y1","X3","Y2","Y3","Y4"]



#Notes:
#step 1) measure all ancillary qubits ##dynamic_error_correction
#step 2) run error search through dynamic circuit
#step 3) convert error code to circuit instruction ##single_qb_EC
#step 4) apply error correction
# 

def single_qb_EC(qc,e_code):
    """
    Given an error code string, for example 'X1' for a bit-flip error on the first qubit, applies the correct error correcting gate and returns the QuantumCircuit object.
    This method is used in conjunction with Qiskit's 'QuantumCircuit.if_test()' method to make this conditional on syndrome measurement through ancillary qubits.
    
    :param qc: QuantumCircuit object to which the single qubit error correcting gate should be applied.
    :param e_code: The error code string defining the single qubit error to be corrected. e.g. 'X1' is a bit-flip error on the first qubit.
    :return qc: The error corrected QuantumCircuit object.
    """
    ind = error_binaries.index(e_code)
    EC_str = error_paulis[ind]
    offset = qc.num_clbits-1
    if EC_str == None:
        None
    elif EC_str[0] == 'X':
        qc.x(offset+int(EC_str[1]))
    elif EC_str[0] == 'Y':
        qc.y(offset+int(EC_str[1]))
    elif EC_str[0] == 'Z':
        qc.z(offset+int(EC_str[1]))

    return qc


def dynamic_error_correction(qc,clr):
    """
    Taking in a circuit produced by the 'syndrome_detection' method and the ClassicalRegister that will contain the results of the syndrome measurment,
    Qiskits dynamic programming capabilities are used to test for all possible single qubit errors and apply the corresponding error correcting gate.
    
    :param qc: Post syndrome measurement QuantumCircuit object within which single qubit errors may need to be corrected.
    :param clr: The classical register holding the results of previously applied syndrome measurments.
    :return qc: The complete single qubit error correcting circuit from which the logical qubit statevcan be re-extracted.
    """
    for e_code in error_binaries:
        with qc.if_test((clr,e_code)):
            qc = single_qb_EC(qc,e_code)

    return qc



def QEC(logical_qubit):
    """
    Compilation of the above 'syndrome_measurment' and 'dynamic_error_correction' functions to perform single qubit error correction on a single 
    logical qubit of the 5 qubit code.
    
    :param logical_qubit: 5 qubit QuantumCircuit object containing the current state of a single logical qubit.
    :return qc: The complete single qubit error correcting circuit for 'logical_qubit' from which the logical qubit state can be re-extracted.
    """
    
    qc, clr = syndrome_detection(logical_qubit)
    qc = dynamic_error_correction(qc,clr)

    return qc



#Examples:

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram

#No Error Example:------------------------------------------------------------------------------------------------------


zero_L = logi_0() #Initialise a logical qubit in the state |0>_L
qc = QEC(zero_L) #Create the single qubit QEC circuit for this logical qubit
qc.draw(output = 'mpl', filename = 'NoErrorExampleCircuit') #Complete error correcting circuit visualization.

simulator = AerSimulator()
circ = transpile(qc, simulator)

result = simulator.run(circ).result() #Run and get counts
counts = result.get_counts(circ)
plot_histogram(counts, title='Syndrome Measurement Result', filename = 'NoErrorExampleCounts') 
#All counts will be in the state |0000> of the ancillary qubits, corresponding to no error.


#Single Qubit Error Example:--------------------------------------------------------------------------------------------

zero_L = logi_0() #Initialise a logical qubit in the state |0>_L
#-----------------------ERROR----------------------------
zero_L.x(3) #Manually add a bit-flip error on the 4'th qubit. 
#--------------------------------------------------------
qc = QEC(zero_L) #Create the single qubit QEC circuit for this logical qubit
qc.draw(output = 'mpl', filename = 'X4ErrorExampleCircuit') #Complete error correcting circuit visualization.

simulator = AerSimulator()
circ = transpile(qc, simulator)

result = simulator.run(circ).result() #Run and get counts
counts = result.get_counts(circ)
plot_histogram(counts, title='Syndrome Measurement Result', filename = 'X4ErrorExampleCounts') 
#All counts will be in the state |0110> of the ancillary qubits, corresponding to no a bit flip error on the 4'th qubit. 
#This error is then corrected by application of a further bit-flip on the 4'th qubit (i.e. X4).








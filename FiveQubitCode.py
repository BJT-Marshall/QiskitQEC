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
    Creates a 5 qubit QuantumCircuit object initialised in the 5-qubit-code logical zero state, :math:`\ket{0}_L`. 
    Used for initialising circuits in the 5-qubit code, preparing them for logical operations and syndrome measurements,
    Appending the returned QuantumCircuit to an existing instance of a logical qubit resets the qubit to the :math:`\ket{0}_L` state.

    :returns: Initialised :math:`\ket{0}_L` state.
    :rtype: QuantumCircuit

    """
    qc = QuantumCircuit(5)
    qc.initialize(logi_0_vector,range(5))

    return qc

def logi_1():
    """
    Creates a 5 qubit QuantumCircuit object initialised in the 5-qubit-code logical zero state, :math:`\ket{1}_L`. 
    Used for initialising circuits in the 5-qubit code, preparing them for logical operations and syndrome measurements,
    Appending the returned QuantumCircuit to an existing instance of a logical qubit resets the qubit to the :math:`\ket{1}_L` state.

    :returns: Initialised :math:`\ket{1}_L` state.
    :rtype: QuantumCircuit

    """
    qc = QuantumCircuit(5)
    qc.initialize(logi_1_vector,range(5))

    return qc


#Syndrome Detection on a single logical qubit

def syndrome_detection(logical_qubit):
    """
    Takes a single logical qubit encoded in the 5-qubit-code possibly containing a single qubit error and creates a syndrome detection circuit to identify the error. 
    Resulting syndrome measurements are stored in a ClassicalRegister in the order :math:`(result, clbit) = \{(M1,3),(M2,2),(M3,1),(M4,0)\}`. 
    This reverse ordering is done to produce the correct binary keys for error identification.
    
    
    :param logical_qubit: Logical qubit containing the current possibly error prone state.
    :type logical_qubit: QuantumCircuit
    :returns qc: The complete syndrome detection circuit including syndrome detection via ancillary qubits and subsequent syndrome measurment.
    :returns clr: The classical register in which the syndrome measurement results are stored.
    :rtype qc: QuantumCircuit
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
    Applies the required single qubit error correcting gate given an inputted quantum circuit and a single qubit error code.

    Example:
    ```python
    from qiskit import QuantumCircuit

    logical_qubit = logi_0() #error free logical qubit state.
    logical_qubit.x(0) #manually introduce bit flip error on the 1st physical qubit.

    error_corrected_logical_qubit = single_qb_EC(logical_qubit,'X1') #correct error through application of an additional bit flip gate on the 1st qubit.
    ```
    
    :param qc: The quantum circuit to which the single qubit error correcting gate should be applied.
    :type qc: QuantumCircuit
    :param e_code: The error code string defining the single qubit error to be corrected.
    :type e_code: str
    :returns: The quantum circuit corrected for the inputted single qubit error.
    :rtype: QuantumCircuit


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
    Takes in a circuit produced by the 'syndrome_detection' method and the ClassicalRegister that will contain the results of the syndrome measurment. Subsequently tests the contents
    of the classical register against a map of known single qubit error syndrome measurements to determine the single qubit error and applies the corresponding gate.
    
    :param qc: Quantum circuit produced by 'syndrome_detection' within which a single qubit error may need to be corrected.
    :type qc: QuantumCircuit
    :param clr: The classical register holding the results of previously applied syndrome measurments.
    :type clr: ClassicalRegister
    :returns qc: The complete single qubit error correcting circuit from which the logical qubit state can be further manipulated or extracted.
    :rtype: QuantumCircuit
    """
    for e_code in error_binaries:
        with qc.if_test((clr,e_code)):
            qc = single_qb_EC(qc,e_code)

    return qc



def QEC(logical_qubit):
    """
    Given a logical qubit encoded in the 5-qubit-code possibly including a single qubit error, single qubit error correction is performed through syndrome measurements using ancillary qubits 
    and the single qubit error corrected. 
    
    :param logical_qubit: Quantum circuit containing the current state of a single logical qubit, possibly containing a single qubit error.
    :type logical_qubit: QuantumCircuit
    :returns qc: The complete single qubit error correcting circuit for 'logical_qubit' from which the logical qubit state can be further manipulated or extracted.
    :returns clr: The classical register to which syndrome measurements are saved. The contents of this classical register are reset and is only returned for mapping 'qc' to larger quantum circuits.
    :rtype qc: QuantumCircuit
    """
    
    qc, clr = syndrome_detection(logical_qubit)
    qc = dynamic_error_correction(qc,clr)
    qc.reset(range(clr.size)) #Reset the Ancillary Qubits for future QEC cycles.

    return qc, clr



#Examples:

#Bool to change to run examples:
example = False

if example:

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








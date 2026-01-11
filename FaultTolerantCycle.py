from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

#Local Imports
from FiveQubitCode import QEC, logi_0, logi_1
from TransversalLogicGates import X_L, Z_L, Y_L, transversal_CX, transversal_CZ, transversal_CY, K_gate

#Goal is to make some sort of framework that puts things together
#At the moment EACH logical qubit uses a seperate set of ancillary qubits. This is good if qubit coherence time is a limiting 
#factor as resuing the same ancillary qubits would increase computation time between each logical operation, however it does nearly double the amount
#of qubits needed in a circuit to 9n for n logical qubits as opposed to 5n+4 for n logical qubits all using the same ancillas.

#If all logical qubits used the same ancillas then the qubit distance would be another problem. I think this is one of the motivation for
#topological alternative codes like Surface Codes.

#Using one set of ancillas would allow for the easier construction of circuits though... Actually it wouldnt you can just do all indexing mod 9 
#instead of mod 5. But this increases qubit distance for transversal operations between logocal qubits.

#Lets to both and look at the circuits, then can come up with reasoning for a decision.

def cycle(n_logical_qubits, operation_list, structure = '9n'):

    #Initializes logical qubits in the |0>_L state.

    #operation_list = [['CX',0,1],['Z',2],...]

    #structure = '9n' for structure of n(4ancilla5calcqubit) blocks
    #structure = '4n5n' for structure of (4n)ancilla(5n)calcqubit blocks #NOT SUPPORTED

    #Qubit indexing for structure of n(4ancilla5calcqubit) (Logical Qubit Indexing to start at zero like with physical qubits)
    #logical_qubit_n will be on the last 5 physical qubits of block [9*n,9(n+1)-1] i.e. the block [9n+4,9n+8]
    #Ancillas for logical_qubit_n will be [9n,9n+3]
    qubit_map = []
    if structure == '9n':
        logical_qubit_map = [[x for x in range(9*n+4,9*n+9)] for n in range(n_logical_qubits)]
        ancillary_qubit_map = [[x for x in range(9*n,9*n+4)] for n in range(n_logical_qubits)]
        qubit_map = [[x for x in range(9*n,9*n+9)] for n in range(n_logical_qubits)]
    elif structure == '4n5n':
        logical_qubit_map = [[x for x in range(4*n_logical_qubits,4*n_logical_qubits+5*n)] for n in range(n_logical_qubits)]
        ancillary_qubit_map = [[x for x in range(4*n,4*n+4)] for n in range(n_logical_qubits)]
        for i in range(n_logical_qubits):
            temp = []
            for element in ancillary_qubit_map[i]:
                temp.append(element)
            for element in logical_qubit_map[i]:
                temp.append(element)
            qubit_map.append(temp)
    else:
        print("Unsupported Structure")
    
    #Create the QuantumCircuit object to add the logical qubits to.
    qc = QuantumCircuit()
    for i in range(n_logical_qubits):
        qc.add_register(QuantumRegister(9, 'q'+str(i)))


    initial_logical_qubits = [logi_0() for i in range(n_logical_qubits)]
    #Initial error correction for single qubit errors in state preperation
    #initial_QEC = [QEC(Lq) for Lq in initial_logical_qubits]
    clrs = []
    initial_QEC = []
    for Lq in initial_logical_qubits:
        QC, clr = QEC(Lq)
        initial_QEC.append(QC)
        qc.add_register(clr)
        clrs.append(clr)
    
    for i in range(len(initial_QEC)):
        initial_QEC[i].draw(output = 'mpl', filename = str(i)+'Single')
        qc = qc.compose(initial_QEC[i],qubit_map[i],clrs[i])
        #qc.draw(output = 'mpl', filename = str(i)+'comp', fold = -1)

    #Load the current logical qubit states
    #Perform Transversal Logic Gates
    for op in operation_list:
        qc.barrier()
        if op[0] == 'CX':
            op_circ = transversal_CX()
            map = [i for i in logical_qubit_map[op[1]]]
            for i in logical_qubit_map[op[2]]:
                map.append(i)
            qc.append(op_circ,map)
        elif op[0] == 'CZ':
            op_circ = transversal_CZ()
            map = [i for i in logical_qubit_map[op[1]]]
            for i in logical_qubit_map[op[2]]:
                map.append(i)
            qc.append(op_circ,map)
        elif op[0] == 'X':
            op_circ = X_L
            qc.append(op_circ,logical_qubit_map[op[1]])
        elif op_circ == 'Z':
            op_circ = Z_L
            qc.append(op_circ,logical_qubit_map[op[1]])
        else:
            print("Unsupported Operation")
        qc.barrier()

    #Perform QEC
    
    #A 'logical qubit' must be passes into the QEC function to perform QEC on. As we want to just perform QEC on out existing circuit 
    #(i.e. append our circuit with the QEC procedure) we perform QEC on an empty logical qubit and append this to our circuit.
        empty_logical_qubits = [QuantumCircuit(5) for i in range(n_logical_qubits)]

        clrs = []
        QEC_circs = []
        for Lq in empty_logical_qubits:
            QC, clr = QEC(Lq)
            QEC_circs.append(QC)
            qc.add_register(clr)
            clrs.append(clr)
        
        for i in range(len(initial_QEC)):
            QEC_circs[i].decompose().draw(output = 'mpl', filename = str(i)+'SingleEmpty')
            qc = qc.compose(QEC_circs[i],qubit_map[i],clrs[i])
            #qc.draw(output = 'mpl', filename = str(i)+'compEmpty', fold = -1)






    return qc


#test = cycle(2,[['X',1],['CZ',0,1]])
#test.draw(output = 'mpl', filename = '00', fold = -1)

def measure_all_logical_qubits(quantum_circuit):

    logical_qubit_map = [[x for x in range(9*n+4,9*n+9)] for n in range(quantum_circuit.num_qubits//9)]
    for i in range(quantum_circuit.num_qubits//9):
        measure_reg = ClassicalRegister(5, 'meas'+str(i))
        quantum_circuit.add_register(measure_reg)
        quantum_circuit.measure(logical_qubit_map[i],measure_reg)

    return quantum_circuit

#test_2 = measure_all_logical_qubits(test)
#test_2.draw(output = 'mpl', filename = '000')

def measure_logical_qubits(quantum_circuit, logical_qubit_list):
    logical_qubit_map = [[x for x in range(9*n+4,9*n+9)] for n in range(quantum_circuit.num_qubits//9)]
    for i in range(len(logical_qubit_list)):
        measure_reg = ClassicalRegister(5, 'meas'+str(i))
        quantum_circuit.add_register(measure_reg)
        quantum_circuit.measure(logical_qubit_map[logical_qubit_list[i]],measure_reg)

    return quantum_circuit

#-------------------------------------------------------------Refactoring----------------------------------------------------------------

class LogicalQuantumCircuit():


    def __init__(self,n_logical_qubits,structure = '9n'):
        qubit_map = []
        if structure == '9n':
            logical_qubit_map = [[x for x in range(9*n+4,9*n+9)] for n in range(n_logical_qubits)]
            ancillary_qubit_map = [[x for x in range(9*n,9*n+4)] for n in range(n_logical_qubits)]
            qubit_map = [[x for x in range(9*n,9*n+9)] for n in range(n_logical_qubits)]
        elif structure == '4n5n':
            logical_qubit_map = [[x for x in range(4*n_logical_qubits,4*n_logical_qubits+5*n)] for n in range(n_logical_qubits)]
            ancillary_qubit_map = [[x for x in range(4*n,4*n+4)] for n in range(n_logical_qubits)]
            for i in range(n_logical_qubits):
                temp = []
                for element in ancillary_qubit_map[i]:
                    temp.append(element)
                for element in logical_qubit_map[i]:
                    temp.append(element)
                qubit_map.append(temp)
        else:
            print("Unsupported Structure")

        #Qubit Map Attributes
        self.qubit_map = qubit_map
        self.logical_qubit_map = logical_qubit_map
        self.ancillary_qubit_map = ancillary_qubit_map
        
        #Create the QuantumCircuit object to add the logical qubits to.
        qc = QuantumCircuit()
        for i in range(n_logical_qubits):
            qc.add_register(QuantumRegister(9, 'q'+str(i)))

        initial_logical_qubits = [logi_0() for i in range(n_logical_qubits)]
        
        #Initial error correction for single qubit errors in state preperation.
        i=0
        for Lq in initial_logical_qubits:
            QC, clr = QEC(Lq)
            qc.add_register(clr)
            qc = qc.compose(QC,qubit_map[i],clr)
            i+=1
    
        #Quantum Circuit Attribute
        self.qc = qc
        self.num_logical_qubits = n_logical_qubits

        return None
    
    #---------------------------------------------------------Logic Gates--------------------------------------------------------------------
    


    #-----------------------------Single Qubit Gates------------------------------------
    def x(self, logical_qubit):

        self.qc.barrier()
        op_circ = X_L
        self.qc.append(op_circ,self.logical_qubit_map[logical_qubit])
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def z(self, logical_qubit):
        
        self.qc.barrier()
        op_circ = X_L
        self.qc.append(op_circ,self.logical_qubit_map[logical_qubit])
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def y(self, logical_qubit):
        
        self.qc.barrier()
        op_circ = Y_L
        self.qc.append(op_circ,self.logical_qubit_map[logical_qubit])
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def K(self, logical_qubit, params):
        K = K_gate(params[0],params[1],params[2])
        map = [i for i in self.logical_qubit_map[logical_qubit]]
        self.qc.barrier()
        self.qc.append(K,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    

    #--------------------------------------------Two Qubit Contol Gates----------------------------------------
    
    def cx(self, logical_control_qubit, logical_target_qubit):

        self.qc.barrier()
        op_circ = transversal_CX()
        map = [i for i in self.logical_qubit_map[logical_control_qubit]]
        for i in self.logical_qubit_map[logical_target_qubit]:
            map.append(i)
        self.qc.append(op_circ,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def cz(self, logical_control_qubit, logical_target_qubit):

        self.qc.barrier()
        op_circ = transversal_CZ()
        map = [self.logical_qubit_map[logical_control_qubit],self.logical_qubit_map[logical_target_qubit]]
        self.qc.append(op_circ,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def cy(self, logical_control_qubit, logical_target_qubit):

        self.qc.barrier()
        op_circ = transversal_CY()
        map = [self.logical_qubit_map[logical_control_qubit],self.logical_qubit_map[logical_target_qubit]]
        self.qc.append(op_circ,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    #-------------------------------------------Helpers---------------------------------------------------
    

    def post_op_QEC(self):
        empty_logical_qubits = [QuantumCircuit(5) for i in range(self.num_logical_qubits)]
        i=0
        for Lq in empty_logical_qubits:
            QC, clr = QEC(Lq)
            self.qc.add_register(clr)
            self.qc = self.qc.compose(QC,self.qubit_map[i],clr)
            i+=1
    
        return None
    

    #------------------------------------------------------------Measurement----------------------------------------------------------------

    def measure_all(self):
        
        for i in range(self.num_logical_qubits):
            measure_reg = ClassicalRegister(5, 'meas'+str(i))
            self.qc.add_register(measure_reg)
            self.qc.measure(self.logical_qubit_map[i],measure_reg)

        return None
    
    def measure(self, logical_qubits):

        if type(logical_qubits) == int:
            measure_reg = ClassicalRegister(5, 'meas'+str(logical_qubits))
            self.qc.add_register(measure_reg)
            self.qc.measure(self.logical_qubit_map[logical_qubits],measure_reg)
        else:
            for i in range(len(logical_qubits)):
                measure_reg = ClassicalRegister(5, 'meas'+str(i))
                self.qc.add_register(measure_reg)
                self.qc.measure(self.logical_qubit_map[logical_qubits[i]],measure_reg)

        return None




#Test if action of X_L is as expected and general framwork works.

test = LogicalQuantumCircuit(1)
test.x(0)
test.measure_all()
test.qc.draw(output = 'mpl', filename = '3', fold = -1)

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram

simulator = AerSimulator()
circ = transpile(test.qc, simulator)

result = simulator.run(circ).result() #Run and get counts
print(result)
counts = result.get_counts(circ)
print(counts)
plot_histogram(counts, title='Logical X counts on |0>_L', filename = 'X_L_Test') 

#Works, but measurement in the computational basis isnt ideal. Need to figure out a way to measure in the 'logical qubit basis'. 
#Encoding projector is a singular matrix so cant apply an inverse to re-map to the computational basis before measurement. Need to think about this. 
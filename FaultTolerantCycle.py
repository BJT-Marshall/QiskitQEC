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
    """
    A quantum circuit framework for performing transversal operations and single qubit error correction in the 5-qubit-code. Methods and attributes meant to closely follow the framework of a 
    :python:`qiskit.QuantumCircuit` object for ease of use.

    Attributes:
    -----------

    qubit_map : list
        A map grouping the physical qubits into the logical and ancillary qubits they represent. 
    logical_qubit_map : list
        A map grouping the physical qubits into the logical qubits they represent.
    ancillary_qubit_map : list
        A map grouping the physical qubits into the ancillary qubit sets they represent.
    qc : QuantumCircuit
        The :python:`qiskit.QuantumCircuit` object representing the logical quantum circuit. It is this attribute object that is interfaced with to perform logical operations, error correction and measurements.
    num_logical_qubits : int
        The number of logical qubits encoded in the 5-qubit-code contained within the :python:`LogicalQuantumCircuit`

    Methods:
    --------

    Logical Operations:
    -------------------
    
    x : Applies the logical :math:`X_L` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the single qubit :math:`X` gate to a :python:`qiskit.QuantumCircuit` object.
    y : Applies the logical :math:`Y_L` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the single qubit :math:`Y` gate to a :python:`qiskit.QuantumCircuit` object.
    z : Applies the logical :math:`Z_L` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the single qubit :math:`Z` gate to a :python:`qiskit.QuantumCircuit` object.
    k : Applies the logical :math:`K_{s_x,s_y,s_z}` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`.
    cx : Applies the logical :math:`CX_L` quantum gate to a pair of indexed logical qubits of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the two qubit :math:`CX` gate to a :python:`qiskit.QuantumCircuit` object.
    cy : Applies the logical :math:`CY_L` quantum gate to a pair of indexed logical qubits of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the two qubit :math:`CY` gate to a :python:`qiskit.QuantumCircuit` object.
    cz : Applies the logical :math:`CZ_L` quantum gate to a pair of indexed logical qubits of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the two qubit :math:`CZ` gate to a :python:`qiskit.QuantumCircuit` object.

    Measurement:
    ------------

    measure_all : Measures all logical qubits of :python:`LogicalQuantumCircuit.qc` out to classical registers in the computational basis with labels :python:`'meas+str(i)'` for logical qubit index :python:`i`. Analogous to the :python:`QuantumCircuit.measure_all()` method of the standard Qiskit :python:`qiskit.QuantumCircuit` object.
    measure : Measures the indexed logical qubits of :python:`LogicalQuantumCircuit.qc` out to classical registers in the computational basis with labels :python:`'meas+str(i)'` for logical qubit index :python:`i`. Analogous to the :python:`QuantumCircuit.measure(i)` method of the standard Qiskit :python:`qiskit.QuantumCircuit` object.

    Helper Methods:
    ---------------

    post_op_QEC : Performs single qubit error correction on the current state of the logical qubits of :python:`LogicalQuantumCircuit.qc`. Called during the application of logical operations such as :python:`LogicalQuantumCircuit.x(0)` to correct for any single qubit errors accumulated during the operation.

    """


    def __init__(self,n_logical_qubits,structure = '9n'):
        """
        Generates an instance of :python:`LogicalQuantumCircuit` containing :python:`n_logical_qubits` with all the required attributes.
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param n_logical_qubits: Number of logical qubits required in :python:`LogicalQuantumCircuit.qc`
        :type n_logical_qubits: int
        :param structure: Optional string defining the structure of :python:`LogicalQuantumCircuit.qc` with respect to the ordering of logical qubits and sets of ancillary qubits used for single qubit error correction. Default: :python:`9n`.
        :type structure: str
        """
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
    
        #Quantum Circuit Attributes
        self.qc = qc
        self.num_logical_qubits = n_logical_qubits

        return None
    
    #---------------------------------------------------------Logic Gates--------------------------------------------------------------------
    


    #-----------------------------Single Qubit Gates------------------------------------
    def x(self, logical_qubit):
        """
        Applies the logical :math:`X_L` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the single qubit :math:`X` gate to a :python:`qiskit.QuantumCircuit` object.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(1)
        logical_circuit.x(0) #Applies the logical X operation to the 1st qubit of the LogicalQuantumCircuit.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit_1.x(3) #Applies the logical X operation to the 4th qubit of the LogicalQuantumCircuit.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_qubit: The index of the logical qubit to which the operation should be applied.
        :type logical_qubit: int
        """

        self.qc.barrier()
        op_circ = X_L
        self.qc.append(op_circ,self.logical_qubit_map[logical_qubit])
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def z(self, logical_qubit):
        """
        Applies the logical :math:`Z_L` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the single qubit :math:`Z` gate to a :python:`qiskit.QuantumCircuit` object.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(1)
        logical_circuit.z(0) #Applies the logical Z operation to the 1st qubit of the LogicalQuantumCircuit.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit_1.z(3) #Applies the logical Z operation to the 4th qubit of the LogicalQuantumCircuit.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_qubit: The index of the logical qubit to which the operation should be applied.
        :type logical_qubit: int
        """
        
        self.qc.barrier()
        op_circ = X_L
        self.qc.append(op_circ,self.logical_qubit_map[logical_qubit])
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def y(self, logical_qubit):
        """
        Applies the logical :math:`Y_L` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the single qubit :math:`Y` gate to a :python:`qiskit.QuantumCircuit` object.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(1)
        logical_circuit.y(0) #Applies the logical Y operation to the 1st qubit of the LogicalQuantumCircuit.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit_1.y(3) #Applies the logical Y operation to the 4th qubit of the LogicalQuantumCircuit.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_qubit: The index of the logical qubit to which the operation should be applied.
        :type logical_qubit: int
        """
        
        self.qc.barrier()
        op_circ = Y_L
        self.qc.append(op_circ,self.logical_qubit_map[logical_qubit])
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def K(self, logical_qubit, params):
        """
        Applies the logical :math:`K_{s_x,s_y,s_z}` quantum gate to an indexed logical qubit of the :python:`LogicalQuantumCircuit.qc`.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(1)
        logical_circuit.k(0, [1,1,1]) #Applies the logical K_{+,+,+} operation to the 1st qubit of the LogicalQuantumCircuit.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit_1.k(4, [-1,1,-1]) #Applies the logical K_{-,+,-} operation to the 4th qubit of the LogicalQuantumCircuit.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_qubit: The index of the logical qubit to which the operation should be applied.
        :type logical_qubit: int
        :param params: List of parameters [sx,sy,sz] used in construction of the :math:`K_{s_x,s_y,s_z}` gate. See documentation of :python:`TransversalLogicGates.K_gate` for more details.
        :type params: list
        """
        K = K_gate(params[0],params[1],params[2])
        map = [i for i in self.logical_qubit_map[logical_qubit]]
        self.qc.barrier()
        self.qc.append(K,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    

    #--------------------------------------------Two Qubit Contol Gates----------------------------------------
    
    def cx(self, logical_control_qubit, logical_target_qubit):
        """
        Applies the logical :math:`CX_L` quantum gate to a pair of indexed logical qubits of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the two qubit :math:`CX` gate to a :python:`qiskit.QuantumCircuit` object.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(2)
        logical_circuit.cx(0,1) #Applies the logical CX operation to the LogicalQuantumCircuit with the 1st qubit acting as the control qubit and the 2nd as the target.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit.cx(1,3) #Applies the logical CX operation to the LogicalQuantumCircuit with the 2nd qubit acting as the control qubit and the 4th as the target.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_control_qubit: The index of the control logical qubit to which the operation should be applied.
        :type logical_control_qubit: int
        :param logical_control_qubit: The index of the target logical qubit to which the operation should be applied.
        :type logical_control_qubit: int
        """

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
        """
        Applies the logical :math:`CZ_L` quantum gate to a pair of indexed logical qubits of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the two qubit :math:`CZ` gate to a :python:`qiskit.QuantumCircuit` object.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(2)
        logical_circuit.cz(0,1) #Applies the logical CZ operation to the LogicalQuantumCircuit with the 1st qubit acting as the control qubit and the 2nd as the target.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit.cz(1,3) #Applies the logical CZ operation to the LogicalQuantumCircuit with the 2nd qubit acting as the control qubit and the 4th as the target.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_control_qubit: The index of the control logical qubit to which the operation should be applied.
        :type logical_control_qubit: int
        :param logical_control_qubit: The index of the target logical qubit to which the operation should be applied.
        :type logical_control_qubit: int
        """

        self.qc.barrier()
        op_circ = transversal_CZ()
        map = [self.logical_qubit_map[logical_control_qubit],self.logical_qubit_map[logical_target_qubit]]
        self.qc.append(op_circ,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    def cy(self, logical_control_qubit, logical_target_qubit):
        """
        Applies the logical :math:`CY_L` quantum gate to a pair of indexed logical qubits of the :python:`LogicalQuantumCircuit.qc`. Analogous to the application of the two qubit :math:`CY` gate to a :python:`qiskit.QuantumCircuit` object.

        Example:
        ```python
        logical_circuit = LogicalQuantumCircuit(2)
        logical_circuit.cy(0,1) #Applies the logical CY operation to the LogicalQuantumCircuit with the 1st qubit acting as the control qubit and the 2nd as the target.

        logical_circuit_1 = LogicalQuantumCircuit(5)
        logical_circuit.cy(1,3) #Applies the logical CY operation to the LogicalQuantumCircuit with the 2nd qubit acting as the control qubit and the 4th as the target.
        ```
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_control_qubit: The index of the control logical qubit to which the operation should be applied.
        :type logical_control_qubit: int
        :param logical_control_qubit: The index of the target logical qubit to which the operation should be applied.
        :type logical_control_qubit: int
        """

        self.qc.barrier()
        op_circ = transversal_CY()
        map = [self.logical_qubit_map[logical_control_qubit],self.logical_qubit_map[logical_target_qubit]]
        self.qc.append(op_circ,map)
        self.qc.barrier()
        self.post_op_QEC()

        return None
    
    #-------------------------------------------Helpers---------------------------------------------------
    

    def post_op_QEC(self):
        """
        Performs single qubit error correction on the current state of the logical qubits of :python:`LogicalQuantumCircuit.qc`. Called during the application of logical operations such as :python:`LogicalQuantumCircuit.x(0)`
        to correct for any single qubit errors accumulated during the operation.
        
        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        """
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
        """
        Measures all logical qubits of :python:`LogicalQuantumCircuit.qc` out to classical registers in the computational basis with labels :python:`'meas+str(i)'` for logical qubit index :python:`i`.
        Analogous to the :python:`QuantumCircuit.measure_all()` method of the standard Qiskit :python:`qiskit.QuantumCircuit` object.

        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        """
        
        for i in range(self.num_logical_qubits):
            measure_reg = ClassicalRegister(5, 'meas'+str(i))
            self.qc.add_register(measure_reg)
            self.qc.measure(self.logical_qubit_map[i],measure_reg)

        return None
    
    def measure(self, logical_qubits):
        """
        Measures the indexed logical qubits of :python:`LogicalQuantumCircuit.qc` out to classical registers in the computational basis with labels :python:`'meas+str(i)'` for logical qubit index :python:`i`.
        Analogous to the :python:`QuantumCircuit.measure(i)` method of the standard Qiskit :python:`qiskit.QuantumCircuit` object.

        :param self: LogicalQuantumCircuit
        :type self: LogicalQuantumCircuit
        :param logical_qubits: The logical qubit index or indices that should be measured.
        :type logical_qubits: list or int
        """

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




#Testing:

test = False

if test:
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
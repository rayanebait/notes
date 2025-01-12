from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.all import CC

from quantum.qubits import Qubits
from quantum.qubit import Qubit, ket, index
from quantum.qubit import normalize_positions
from quantum.qubit import decompose_qubit, recompose_qubit

class QuantumGate_generic(SageObject):
    def __init__(self, length, name=""):
        self._name = str(name)
        self._length = length
        self._type = "gate"

    def __call__(self, qubit, positions=None):
        if not isinstance(qubit, Qubit):
            raise TypeError("argument must be a qubit")
        s = self._length
        t = qubit.length() - s
        J, I = normalize_positions(s, t, positions, None)
        Q = decompose_qubit(qubit, I, J)
        R = [ self._call_(q) for q in Q ]
        return recompose_qubit(R, I, J)

    def __hash__(self):
        return hash((self._length, self._name))

    def __repr__(self):
        s = "Quantum "
        if len(self._name) > 0: s += self._name + " "
        s += "%s acting on %s qubit" % (self._type, self._length)
        if self._length > 1: s += "s"
        return s

    def _call_(self, q):
        raise NotImplementedError

    def length(self):
        return self._length

    def matrix(self):
        n = self._length
        parent = Qubits(n)
        rows = [ ]
        for u in range(2**n):
            q = parent(ket(u,n))
            q = self._call_(q)
            rows.append(q.list())
        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(CC,2**n)(rows)

    def complexity(self):
        return { self: 1 }


# Standard gates
################

# Hadamard gate

class HadamardGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 1, "Hadamard")
        self._sq = CC(0.5).sqrt()

    def _call_(self, q):
        return q.parent()([ self._sq*(q[0] + q[1]), self._sq*(q[0] - q[1]) ])

HadamardGate = HadamardGate_class()


# Phase shift

class PhaseShiftGate(QuantumGate_generic, UniqueRepresentation):
    def __init__(self, angle):
        QuantumGate_generic.__init__(self, 1, "phase shift")
        self._angle = angle
        self._exp = CC((2*ZZ(-1).sqrt()*self._angle).exp())

    def __hash__(self):
        return hash((self._length, self._name, self._angle))

    def angle(self):
        return self._angle

    def _call_(self, q):
        return q.parent()([ q[0], self._exp*q[1] ])


# Swap gate

class SwapGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 2, "swap")

    def _call_(self, q):
        return q.parent()([ q[0], q[2], q[1], q[3] ])

SwapGate = SwapGate_class()


# Pauli gates

class PauliXGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 1, "Pauli X")

    def _call_(self, q):
        return q.parent()([ q[1], q[0] ])

class PauliYGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 1, "Pauli Y")
        self._i = CC(-1).sqrt()

    def _call_(self, q):
        return q.parent()([ -self._i*q[1], self._i*q[0] ])

class PauliZGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 1, "Pauli Z")

    def _call_(self, q):
        return q.parent()([ q[0], -q[1] ])

PauliXGate = NOTGate = PauliXGate_class()
PauliYGate = PauliYGate_class()
PauliZGate = PauliZGate_class()


# Controlled gates
##################

class ControlledGate(QuantumGate_generic):
    def __init__(self, gate):
        QuantumGate_generic.__init__(self, gate.length() + 1)
        self._name = "controlled " + gate._name
        self._gate = gate

    def __hash__(self):
        return hash((self._length, self._name, self._gate))

    def _call_(self, q):
        I = [0]
        J = range(1, self._length)
        Q = decompose_qubit(q, I, J)
        R = [ Q[0], self._gate(Q[1]) ]
        return recompose_qubit(R, I, J)


class CNOTGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 2, "CNOT")

    def _call_(self, q):
        return q.parent()([ q[0], q[1], q[3], q[2] ])

class ToffoliGate_class(QuantumGate_generic, UniqueRepresentation):
    def __init__(self):
        QuantumGate_generic.__init__(self, 3, "Toffoli")

    def _call_(self, q):
        return q.parent()([ q[0], q[1], q[2], q[3], q[4], q[5], q[7], q[6] ])

CNOTGate = CNOTGate_class()
ToffoliGate = ToffoliGate_class()


# Gates from boolean functions
##############################

class QuantumGateFromBooleanFunction(QuantumGate_generic):
    def __init__(self, n, m, f, type_input='integer', type_output='integer'):
        self._n = n
        self._m = m
        self._mask = ZZ(2**m - 1)
        self._function = f
        if not type_input in [ 'integer', 'string' ]:
            raise ValueError("type_input must be 'integer' or 'string'")
        self._convert_input = (type_input == 'string')
        if not type_output in [ 'integer', 'string' ]:
            raise ValueError("type_output must be 'integer' or 'string'")
        self._convert_output = (type_output == 'string')
        QuantumGate_generic.__init__(self, n+m)

    def _call_(self, q):
        n = self._n
        m = self._m
        mask = self._mask
        convert_input = self._convert_input
        convert_output = self._convert_output
        coords = [ ]
        for x in range(2**n):
            if convert_input:
                arg = ket(x, n)
            else:
                arg = ZZ(x)
            value = self._function(arg)
            if convert_output:
                value = index(value)
            else:
                value = ZZ(value)
            value &= mask
            u = x << m
            for y in range(2**m):
                j = u + value.__xor__(y)
                coords.append(q[j])
        return q.parent()(coords)



# Circuits
##########

class QuantumCircuit(QuantumGate_generic):
    def __init__(self, length, gates, name=""):
        QuantumGate_generic.__init__(self, length, name)
        self._gates = gates  # do checkings
        self._type = "circuit"

    def _call_(self, q):
        for gate, positions in self._gates:
            q = gate(q, positions)
        return q

    def is_deterministic(self):
        for gate, _ in self._gates:
            if not gate.is_deterministic(): return False
        return True

    def complexity(self):
        c = { }
        for gate, _ in self._gates:
            if gate in c:
                c[gate] += 1
            else:
                c[gate] = 1
        return c

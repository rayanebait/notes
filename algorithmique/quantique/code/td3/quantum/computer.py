from sage.structure.sage_object import SageObject
from sage.arith.srange import srange
from sage.misc.prandom import random

from sage.rings.integer_ring import ZZ
from sage.rings.all import CC

from quantum.qubits import Qubits
from quantum.qubit import ket
from quantum.qubit import normalize_positions
from quantum.qubit import decompose_qubit, recompose_qubit
from quantum.gate import QuantumGate_generic
from quantum.gate import HadamardGate, NOTGate, CNOTGate, ToffoliGate


class QuantumRegister(SageObject):
    def __init__(self, computer, address, size):
        self._computer = computer
        self._address = address
        self._size = size

    def __repr__(self):
        s = "Quantum register of %s qubit" % self._size
        if self._size > 1:
            s += "s"
        return s

    def computer(self):
        return self._computer

    def address(self):
        return self._address

    def size(self):
        return self._size

    def __getitem__(self, key, reverse=True):
        if isinstance(key, slice):
            if key.step is not None:
                raise TypeError("specifying an increment is not allowed")
            start = key.start
            stop = key.stop
        else:
            start = ZZ(key)
            stop = start + 1
        if start < 0 or stop > self._size:
            raise ValueError("segmentation fault")
        if reverse:
            return QuantumRegister(self._computer, self._address+self._size-stop, stop-start)
        else:
            return QuantumRegister(self._computer, self._address+start, stop-start)



class QuantumComputer(SageObject):
    def __init__(self, name=None):
        self._length = 0
        self._maxlength = 0
        self._state = Qubits(0)('')
        if name is None:
            self._name = None
        else:
            self._name = str(name)
        self._positions = [ ]
        self._operations = 0
        self._oracles = 0

    def __repr__(self):
        s = "Quantum computer"
        if self._name is not None:
            s += " '%s'" % self._name
        return s

    def internal_state(self):
        return self._state

    def used_memory(self):
        return self._length

    def cputime(self):
        return self._operations

    def walltime(self):
        return self._operations + self._oracles

    def _get_positions(self, *registers):
        I = [ ]
        for register in registers:
            address = register.address()
            size = register.size()
            if address+size > len(self._positions):
                raise ValueError("segmentation fault")
            I += self._positions[address:address+size]
        if None in I:
            raise ValueError("segmentation fault")
        return I

    def _measure(self, I, J):
        self._operations += len(I)
        self._state = self._state.normalize()
        Q = decompose_qubit(self._state, I, J)
        rand = random()
        thresold = 0
        for u in range(2**len(I)):
            thresold += Q[u].normsquare()
            if rand <= thresold:
                return u, Q[u].normalize()

    def malloc(self, size):
        address = len(self._positions)
        self._positions += [ ZZ(i) for i in range(self._length, self._length+size) ]
        self._length += size
        self._state = self._state.tensor_with_ket(size*"0")
        return QuantumRegister(self, address, size)

    def free(self, register):
        I = self._get_positions(register)
        J = [ i for i in range(self._length) if i not in I ]
        u, self._state = self._measure(I, J)
        address = register.address()
        size = register.size()
        self._length -= size
        for addr in range(address, address+size):
            self._positions[addr] = None
        for addr in range(address+size, len(self._positions)):
            if self._positions[addr] is not None:
                self._positions[addr] -= size
        while len(self._positions) > 0 and self._positions[-1] is None:
            del self._positions[-1]
        return u

    def reset(self, register):
        I = self._get_positions(register)
        J = [ i for i in range(self._length) if i not in I ]
        _, q = self._measure(I, J)
        size = register.size()
        self._state = q.tensor_with_ket(size*"0", J, I)

    def set(self, register, value):
        I = self._get_positions(register)
        J = [ i for i in range(self._length) if i not in I ]
        _, q = self._measure(I, J)
        size = register.size()
        self._state = q.tensor_with_ket(ket(value,size), J, I)

    def measure(self, register):
        I = self._get_positions(register)
        J = [ i for i in range(self._length) if i not in I ]
        u, q = self._measure(I, J)
        size = register.size()
        self._state = q.tensor_with_ket(ket(u,size), J, I)
        return ZZ(u)

    def hadamard(self, x):
        if not isinstance(x, QuantumRegister):
            raise ValueError("x must be a quantum register")
        self._operations += x.size()
        for i in range(x.size()):
            I = self._get_positions(x[i])
            self._state = HadamardGate(self._state, I)
        self._state = self._state.normalize()

    def X(self, x):
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        self._operations += x.size()
        sx = self._length - 1 - self._get_positions(x)[-1]
        mask = ((1 << x.size()) - 1) << sx
        coeffs = self._state.list()
        coeffs = [ coeffs[i ^ mask] for i in range(len(coeffs)) ]
        self._state = Qubits(self._length)(coeffs)

    def CX(self, c, x):
        if not (isinstance(c, QuantumRegister) and c.computer() == self and c.size() == 1):
            raise ValueError("c2 must be a quantum register of size 1 for this computer")
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        self._operations += 1
        sc = self._length - 1 - self._get_positions(c)[0]
        sx = self._length - 1 - self._get_positions(x)[-1]
        mx = (1 << x.size()) - 1
        coeffs = self._state.list()
        coeffs = [ self._state[i ^ ((((i >> sc) & 1)*mx) << sx)] for i in range(len(coeffs)) ]
        self._state = Qubits(self._length)(coeffs)

    def swap(self, x, y):
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        if not (isinstance(y, QuantumRegister) and y.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        size = x.size()
        if y.size() != size:
            raise ValueError("the two registers must have the same size")
        self._operations += 1
        sx = self._length - 1 - self._get_positions(x)[-1]
        sy = self._length - 1 - self._get_positions(y)[-1]
        ds = sx - sy
        if ds < 0:
            ds = -ds
            sx, sy = sy, sx
        mx = ((1 << size) - 1) << sx
        my = ((1 << size) - 1) << sy
        coeffs = self._state.list()
        newcoeffs = [ ]
        for i in range(len(coeffs)):
            vx = i & mx
            vy = i & my
            j = (i - vx - vy) + (vx >> ds) + (vy << ds)
            newcoeffs.append(coeffs[j])
        self._state = Qubits(self._length)(newcoeffs)

    def CCX(self, c1, c2, x):
        if not (isinstance(c1, QuantumRegister) and c1.computer() == self and c1.size() == 1):
            raise ValueError("c1 must be a quantum register of size 1 for this computer")
        if not (isinstance(c2, QuantumRegister) and c2.computer() == self and c2.size() == 1):
            raise ValueError("c2 must be a quantum register of size 1 for this computer")
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        self._operations += 1
        s1 = self._length - 1 - self._get_positions(c1)[0]
        s2 = self._length - 1 - self._get_positions(c2)[0]
        sx = self._length - 1 - self._get_positions(x)[-1]
        mx = (1 << x.size()) - 1
        coeffs = self._state.list()
        coeffs = [ coeffs[i ^ ((((i >> s1) & (i >> s2) & 1)*mx) << sx)] for i in range(len(coeffs)) ]
        self._state = Qubits(self._length)(coeffs)

    def phase_shift(self, x, ang):
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        self._operations += 1
        I = CC(-1).sqrt()
        scalar = (I*ang).exp()
        I = self._get_positions(x)
        I = [ self._length - s - 1 for s in I ]
        coeffs = self._state.list()
        coeffs = [ coeffs[i] * scalar**sum((i >> s) & 1 for s in I) for i in range(len(coeffs)) ]
        self._state = Qubits(self._length)(coeffs)

    def apply(self, f, *registers):
        pos = [ ]
        masks = [ ]
        if len(registers) == 1 and isinstance(registers[0], (list, tuple)):
            registers = registers[0]
        elif len(registers) == 2 and isinstance(registers[0], (list, tuple)):
            registers = registers[0] + [ registers[1] ]
        for x in registers:
            if not (isinstance(x, QuantumRegister) and x.computer() == self):
                raise ValueError("you must pass in quantum registers for this computer")
            pos.append(self._length - 1 - self._get_positions(x)[-1])
            masks.append((1 << x.size()) - 1)
        self._oracles += 1
        coeffs = self._state.list()
        newcoeffs = [ ]
        for i in range(len(coeffs)):
            args = [ ]
            for m in range(len(registers) - 1):
                args.append((i >> pos[m]) & masks[m])
            j = i ^ ((f(*args) & masks[-1]) << pos[-1])
            newcoeffs.append(coeffs[j])
        self._state = Qubits(self._length)(newcoeffs)

    def ciaddmod(self, c, x, a, N):
        if not (isinstance(c, QuantumRegister) and c.computer() == self and c.size() == 1):
            raise ValueError("c must be a quantum register of size 1 for this computer")
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        sc = self._length - 1 - self._get_positions(c)[0]
        sx = self._length - 1 - self._get_positions(x)[-1]
        coeffs = self._state.list()
        maskc = 1 << sc
        mask1 = ((1 << x.size()) - 1) << sx
        mask2 = ~mask1
        aa = a << sx
        NN = N << sx
        newcoeffs = [ ]
        for i in range(len(coeffs)):
            if i & maskc:
                v1 = i & mask1
                v2 = i & mask2
                if v1 < NN:
                    v1 = (v1 - aa) % NN
                j = v1 | v2
            else:
                j = i
            newcoeffs.append(coeffs[j])
        self._state = Qubits(self._length)(newcoeffs)

    def cimulmod(self, c, x, a, N):
        if not (isinstance(c, QuantumRegister) and c.computer() == self and c.size() == 1):
            raise ValueError("c must be a quantum register of size 1 for this computer")
        if not (isinstance(x, QuantumRegister) and x.computer() == self):
            raise ValueError("x must be a quantum register for this computer")
        d, b, _ = a.xgcd(N)
        if d != 1:
            raise ValueError("the multiplier must be coprime with the modulus")
        sc = self._length - 1 - self._get_positions(c)[0]
        sx = self._length - 1 - self._get_positions(x)[-1]
        coeffs = self._state.list()
        maskc = 1 << sc
        mask1 = ((1 << x.size()) - 1) << sx
        mask2 = ~mask1
        NN = N << sx
        newcoeffs = [ ]
        for i in range(len(coeffs)):
            if i & maskc:
                v1 = i & mask1
                v2 = i & mask2
                if v1 < NN:
                    v1 = (v1 * b) % NN
                j = v1 | v2
            else:
                j = i
            newcoeffs.append(coeffs[j])
        self._state = Qubits(self._length)(newcoeffs)

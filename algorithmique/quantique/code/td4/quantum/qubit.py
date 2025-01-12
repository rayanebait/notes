from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.rings.all import CC
from sage.misc.prandom import random


concat = "".join

def index(x):
    if x == "": return 0
    return ZZ(x, base=2)

def ket(u, n):
    if n == 0: return ""
    x = ZZ(u).str(base=2)
    return (n - len(x))*"0" + x

def check_positions(pos, s, n):
    for i in pos:
        if not (i in ZZ and i >= 0 and i < n):
            raise ValueError("incorrect positions")
    if len(pos) != s:
        raise ValueError("the number of positions does not agree with the length of the qubit")
    if len(set(pos)) != s:
        raise ValueError("positions contain duplicates items")

def normalize_positions(s, t, I, J):
    n = s + t
    if I is None and J is None:
        return range(s), range(s,n)
    swap = False
    if I is None:
        swap = True
        I, J = J, I
        s, t = t, s
    check_positions(I, s, n)
    if J is None:
        J = [ i for i in range(n) if i not in I ]
    else:
        check_positions(J, t, n)
        if len(set(I + J)) != n:
            raise ValueError("positions intersect")
    if swap:
        return J, I
    else:
        return I, J


class Qubit(Element):
    def __init__(self, parent, deftn=None):
        Element.__init__(self, parent)
        n = parent.length()
        if deftn is None:
            self._coordinates = [CC(1)] + (2**n-1) * [CC(0)]
        elif isinstance(deftn, str):
            if len(deftn) > n:
                raise ValueError("cannot convert a |%s> into a %s-qubit" % (deftn, n))
            deftn += (n - len(deftn)) * "0"
            self._coordinates = (2**n) * [CC(0)]
            self._coordinates[index(deftn)] = CC(1)
        elif isinstance(deftn, (list,tuple)) and len(deftn) == 2**n:
            self._coordinates = [ CC(c) for c in deftn ]
        elif isinstance(deftn, Qubit):
            m = deftn.length()
            if m > n:
                raise ValueError("cannot convert a %s-qubit into a %s-qubit" % (m, n))
            self._coordinates = (2**n) * [CC(0)]
            for i in range(2**m):
                self._coordinates[i << (n-m)] = deftn[i]
        else:
            raise TypeError

    def length(self):
        return self.parent().length()

    def _repr_(self):
        n = self.length()
        s = ""
        for u in range(2**n):
            c = self._coordinates[u]
            if c == 0: continue
            if c == 1:
                s += " + |%s>" % ket(u,n)
            elif c == -1:
                s += " - |%s>" % ket(u,n)
            elif c._is_atomic():
                s += " + %s*|%s>" % (c, ket(u,n))
            elif (-c)._is_atomic():
                s += " - %s*|%s>" % (-c, ket(u,n))
            else:
                s += " + (%s)*|%s>" % (c, ket(u,n))
        if s[:3] == "-":
            s = "-" + s[3:]
        else:
            s = s[3:]
        return s

    def _latex_(self):
        n = self.length()
        s = ""
        for u in range(2**n):
            c = self._coordinates[u]
            if c == 0: continue
            if c == 1:
                s += " + "
            elif c == -1:
                s += " - "
            elif c._is_atomic():
                s += " + %s\\cdot" % c._latex_()
            elif (-c)._is_atomic():
                s += " - %s\\cdot" % (-c)._latex_()
            else:
                s += " + \\left(%s\\right)\\cdot" % c._latex_()
            s += "\\left|%s\\right>" % ket(u,n)
        if s[:3] == "-":
            s = "-" + s[3:]
        else:
            s = s[3:]
        return s

    def list(self):
        return self._coordinates

    def dict(self, nonzero=False):
        n = self.length()
        if nonzero:
            return { ket(u,n): self._coordinates[u] for u in range(2**n) if self._coordinates[u] != 0 }
        else:
            return { ket(u,n): self._coordinates[u] for u in range(2**n) }

    def __getitem__(self, key):
        if isinstance(key, str):
            n = self.length()
            if len(key) <= n:
               key += (n - len(key)) * "0"
               return self._coordinates[index(key)]
        elif key in ZZ and key >= 0:
            return self._coordinates[key]
        raise IndexError("key '%s' is invalid for a %s-qubit" % (key, n))

    def norm(self):
        return sum(c.norm() for c in self._coordinates).sqrt()

    def normsquare(self):
        return sum(c.norm() for c in self._coordinates)

    def normalize(self):
        norm = self.norm()
        return self.parent()([ c/norm for c in self._coordinates ])

    def _add_(self, other):
        parent = self.parent()
        n = parent.length()
        return parent([ self._coordinates[u] + other._coordinates[u] for u in range(2**n) ])

    def _sub_(self, other):
        parent = self.parent()
        n = parent.length()
        return parent([ self._coordinates[u] - other._coordinates[u] for u in range(2**n) ])

    def tensor(self, other, I=None, J=None):
        from quantum.qubits import Qubits
        if not isinstance(other, Qubit):
            raise TypeError
        s = self.length()
        t = other.length()
        n = s + t
        I, J = normalize_positions(s, t, I, J)
        q = [ ]
        for u in range(2**n):
            z = ket(u,n)
            x = concat(z[i] for i in I)
            y = concat(z[i] for i in J)
            q.append(self[x] * other[y])
        return Qubits(n)(q)

    def tensor_with_ket(self, y, I=None, J=None):
        from quantum.qubits import Qubits
        if not isinstance(y, str):
            raise TypeError
        s = self.length()
        t = len(y)
        n = s + t
        I, J = normalize_positions(s, t, I, J)
        q = (2**n) * [ CC(0) ]
        z = n * [ "" ]
        for k in range(t):
            z[J[k]] = y[k]
        for u in range(2**s):
            x = ket(u,s)
            for k in range(s):
                z[I[k]] = x[k]
            q[index(concat(z))] = self._coordinates[u]
        return Qubits(n)(q)

    def __mul__(self, other):
        return self.tensor(other)


# qubit decomposition
# (helper functions)

def decompose_qubit(q, I, J):
    from quantum.qubits import Qubits
    s = len(I); t = len(J)
    z = (s+t) * [""]
    Q = [ ]
    parent = Qubits(t)
    for u in range(2**s):
        x = ket(u,s)
        for k in range(s):
            z[I[k]] = x[k]
        qx = [ ]
        for v in range(2**t):
            y = ket(v,t)
            for k in range(t):
                z[J[k]] = y[k]
            qx.append(q[concat(z)])
        Q.append(parent(qx))
    return Q

def recompose_qubit(Q, I, J):
    from quantum.qubits import Qubits
    s = len(I); t = len(J)
    n = s + t
    q = [ ]
    for u in range(2**n):
        z = ket(u,n)
        x = concat(z[i] for i in I)
        y = concat(z[j] for j in J)
        q.append(Q[index(x)][y])
    return Qubits(n)(q)

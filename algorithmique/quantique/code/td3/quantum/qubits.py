from sage.structure.parent import Parent
from sage.rings.all import CC


class Qubits(Parent):
    def __init__(self, n):
        Parent.__init__(self)
        self._n = n

    def _element_constructor_(self, *args, **kwargs):
        from quantum.qubit import Qubit
        return Qubit(self, *args, **kwargs)

    def __repr__(self):
        return "Set of %s-qubits" % self._n

    def length(self):
        return self._n

    def random_element(self, normalize=True):
        q = self([ CC.random_element() for _ in range(2**self._n) ])
        if normalize:
            return q.normalize()
        else:
            return q

    def uniform_superposition(self):
        c = CC(2 ** (-self._n / 2))
        return self([ c for _ in range(2**self._n) ])

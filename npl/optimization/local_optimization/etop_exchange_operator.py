import itertools

import numpy as np
from sortedcontainers import SortedKeyList


class ETOPExchangeOperator:
    """Guided chemical-ordering exchange for an arbitrary number of metals.

    For each ordered species pair ``(i, j)`` (i != j) keeps a ``SortedKeyList`` of
    atoms currently of species ``i`` ranked by the flip-gain estimate
    ``g_{i->j}(A) = coef . (f_A^(j) - f_A^(i))``, where ``f_A^(s)`` is A's ETOP
    per-atom feature as if A were species ``s``. A canonical swap of an i-atom A
    with a j-atom B has estimated energy change ``g_{i->j}(A) + g_{j->i}(B)``;
    this is a proposal heuristic (it ignores the i-j cross-bond and shared-neighbor
    coupling). Exact accept/reject is the driver's responsibility.

    Requires a linear energy model: per-atom energy = ``coef . atom_feature``.
    """

    def __init__(self, coefficients, feature_classifier):
        self.coeffs = np.asarray(coefficients)
        self.classifier = feature_classifier
        self.symbols = list(feature_classifier.symbols)
        assert len(self.coeffs) == feature_classifier.n_features, (
            "coefficient vector length {} != classifier.n_features {}".format(
                len(self.coeffs), feature_classifier.n_features))

        self.flip_gain = {}     # (sym_i, sym_j) -> {atom_index: gain}
        self.indices = {}       # (sym_i, sym_j) -> SortedKeyList(atom_index)
        self.atom_symbol = {}   # atom_index -> current symbol

    def _pairs(self):
        for i in self.symbols:
            for j in self.symbols:
                if i != j:
                    yield (i, j)

    def _set_gains(self, particle, index):
        sym_i = particle.get_symbol(index)
        self.atom_symbol[index] = sym_i
        base = np.dot(self.coeffs,
                      self.classifier.compute_atom_feature_for_symbol(particle, index, sym_i))
        for sym_j in self.symbols:
            if sym_j == sym_i:
                continue
            f_j = self.classifier.compute_atom_feature_for_symbol(particle, index, sym_j)
            self.flip_gain[(sym_i, sym_j)][index] = float(np.dot(self.coeffs, f_j) - base)
            self.indices[(sym_i, sym_j)].add(index)

    def bind_particle(self, particle):
        pairs = list(self._pairs())
        self.flip_gain = {pair: {} for pair in pairs}
        self.indices = {
            pair: SortedKeyList(key=lambda idx, p=pair: self.flip_gain[p][idx])
            for pair in pairs
        }
        self.atom_symbol = {}
        for index in particle.get_indices():
            self._set_gains(particle, index)

    def guided_exchange(self, particle):
        best = None  # (gain, A, B)
        for i, j in itertools.combinations(self.symbols, 2):
            la, lb = self.indices[(i, j)], self.indices[(j, i)]
            if la and lb:
                a, b = la[0], lb[0]
                gain = self.flip_gain[(i, j)][a] + self.flip_gain[(j, i)][b]
                if best is None or gain < best[0]:
                    best = (gain, a, b)
        if best is None:
            raise ValueError(
                "ETOP guided exchange needs at least two species present in the particle.")
        gain, index_a, index_b = best
        particle.swap_symbols([(index_a, index_b)])
        return index_a, index_b, gain

    def basin_hop_step(self, particle):
        candidates = []  # (gain, A, B)
        for i, j in itertools.combinations(self.symbols, 2):
            la, lb = self.indices[(i, j)], self.indices[(j, i)]
            if la and lb:
                a, b = la[0], lb[0]
                candidates.append((self.flip_gain[(i, j)][a] + self.flip_gain[(j, i)][b], a, b))
        if not candidates:
            raise ValueError(
                "ETOP basin hop needs at least two species present in the particle.")
        uphill = [c for c in candidates if c[0] > 0]
        gain, index_a, index_b = (min(uphill, key=lambda c: c[0]) if uphill
                                  else max(candidates, key=lambda c: c[0]))
        particle.swap_symbols([(index_a, index_b)])
        return index_a, index_b

    def update(self, particle, indices, exchange_indices):
        # exchange_indices is accepted for interface parity with the binary
        # operators; ETOP refreshes every atom in `indices` from self.atom_symbol.
        for index in indices:
            sym_old = self.atom_symbol[index]
            for sym_j in self.symbols:
                if sym_j != sym_old:
                    self.indices[(sym_old, sym_j)].remove(index)
                    del self.flip_gain[(sym_old, sym_j)][index]
        for index in indices:
            self._set_gains(particle, index)

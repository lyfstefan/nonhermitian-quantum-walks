"""Site-major (up, down) amplitudes; all products act from right to left."""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eig


def rotation(model, theta):
    c, s = np.cos(theta), np.sin(theta)
    if model == "U1":
        return np.array([[c, 1j*s], [1j*s, c]])
    if model == "U2":
        return np.array([[c, -s], [s, c]], dtype=complex)
    raise ValueError("model must be U1 or U2")


def gain(model, gamma):
    if not np.isfinite(gamma) or gamma < 0 or (model == "U2" and gamma == 0):
        raise ValueError("U1 requires gamma >= 0; U2 requires gamma > 0")
    if model == "U1":
        return np.diag(np.exp([-gamma, gamma]))
    if model == "U2":
        return np.array([[gamma+1/gamma, gamma-1/gamma],
                         [gamma-1/gamma, gamma+1/gamma]]) / 2
    raise ValueError("model must be U1 or U2")


def bulk_discriminant(model, gamma, theta1, theta2, k):
    """Half trace d0(k), Appendix A Eqs. (A2), (A8); phase defect absent."""
    gain(model, gamma)  # validate the convention and domain
    a = np.cos(theta1)*np.cos(theta2)*np.cos(2*np.asarray(k))
    b = np.sin(theta1)*np.sin(theta2)
    return a - np.cosh(2*gamma)*b if model == "U1" else (gamma+1/gamma)/2*(a-b)


def bulk_regime(model, gamma, theta1, theta2, tolerance=1e-10):
    """Classify the entire continuous Brillouin zone, not a sampled grid.

    |d0| < 1 has real quasienergy. Report a separate critical boundary
    whenever an endpoint touches |d0|=1, including a transition into a
    completely broken band. This describes the defect-free bulk only.
    """
    ends = bulk_discriminant(model, gamma, theta1, theta2, [0, np.pi/2])
    lo, hi = float(min(ends)), float(max(ends))
    minimum = 0.0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
    maximum = max(abs(lo), abs(hi))
    if abs(maximum-1) <= tolerance or abs(minimum-1) <= tolerance:
        regime = "critical"
    elif maximum < 1:
        regime = "unbroken"
    elif minimum > 1:
        regime = "completely broken"
    else:
        regime = "partially broken"
    return {"regime": regime, "min_abs_d0": minimum, "max_abs_d0": maximum}


@dataclass(frozen=True)
class Walk:
    model: str
    gamma: float
    half_size: int = 80
    theta1: float = -4*np.pi/9
    theta2: float = 5*np.pi/9 + 0.1113
    phase: float = np.pi
    defect: int = 0

    def __post_init__(self):
        gain(self.model, self.gamma)
        if not isinstance(self.half_size, int) or self.half_size < 2:
            raise ValueError("half_size must be an integer >= 2")
        if not isinstance(self.defect, int) or abs(self.defect) > self.half_size:
            raise ValueError("defect must be an integer within the lattice")
        if not np.all(np.isfinite([self.theta1, self.theta2, self.phase])):
            raise ValueError("angles must be finite")

    @property
    def sites(self):
        return np.arange(-self.half_size, self.half_size+1)

    def _coin(self, state, matrix):
        return np.einsum("ab,nb...->na...", matrix, state)

    def _shift(self, state):
        # S = S_free D: the phase is at the departure site, at EACH shift.
        phased = state.copy()
        phased[self.half_size+self.defect] *= np.exp(1j*self.phase)
        direction = 1 if self.model == "U1" else -1
        result = np.empty_like(phased)
        result[:, 0] = np.roll(phased[:, 0], direction, axis=0)
        result[:, 1] = np.roll(phased[:, 1], -direction, axis=0)
        return result

    def step(self, state):
        """Apply the linear operator without normalization; optional batch axis."""
        state = np.asarray(state, dtype=complex)
        if state.shape[:2] != (len(self.sites), 2):
            raise ValueError("expected (2*half_size+1, 2, ...) amplitudes")
        g = gain(self.model, self.gamma)
        if self.model == "U1":
            state = self._coin(state, rotation("U1", self.theta1))
            state = self._coin(state, g)
            state = self._shift(state)
            state = self._coin(state, rotation("U1", self.theta2))
            state = self._coin(state, np.diag(1/np.diag(g)))
            return self._shift(state)
        r1, r2 = rotation("U2", self.theta1/2), rotation("U2", self.theta2/2)
        state = self._coin(state, r1)
        state = self._shift(state)
        state = self._coin(state, r2)
        state = self._coin(state, g)
        state = self._coin(state, r2)
        state = self._shift(state)
        return self._coin(state, r1)

    def matrix(self):
        dimension = 2*len(self.sites)
        basis = np.eye(dimension, dtype=complex).reshape(len(self.sites), 2, dimension)
        return self.step(basis).reshape(dimension, dimension)

    def initial_state(self):
        state = np.zeros((len(self.sites), 2), dtype=complex)
        state[self.half_size] = np.array([1, 1j])/np.sqrt(2)
        return state

    def evolve(self, steps=300, initial=None):
        """Normalize amplitudes after every FULL step; return P(T,n), T=0...steps."""
        if not isinstance(steps, int) or steps < 0:
            raise ValueError("steps must be a nonnegative integer")
        state = self.initial_state() if initial is None else np.array(initial, dtype=complex)
        if state.shape != (len(self.sites), 2) or not np.isfinite(state).all():
            raise ValueError("invalid initial state")
        norm = np.linalg.norm(state)
        if norm == 0:
            raise ValueError("initial state must be nonzero")
        state /= norm
        probability = [np.sum(abs(state)**2, axis=1)]
        for _ in range(steps):
            state = self.step(state)
            norm = np.linalg.norm(state)
            if norm == 0 or not np.isfinite(norm):
                raise FloatingPointError("nonfinite or zero norm during evolution")
            state /= norm
            probability.append(np.sum(abs(state)**2, axis=1))
        return np.asarray(probability)

    def eigensystem(self):
        """Right eigenvectors, ordinary Euclidean normalization (not biorthogonal)."""
        matrix = self.matrix()
        values, vectors = eig(matrix)
        vectors /= np.linalg.norm(vectors, axis=0)
        profiles = np.sum(abs(vectors.reshape(len(self.sites), 2, -1))**2, axis=1)
        residuals = np.linalg.norm(matrix@vectors-vectors*values, axis=0)
        return values, profiles, residuals

    def localized_modes(self, radius=8, minimum_mass=0.95):
        """Finite-ring localization diagnostic, with explicit selection criteria.

        Not a proof of localization on the infinite line. Radius is in physical
        site coordinates; boundary-convergence checks are supplied separately.
        """
        values, profiles, residuals = self.eigensystem()
        mass = profiles[abs(self.sites-self.defect) <= radius].sum(axis=0)
        selected = np.flatnonzero(mass > minimum_mass)
        return values[selected], profiles[:, selected], residuals[selected]

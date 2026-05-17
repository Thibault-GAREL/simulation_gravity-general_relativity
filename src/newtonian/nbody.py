"""Simulation gravitationnelle newtonienne (2-corps et N-corps) en 2D.

Force entre deux masses i et j :
    F_ij = -G m_i m_j / |r_ij|² · r̂_ij

On intègre avec scipy.integrate.solve_ivp en méthode DOP853 (Runge-Kutta
d'ordre 8) qui préserve bien l'énergie et le moment angulaire sur des
durées longues.
"""

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy.integrate import solve_ivp

from ..config import G


@dataclass
class OrbitResult:
    """Résultat 2D d'une simulation à un corps en orbite.

    t : temps (s)
    x, y, vx, vy : trajectoire et vitesse (m, m/s)
    t_perihelia, phi_perihelia : instants et angles polaires (event-detected)
    des passages au périhélie. Permet de mesurer une éventuelle précession.
    """
    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    t_perihelia: np.ndarray = field(default_factory=lambda: np.array([]))
    phi_perihelia: np.ndarray = field(default_factory=lambda: np.array([]))

    @property
    def r(self) -> np.ndarray:
        return np.hypot(self.x, self.y)


# ---------------------------------------------------------------------------
# Cas 2-corps : masse centrale fixée à l'origine, corps test en orbite
# ---------------------------------------------------------------------------
def simulate_two_body(
    M_central: float,
    r0: float,
    v0: float,
    t_max: float,
    n_steps: int = 5000,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> OrbitResult:
    """Orbite képlérienne d'un corps test autour d'une masse centrale.

    Conditions initiales : périhélie en (r0, 0), vitesse (0, v0).
    Pour une orbite elliptique au périhélie : v0 = sqrt(GM(1+e)/(a(1-e))).

    Le périhélie est détecté via event sur dr/dt = (x·vx + y·vy)/r = 0
    dans le sens montant (dr/dt passe de - à +), avec précision machine.
    """
    mu = G * M_central

    def rhs(_t, state):
        x, y, vx, vy = state
        r3 = (x * x + y * y) ** 1.5
        ax = -mu * x / r3
        ay = -mu * y / r3
        return [vx, vy, ax, ay]

    # Événement : passage au périhélie ⟺ dr/dt = 0 dans le sens montant
    def perihelion_event(_t, state):
        x, y, vx, vy = state
        return x * vx + y * vy
    perihelion_event.direction = 1
    perihelion_event.terminal = False

    state0 = [r0, 0.0, 0.0, v0]
    t_eval = np.linspace(0.0, t_max, n_steps)

    sol = solve_ivp(
        rhs, (0.0, t_max), state0,
        t_eval=t_eval,
        method="DOP853",
        rtol=rtol, atol=atol,
        events=perihelion_event,
    )

    x, y, vx, vy = sol.y
    # Angles cumulés (déroulés) aux instants des périhélies
    t_peri = sol.t_events[0]
    if len(t_peri) > 0:
        xy_peri = sol.y_events[0][:, :2]
        # Angle déroulé : on suit la trajectoire pour ajouter les tours
        phi_traj = np.unwrap(np.arctan2(y, x))
        phi_peri = np.interp(t_peri, sol.t, phi_traj)
    else:
        phi_peri = np.array([])

    return OrbitResult(
        t=sol.t, x=x, y=y, vx=vx, vy=vy,
        t_perihelia=t_peri,
        phi_perihelia=phi_peri,
    )


# ---------------------------------------------------------------------------
# Cas N-corps général (toutes les masses sont libres)
# ---------------------------------------------------------------------------
def simulate_n_body(
    masses: Sequence[float],
    positions: Sequence[Sequence[float]],
    velocities: Sequence[Sequence[float]],
    t_max: float,
    n_steps: int = 5000,
    rtol: float = 1e-10,
    atol: float = 1e-12,
):
    """N-corps 2D avec interactions mutuelles.

    Returns t (T,), pos (T, N, 2), vel (T, N, 2).
    """
    masses = np.asarray(masses, dtype=float)
    pos0 = np.asarray(positions, dtype=float)
    vel0 = np.asarray(velocities, dtype=float)
    n = len(masses)

    def rhs(_t, state):
        p = state[: 2 * n].reshape(n, 2)
        v = state[2 * n :].reshape(n, 2)
        a = np.zeros_like(p)
        for i in range(n):
            d = p - p[i]
            dist2 = np.einsum("ij,ij->i", d, d)
            dist2[i] = 1.0
            inv_r3 = dist2 ** (-1.5)
            inv_r3[i] = 0.0
            a[i] = G * np.sum((masses * inv_r3)[:, None] * d, axis=0)
        return np.concatenate([v.ravel(), a.ravel()])

    state0 = np.concatenate([pos0.ravel(), vel0.ravel()])
    t_eval = np.linspace(0.0, t_max, n_steps)

    sol = solve_ivp(
        rhs, (0.0, t_max), state0,
        t_eval=t_eval, method="DOP853",
        rtol=rtol, atol=atol,
    )

    p = sol.y[: 2 * n].T.reshape(-1, n, 2)
    v = sol.y[2 * n :].T.reshape(-1, n, 2)
    return sol.t, p, v


# ---------------------------------------------------------------------------
# Diagnostics : énergie et moment angulaire
# ---------------------------------------------------------------------------
def specific_energy(result: OrbitResult, M_central: float) -> np.ndarray:
    """Énergie spécifique E/m = v²/2 - GM/r (J/kg). Doit être ~constante."""
    v2 = result.vx**2 + result.vy**2
    return 0.5 * v2 - G * M_central / result.r


def specific_angular_momentum(result: OrbitResult) -> np.ndarray:
    """Moment angulaire spécifique L/m = x·vy - y·vx (m²/s)."""
    return result.x * result.vy - result.y * result.vx

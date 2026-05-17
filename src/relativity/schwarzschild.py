"""Géodésiques dans la métrique de Schwarzschild (plan équatorial).

Métrique (coordonnées (t, r, θ, φ), θ=π/2) :

    ds² = -(1 - rs/r) c² dt² + (1 - rs/r)^-1 dr² + r² dφ²
    avec rs = 2GM/c²

On utilise l'**équation de Binet-Einstein** qui décrit la trajectoire spatiale
u(φ) avec u = 1/r :

  - Particule massive :   d²u/dφ² + u = GM/h² + (3GM/c²) u²
  - Photon (lumière)  :   d²u/dφ² + u =          (3GM/c²) u²

où h = L/m = moment angulaire spécifique (m²/s), constante du mouvement.

Avantages de cette formulation :
- Élimine le temps → trajectoire purement géométrique.
- Le terme `3GM·u²/c²` est LA correction relativiste pure. À c → ∞ on retrouve
  l'orbite képlérienne newtonienne (terme nul).
- Intégrable directement par scipy.solve_ivp.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ..config import G, C


def schwarzschild_radius(mass: float) -> float:
    """Rayon de Schwarzschild rs = 2GM/c² (m)."""
    return 2.0 * G * mass / C**2


def isco_radius(mass: float) -> float:
    """Innermost Stable Circular Orbit (ISCO) Schwarzschild : r = 3·rs = 6GM/c²."""
    return 3.0 * schwarzschild_radius(mass)


@dataclass
class GeodesicResult:
    """Trajectoire spatiale dans le plan équatorial.

    phi : angle polaire échantillonné (rad)
    r, x, y : coordonnées correspondantes (m)
    phi_perihelia : angles précis (event-detected) des passages au périhélie (rad)
    """
    phi: np.ndarray
    r: np.ndarray
    x: np.ndarray
    y: np.ndarray
    phi_perihelia: np.ndarray = field(default_factory=lambda: np.array([]))


# ---------------------------------------------------------------------------
# Orbite d'une particule massive (Mercure, étoile autour d'un trou noir)
# ---------------------------------------------------------------------------
def simulate_orbit_schwarzschild(
    M_central: float,
    r0: float,
    v0: float,
    n_orbits: float = 1.0,
    n_points: int = 5000,
    rtol: float = 1e-11,
    atol: float = 1e-13,
) -> GeodesicResult:
    """Orbite Schwarzschild d'une particule test (plan équatorial).

    Conditions initiales : périhélie en (r0, φ=0), vitesse tangentielle v0.
    On en déduit h = r0·v0 (moment angulaire spécifique conservé).

    Retourne aussi `phi_perihelia` : angles cumulés des périhélies, détectés
    via event-detection sur du/dφ = 0 (transitions montant→descendant), avec
    une précision largement supérieure à 1 pas d'échantillonnage.

    Parameters
    ----------
    M_central : masse centrale (kg)
    r0        : rayon initial = périhélie (m). Doit être >> rs.
    v0        : vitesse tangentielle initiale (m/s)
    n_orbits  : étendue angulaire totale (Δφ = 2π · n_orbits)
    n_points  : nombre d'échantillons pour l'affichage de la trajectoire
    """
    mu = G * M_central
    h = r0 * v0                       # moment angulaire spécifique

    # Non-dimensionalisation : U = u · r0, sans dimension, U(0) = 1.
    # Équation : d²U/dφ² + U = A + B·U²
    A = mu * r0 / h**2                       # ~1/(1+e) pour une orbite elliptique
    B = 3.0 * mu / (C**2 * r0)               # paramètre relativiste (≪ 1 hors trou noir)

    def rhs(_phi, state):
        U, dU = state
        d2U = A + B * U * U - U
        return [dU, d2U]

    # Événement : dU/dφ = 0 dans le sens descendant → max(U) → min(r) → périhélie
    def perihelion_event(_phi, state):
        return state[1]
    perihelion_event.direction = -1
    perihelion_event.terminal = False

    state0 = [1.0, 0.0]
    phi_max = 2.0 * np.pi * n_orbits
    phi_eval = np.linspace(0.0, phi_max, n_points)

    sol = solve_ivp(
        rhs, (0.0, phi_max), state0,
        t_eval=phi_eval,
        method="DOP853",
        rtol=rtol, atol=atol,
        events=perihelion_event,
    )

    u = sol.y[0] / r0                        # retour à u (1/m)
    # Tronquer si la trajectoire plonge sous rs (cas trou noir)
    rs = schwarzschild_radius(M_central)
    u_max = 1.0 / (1.01 * rs)
    if np.any(u > u_max):
        cut = int(np.argmax(u > u_max))
        sol_phi = sol.t[:cut]
        u = u[:cut]
    else:
        sol_phi = sol.t

    r = 1.0 / u
    x = r * np.cos(sol_phi)
    y = r * np.sin(sol_phi)

    return GeodesicResult(
        phi=sol_phi, r=r, x=x, y=y,
        phi_perihelia=np.asarray(sol.t_events[0]),
    )


# ---------------------------------------------------------------------------
# Géodésique nulle (photon) — pour la déviation de la lumière
# ---------------------------------------------------------------------------
def simulate_photon_schwarzschild(
    M_central: float,
    r_min: float,
    phi_max: float = np.pi,
    n_points: int = 5000,
    rtol: float = 1e-11,
    atol: float = 1e-13,
) -> GeodesicResult:
    """Trajectoire d'un photon dans Schwarzschild.

    Le photon est tangent au cercle r = r_min à φ = 0, on intègre des deux
    côtés. Pour b grand devant rs, l'angle de déviation théorique est
    α ≈ 4GM/(b·c²) (mesuré en 1919, ~1.75" pour b = R_sun).
    """
    mu = G * M_central
    rs_coef = 3.0 * mu / C**2

    def rhs(_phi, state):
        u, du = state
        d2u = rs_coef * u**2 - u
        return [du, d2u]

    state0 = [1.0 / r_min, 0.0]

    sol_fwd = solve_ivp(
        rhs, (0.0, phi_max), state0,
        t_eval=np.linspace(0.0, phi_max, n_points),
        method="DOP853", rtol=rtol, atol=atol,
    )
    sol_bwd = solve_ivp(
        rhs, (0.0, -phi_max), state0,
        t_eval=np.linspace(0.0, -phi_max, n_points),
        method="DOP853", rtol=rtol, atol=atol,
    )

    phi = np.concatenate([sol_bwd.t[::-1], sol_fwd.t[1:]])
    u = np.concatenate([sol_bwd.y[0][::-1], sol_fwd.y[0][1:]])
    valid = u > 0
    phi, u = phi[valid], u[valid]
    r = 1.0 / u
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return GeodesicResult(phi=phi, r=r, x=x, y=y)


def light_deflection_angle(M_central: float, b: float) -> float:
    """Angle de déviation théorique GR pour un photon : α = 4GM/(b·c²) (rad).

    Newton prédirait la moitié (2GM/bc²). Pour le Soleil et b = R_sun, on
    obtient α ≈ 1.75'' (mesuré lors de l'éclipse de 1919).
    """
    return 4.0 * G * M_central / (b * C**2)

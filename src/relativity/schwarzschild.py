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
from scipy.optimize import brentq

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
def critical_impact_parameter(mass: float) -> float:
    """Paramètre d'impact critique : b_crit = (3√3/2)·rs ≈ 2.598·rs.

    Pour b < b_crit, le photon est capturé par le trou noir (plonge sous
    l'horizon). Pour b > b_crit, il s'échappe avec une déviation finie.
    Au seuil exact b = b_crit, il orbite indéfiniment à r = 1.5·rs
    (photon sphère).
    """
    return 1.5 * np.sqrt(3.0) * schwarzschild_radius(mass)


def _periapse_from_impact(mass: float, b: float) -> float:
    """Rayon de plus proche approche r_min en fonction du paramètre d'impact b.

    Relation : b² = r_min² / (1 - rs/r_min)  ⇔  r_min³ - b²·r_min + b²·rs = 0.
    On prend la plus grande racine réelle > 1.5·rs (trajectoire évadée).
    """
    rs = schwarzschild_radius(mass)
    b_crit = critical_impact_parameter(mass)
    if b <= b_crit:
        raise ValueError(
            f"b = {b/rs:.3f}·rs ≤ b_crit = {b_crit/rs:.3f}·rs : photon capturé."
        )
    # Cherche r_min dans [1.5·rs, b] (la plus grande racine y est)
    f = lambda r: r**3 - b * b * r + b * b * rs
    return brentq(f, 1.5 * rs, b)


def simulate_photon_schwarzschild(
    M_central: float,
    b: float,
    phi_max: float = 4.0 * np.pi,
    n_points: int = 5000,
    rtol: float = 1e-11,
    atol: float = 1e-14,
) -> GeodesicResult:
    """Trajectoire d'un photon dans Schwarzschild (plan équatorial).

    Paramètres
    ----------
    M_central : masse de la lentille (kg)
    b         : paramètre d'impact asymptotique (m). Doit être > b_crit pour
                que le photon ne soit pas capturé.
    phi_max   : étendue angulaire maximale d'intégration de chaque côté (rad).

    Le photon est paramétré par son rayon de plus proche approche r_min(b)
    et a son périapse à φ = 0. On intègre l'équation de Binet pour photons
    des deux côtés (avec event detection à u = 0 pour stopper à l'infini).

    En faible champ (b ≫ rs) : α ≈ 4GM/(bc²). Près de b_crit, α diverge
    logarithmiquement.
    """
    r_min = _periapse_from_impact(M_central, b)
    B = 3.0 * G * M_central / (C**2 * r_min)

    def rhs(_phi, state):
        U, dU = state
        return [dU, B * U * U - U]

    def at_infinity(_phi, state):
        return state[0]
    at_infinity.terminal = True
    at_infinity.direction = -1

    state0 = [1.0, 0.0]    # U = u·r_min = 1 au périapse

    sol_fwd = solve_ivp(
        rhs, (0.0, phi_max), state0,
        t_eval=np.linspace(0.0, phi_max, n_points),
        method="DOP853", rtol=rtol, atol=atol,
        events=at_infinity,
    )
    sol_bwd = solve_ivp(
        rhs, (0.0, -phi_max), state0,
        t_eval=np.linspace(0.0, -phi_max, n_points),
        method="DOP853", rtol=rtol, atol=atol,
        events=at_infinity,
    )

    phi = np.concatenate([sol_bwd.t[::-1], sol_fwd.t[1:]])
    U = np.concatenate([sol_bwd.y[0][::-1], sol_fwd.y[0][1:]])
    valid = U > 0
    phi, U = phi[valid], U[valid]
    r = r_min / U
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return GeodesicResult(phi=phi, r=r, x=x, y=y)


def measure_deflection_angle(
    M_central: float,
    b: float,
    phi_max: float = 4.0 * np.pi,
    rtol: float = 1e-12,
    atol: float = 1e-15,
) -> float:
    """Angle de déviation α d'un photon (rad), mesuré avec event detection.

    α = (φ⁺ − φ⁻) − π où φ⁺ et φ⁻ sont les angles asymptotiques (u → 0).
    Lève ValueError si b ≤ b_crit (photon capturé).
    """
    r_min = _periapse_from_impact(M_central, b)
    B = 3.0 * G * M_central / (C**2 * r_min)

    def rhs(_phi, state):
        U, dU = state
        return [dU, B * U * U - U]

    def at_infinity(_phi, state):
        return state[0]
    at_infinity.terminal = True
    at_infinity.direction = -1

    state0 = [1.0, 0.0]

    sol_fwd = solve_ivp(
        rhs, (0.0, phi_max), state0,
        method="DOP853", rtol=rtol, atol=atol,
        events=at_infinity,
    )
    sol_bwd = solve_ivp(
        rhs, (0.0, -phi_max), state0,
        method="DOP853", rtol=rtol, atol=atol,
        events=at_infinity,
    )

    if len(sol_fwd.t_events[0]) == 0 or len(sol_bwd.t_events[0]) == 0:
        raise RuntimeError(
            "L'intégrateur n'a pas atteint u=0 dans la fenêtre angulaire — augmenter phi_max."
        )

    phi_plus = sol_fwd.t_events[0][0]
    phi_minus = sol_bwd.t_events[0][0]
    return float((phi_plus - phi_minus) - np.pi)


def light_deflection_angle(M_central: float, b: float) -> float:
    """Angle de déviation théorique GR (faible champ) : α = 4GM/(b·c²) (rad).

    Newton prédirait la moitié (2GM/bc²). Pour le Soleil et b = R_sun, on
    obtient α ≈ 1.75'' (mesuré lors de l'éclipse de 1919).
    """
    return 4.0 * G * M_central / (b * C**2)

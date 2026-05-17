"""Utilitaires : conversions d'unités, extraction d'éléments orbitaux."""

import numpy as np

from .config import AU, YEAR, CENTURY


# ---------------------------------------------------------------------------
# Conversions d'unités pratiques
# ---------------------------------------------------------------------------
def m_to_au(x):
    return np.asarray(x) / AU


def s_to_year(t):
    return np.asarray(t) / YEAR


def rad_to_arcsec(theta):
    """Radians vers secondes d'arc (1 rad = 206264.806... arcsec)."""
    return np.asarray(theta) * (180.0 / np.pi) * 3600.0


def rad_per_orbit_to_arcsec_per_century(rad_per_orbit, period_seconds):
    """Convertit une précession (rad/orbite) en arcsec/siècle.

    Utile pour comparer à la valeur observée de Mercure (~43"/siècle).
    """
    n_orbits_per_century = CENTURY / period_seconds
    return rad_to_arcsec(rad_per_orbit) * n_orbits_per_century


# ---------------------------------------------------------------------------
# Extraction d'éléments orbitaux à partir d'une trajectoire 2D
# ---------------------------------------------------------------------------
def precession_per_orbit(phi_perihelia):
    """Précession moyenne du périhélie par orbite (rad).

    Prend en entrée la liste des angles cumulés (déroulés) aux passages
    successifs au périhélie. Pour une orbite képlérienne fermée, l'écart
    entre périhélies est exactement 2π → précession = 0. Toute déviation
    représente la précession du périhélie.

    Parameters
    ----------
    phi_perihelia : array-like — angles polaires cumulés (rad) des
        passages au périhélie, mesurés avec event detection.
    """
    phis = np.asarray(phi_perihelia)
    if len(phis) < 2:
        return 0.0
    diffs = np.diff(phis)
    sign = np.sign(np.mean(diffs))   # sens du parcours orbital
    return float(np.mean(diffs - sign * 2.0 * np.pi))


def mercury_precession_theory(a, e, M_central):
    """Précession théorique GR par orbite : Δφ = 6π·GM/(c²·a·(1-e²)) (rad).

    Pour Mercure (a=5.79e10, e=0.206, M=M_sun) → ~5.02e-7 rad/orbite,
    soit ~43 arcsec/siècle.
    """
    from .config import G, C
    return 6.0 * np.pi * G * M_central / (C**2 * a * (1.0 - e**2))

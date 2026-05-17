"""Constantes physiques (SI) et paramètres de référence pour les simulations.

Toutes les valeurs sont en unités SI sauf mention contraire.
"""

# ---------------------------------------------------------------------------
# Constantes fondamentales (CODATA 2018, valeurs exactes ou recommandées)
# ---------------------------------------------------------------------------
G = 6.67430e-11           # Constante de gravitation universelle [m^3 kg^-1 s^-2]
C = 299_792_458.0         # Vitesse de la lumière dans le vide [m s^-1]

# ---------------------------------------------------------------------------
# Masses (kg)
# ---------------------------------------------------------------------------
M_SUN = 1.98892e30
M_MERCURY = 3.3011e23
M_VENUS = 4.8675e24
M_EARTH = 5.9722e24
M_MARS = 6.4171e23
M_JUPITER = 1.89813e27
M_SATURN = 5.6832e26
M_URANUS = 8.6811e25
M_NEPTUNE = 1.02409e26

# ---------------------------------------------------------------------------
# Distances (m)
# ---------------------------------------------------------------------------
AU = 1.495978707e11       # Unité astronomique
R_SUN = 6.957e8           # Rayon physique du Soleil

# Demi-grand axe (a) et excentricité (e) des orbites planétaires
ORBITAL_ELEMENTS = {
    # name:        (a [m],          e)
    "mercury": (5.7909e10,   0.20563),
    "venus":   (1.0821e11,   0.00678),
    "earth":   (1.4960e11,   0.01671),
    "mars":    (2.2794e11,   0.09340),
    "jupiter": (7.7857e11,   0.04839),
    "saturn":  (1.4335e12,   0.05415),
    "uranus":  (2.8725e12,   0.04717),
    "neptune": (4.4951e12,   0.00859),
}

# ---------------------------------------------------------------------------
# Temps (s)
# ---------------------------------------------------------------------------
DAY = 86_400.0
YEAR = 365.25 * DAY
CENTURY = 100.0 * YEAR

# ---------------------------------------------------------------------------
# Helpers dérivés
# ---------------------------------------------------------------------------
def schwarzschild_radius(mass: float) -> float:
    """Rayon de Schwarzschild rs = 2GM/c^2 (m)."""
    return 2.0 * G * mass / C**2


# Rayon de Schwarzschild du Soleil (~2953 m)
RS_SUN = schwarzschild_radius(M_SUN)

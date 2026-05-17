"""Helpers matplotlib pour comparer trajectoires Newton vs Schwarzschild."""

import matplotlib.pyplot as plt
import numpy as np


def plot_orbit_comparison(
    newton_xy,
    gr_xy,
    *,
    central_radius: float = 0.0,
    title: str = "Newton vs Relativité Générale",
    unit_scale: float = 1.0,
    unit_label: str = "m",
    ax=None,
):
    """Trace côte-à-côte la trajectoire newtonienne et celle de Schwarzschild.

    Parameters
    ----------
    newton_xy : tuple (x, y) — trajectoire newtonienne (m)
    gr_xy     : tuple (x, y) — trajectoire Schwarzschild (m)
    central_radius : rayon de la masse centrale à dessiner (m). 0 = invisible.
    unit_scale, unit_label : conversion d'affichage (ex: 1/AU, "UA").
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    xn, yn = newton_xy
    xg, yg = gr_xy

    ax.plot(np.asarray(xn) * unit_scale, np.asarray(yn) * unit_scale,
            label="Newton", color="#1f77b4", lw=1.4)
    ax.plot(np.asarray(xg) * unit_scale, np.asarray(yg) * unit_scale,
            label="Schwarzschild (GR)", color="#d62728", lw=1.4, alpha=0.85)

    # Masse centrale
    if central_radius > 0:
        circle = plt.Circle((0, 0), central_radius * unit_scale,
                            color="black", zorder=5)
        ax.add_patch(circle)
    else:
        ax.scatter([0], [0], marker="*", color="black", s=120, zorder=5,
                   label="Masse centrale")

    ax.set_aspect("equal")
    ax.set_xlabel(f"x [{unit_label}]")
    ax.set_ylabel(f"y [{unit_label}]")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    return ax


def plot_energy_conservation(t, energy, *, ax=None, title="Conservation d'énergie"):
    """Trace E(t)/E(0) - 1 pour vérifier la qualité de l'intégration."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 3.5))

    e0 = energy[0]
    drift = (energy - e0) / np.abs(e0)
    ax.plot(t, drift, color="#2ca02c", lw=1.2)
    ax.axhline(0, color="k", lw=0.5, alpha=0.5)
    ax.set_xlabel("temps (s)")
    ax.set_ylabel(r"$(E - E_0) / |E_0|$")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax

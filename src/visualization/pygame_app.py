"""Application Pygame pédagogique : Newton vs Schwarzschild en temps réel.

Affiche une orbite autour d'une masse centrale, en mode Newton, GR (Schwarzschild)
ou les deux côte-à-côte. L'utilisateur peut modifier la masse, le périapse et
l'excentricité au clavier pour explorer les différences.

Lancer avec :
    python -m src.visualization.pygame_app
ou bien le script run_pygame.py à la racine du projet.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pygame

from ..config import G, C, M_SUN
from ..newtonian import simulate_two_body
from ..relativity import (
    simulate_orbit_schwarzschild,
    schwarzschild_radius,
)


# ----------------------------------------------------------------------------
# Constantes de rendu
# ----------------------------------------------------------------------------
WIDTH, HEIGHT = 1200, 800
PANEL_W = 360                  # bande de droite avec les infos
PLOT_W = WIDTH - PANEL_W       # zone d'affichage de l'orbite
CENTER = (PLOT_W // 2, HEIGHT // 2)

BG       = (12, 14, 22)
PANEL_BG = (22, 26, 38)
GRID     = (40, 46, 60)
TEXT     = (220, 226, 240)
DIM      = (140, 148, 168)
ACCENT   = (255, 200, 80)
NEWTON_COL = (90, 170, 255)
GR_COL     = (240, 90, 110)
BH_COL   = (0, 0, 0)
ISCO_COL = (120, 120, 130)


# ----------------------------------------------------------------------------
# État de la simulation
# ----------------------------------------------------------------------------
@dataclass
class SimState:
    # Paramètres physiques (modifiables au clavier)
    log10_M_in_Msun: float = 6.0       # masse centrale (log10 en M_sun)
    r_peri_in_rs: float = 30.0         # périapse en rayons de Schwarzschild
    eccentricity: float = 0.5

    # Modes d'affichage
    mode: str = "BOTH"     # "NEWTON", "GR", ou "BOTH"
    paused: bool = False
    speed: int = 2         # pas par frame (1..10)
    show_help: bool = True
    n_orbits: int = 20     # durée de la simulation pré-calculée (avant boucle)

    # Données précomputées (remplies par recompute())
    newton_x: np.ndarray = field(default_factory=lambda: np.array([]))
    newton_y: np.ndarray = field(default_factory=lambda: np.array([]))
    gr_x: np.ndarray = field(default_factory=lambda: np.array([]))
    gr_y: np.ndarray = field(default_factory=lambda: np.array([]))
    cursor: int = 0
    rs: float = 1.0
    precession_deg: float = 0.0
    v_over_c: float = 0.0
    pixels_per_rs: float = 10.0

    # --------------------------------------------------------------------
    def mass(self) -> float:
        return (10.0 ** self.log10_M_in_Msun) * M_SUN

    def recompute(self, n_orbits: int | None = None, n_points: int | None = None) -> None:
        """Re-simule Newton et GR avec les paramètres courants.

        Si n_orbits/n_points est None, on prend les valeurs courantes de l'état
        (n_orbits = state.n_orbits, n_points = 600 × n_orbits ≈ 12000 par défaut).
        """
        if n_orbits is None:
            n_orbits = self.n_orbits
        if n_points is None:
            n_points = max(2000, 600 * n_orbits)
        M = self.mass()
        self.rs = schwarzschild_radius(M)
        r_peri = self.r_peri_in_rs * self.rs
        e = self.eccentricity
        a = r_peri / (1.0 - e)
        v_peri = float(np.sqrt(G * M * (1.0 + e) / (a * (1.0 - e))))
        T = 2.0 * np.pi * np.sqrt(a ** 3 / (G * M))
        self.v_over_c = v_peri / C

        # Newton
        res_n = simulate_two_body(M, r_peri, v_peri,
                                   t_max=n_orbits * T, n_steps=n_points)
        self.newton_x = res_n.x
        self.newton_y = res_n.y

        # GR
        res_gr = simulate_orbit_schwarzschild(M, r_peri, v_peri,
                                               n_orbits=n_orbits, n_points=n_points)
        self.gr_x = res_gr.x
        self.gr_y = res_gr.y

        if len(res_gr.phi_perihelia) >= 2:
            dphi = float(np.diff(res_gr.phi_perihelia)[0] - 2.0 * np.pi)
            self.precession_deg = np.degrees(dphi)
        else:
            self.precession_deg = 0.0

        # Échelle d'affichage : on cadre l'orbite Newton (la plus large des deux)
        r_max = max(float(np.hypot(self.newton_x, self.newton_y).max()),
                    float(np.hypot(self.gr_x, self.gr_y).max()))
        margin = 0.92
        plot_half = min(PLOT_W, HEIGHT) / 2.0
        self.pixels_per_rs = (plot_half * margin) / (r_max / self.rs)

        self.cursor = 0


# ----------------------------------------------------------------------------
# Fonctions de rendu
# ----------------------------------------------------------------------------
def world_to_screen(x: float, y: float, state: SimState) -> tuple[int, int]:
    px = CENTER[0] + (x / state.rs) * state.pixels_per_rs
    py = CENTER[1] - (y / state.rs) * state.pixels_per_rs    # y inversé (écran)
    return int(px), int(py)


def draw_trail(surface, xs, ys, color, state: SimState, end_idx: int) -> None:
    if end_idx < 2:
        return
    pts = [world_to_screen(xs[i], ys[i], state) for i in range(end_idx)]
    pygame.draw.lines(surface, color, False, pts, 2)


def draw_orbit_zone(surface, state: SimState, fonts) -> None:
    # Fond + grille
    pygame.draw.rect(surface, BG, (0, 0, PLOT_W, HEIGHT))
    for r_ring in (5, 10, 20, 50, 100):
        radius_px = int(r_ring * state.pixels_per_rs)
        if 20 < radius_px < min(PLOT_W, HEIGHT):
            pygame.draw.circle(surface, GRID, CENTER, radius_px, 1)

    # Horizon (rs) + ISCO (3 rs)
    rs_px = max(3, int(state.pixels_per_rs))
    isco_px = max(5, int(3 * state.pixels_per_rs))
    pygame.draw.circle(surface, BH_COL, CENTER, rs_px)
    pygame.draw.circle(surface, ISCO_COL, CENTER, isco_px, 1)

    # Trails
    cur = state.cursor
    if state.mode in ("NEWTON", "BOTH"):
        draw_trail(surface, state.newton_x, state.newton_y, NEWTON_COL, state, cur)
    if state.mode in ("GR", "BOTH"):
        draw_trail(surface, state.gr_x, state.gr_y, GR_COL, state, cur)

    # Particule courante (la plus avancée des deux modes affichés)
    if cur > 0:
        if state.mode in ("NEWTON", "BOTH"):
            pygame.draw.circle(surface, NEWTON_COL,
                               world_to_screen(state.newton_x[cur - 1], state.newton_y[cur - 1], state),
                               4)
        if state.mode in ("GR", "BOTH"):
            pygame.draw.circle(surface, GR_COL,
                               world_to_screen(state.gr_x[cur - 1], state.gr_y[cur - 1], state),
                               4)


def draw_panel(surface, state: SimState, fonts) -> None:
    x0 = PLOT_W
    pygame.draw.rect(surface, PANEL_BG, (x0, 0, PANEL_W, HEIGHT))

    def line(y, text, color=TEXT, font=fonts["body"]):
        surf = font.render(text, True, color)
        surface.blit(surf, (x0 + 18, y))

    y = 24
    line(y, "Gravité : Newton vs GR", color=ACCENT, font=fonts["title"]); y += 38

    # État physique
    line(y, "Paramètres", color=DIM, font=fonts["bold"]); y += 24
    line(y, f"  Masse centrale : 10^{state.log10_M_in_Msun:.1f} M⊙"); y += 22
    line(y, f"  r_périapse     : {state.r_peri_in_rs:.1f} r_s"); y += 22
    line(y, f"  Excentricité   : {state.eccentricity:.2f}"); y += 22
    line(y, f"  v/c (au péri)  : {state.v_over_c:.3f}"); y += 22
    line(y, f"  r_s            : {state.rs:.2e} m"); y += 30

    # Résultats
    line(y, "Mesures GR", color=DIM, font=fonts["bold"]); y += 24
    prec_str = f"{state.precession_deg:+.2f}°/orbite"
    color = GR_COL if abs(state.precession_deg) > 0.001 else TEXT
    line(y, f"  Précession     : {prec_str}", color=color); y += 22

    if state.r_peri_in_rs < 3.05:
        line(y, "  ⚠ sous l'ISCO", color=GR_COL); y += 22
    else:
        y += 22

    y += 12
    # Mode
    line(y, "Affichage", color=DIM, font=fonts["bold"]); y += 24
    mode_label = {"NEWTON": "Newton seul", "GR": "GR seul", "BOTH": "Newton + GR"}[state.mode]
    line(y, f"  Mode  : {mode_label}", color=ACCENT); y += 22
    line(y, f"  Pause : {'oui' if state.paused else 'non'}"); y += 22
    line(y, f"  Vitesse : x{state.speed}"); y += 22
    line(y, f"  Orbites : {state.n_orbits}"); y += 30

    # Aide
    if state.show_help:
        line(y, "Contrôles", color=DIM, font=fonts["bold"]); y += 24
        for txt in [
            "  N / G / B  : mode Newton/GR/Both",
            "  ↑ / ↓      : masse centrale",
            "  ← / →      : r_périapse",
            "  E / W      : excentricité ± 0.05",
            "  + / -      : vitesse anim.",
            "  P / O      : nb d'orbites (durée)",
            "  Espace     : pause",
            "  R          : reset trajectoire",
            "  H          : cacher cette aide",
            "  Échap      : quitter",
        ]:
            line(y, txt, color=DIM, font=fonts["small"]); y += 18

    # Légende couleurs en bas
    y_legend = HEIGHT - 60
    pygame.draw.line(surface, NEWTON_COL, (x0 + 18, y_legend), (x0 + 60, y_legend), 3)
    surface.blit(fonts["small"].render("Newton", True, TEXT), (x0 + 70, y_legend - 8))
    pygame.draw.line(surface, GR_COL, (x0 + 18, y_legend + 24), (x0 + 60, y_legend + 24), 3)
    surface.blit(fonts["small"].render("Schwarzschild (GR)", True, TEXT), (x0 + 70, y_legend + 16))


# ----------------------------------------------------------------------------
# Boucle principale
# ----------------------------------------------------------------------------
def run() -> None:
    pygame.init()
    pygame.display.set_caption("Newton vs Relativité Générale — Schwarzschild")
    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    clock = pygame.time.Clock()

    fonts = {
        "title": pygame.font.SysFont("Consolas", 22, bold=True),
        "bold":  pygame.font.SysFont("Consolas", 16, bold=True),
        "body":  pygame.font.SysFont("Consolas", 16),
        "small": pygame.font.SysFont("Consolas", 14),
    }

    state = SimState()
    state.recompute()

    running = True
    while running:
        # ---- Events ----
        needs_recompute = False
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                key = event.key
                if key == pygame.K_ESCAPE:
                    running = False
                elif key == pygame.K_SPACE:
                    state.paused = not state.paused
                elif key == pygame.K_h:
                    state.show_help = not state.show_help
                elif key == pygame.K_r:
                    state.cursor = 0
                elif key == pygame.K_n:
                    state.mode = "NEWTON"
                elif key == pygame.K_g:
                    state.mode = "GR"
                elif key == pygame.K_b:
                    state.mode = "BOTH"
                elif key in (pygame.K_PLUS, pygame.K_KP_PLUS, pygame.K_EQUALS):
                    state.speed = min(20, state.speed + 1)
                elif key in (pygame.K_MINUS, pygame.K_KP_MINUS):
                    state.speed = max(1, state.speed - 1)
                elif key == pygame.K_UP:
                    state.log10_M_in_Msun = min(10.0, state.log10_M_in_Msun + 0.5)
                    needs_recompute = True
                elif key == pygame.K_DOWN:
                    state.log10_M_in_Msun = max(0.0, state.log10_M_in_Msun - 0.5)
                    needs_recompute = True
                elif key == pygame.K_RIGHT:
                    state.r_peri_in_rs = min(500.0, state.r_peri_in_rs + 5.0)
                    needs_recompute = True
                elif key == pygame.K_LEFT:
                    state.r_peri_in_rs = max(3.05, state.r_peri_in_rs - 5.0)
                    needs_recompute = True
                elif key == pygame.K_e:
                    state.eccentricity = min(0.85, state.eccentricity + 0.05)
                    needs_recompute = True
                elif key == pygame.K_w:
                    state.eccentricity = max(0.0, state.eccentricity - 0.05)
                    needs_recompute = True
                elif key == pygame.K_p:
                    state.n_orbits = min(200, state.n_orbits + 10)
                    needs_recompute = True
                elif key == pygame.K_o:
                    state.n_orbits = max(2, state.n_orbits - 10)
                    needs_recompute = True

        if needs_recompute:
            state.recompute()

        # ---- Update ----
        if not state.paused:
            n_pts = max(len(state.newton_x), len(state.gr_x))
            state.cursor = min(n_pts, state.cursor + state.speed)
            if state.cursor >= n_pts:
                state.cursor = 0    # loop

        # ---- Render ----
        screen.fill(BG)
        draw_orbit_zone(screen, state, fonts)
        draw_panel(screen, state, fonts)
        pygame.display.flip()

        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    run()

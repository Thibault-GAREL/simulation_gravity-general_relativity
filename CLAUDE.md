# CLAUDE.md — Simulation Gravité Newton vs Relativité Générale

## Objectif du projet

Simuler et **comparer** la gravité selon deux théories :

1. **Mécanique newtonienne** — force attractive `F = -G·M·m/r²`.
2. **Relativité Générale d'Einstein** — métrique de **Schwarzschild** (solution exacte pour une masse sphérique statique), géodésiques calculées numériquement.

Le projet est **pédagogique** : il doit aider à *comprendre* les différences observables entre les deux théories, pas juste à les calculer.

---

## Cas physiques à simuler

Quatre scénarios, du plus faible au plus fort régime relativiste :

| # | Scénario | Effet GR attendu | Pourquoi c'est intéressant |
|---|----------|------------------|----------------------------|
| 1 | **Précession du périhélie de Mercure** | ~43"/siècle de précession en plus | Le test historique de la GR (1915). Effet faible mais mesurable. |
| 2 | **Orbite proche d'un trou noir (Schwarzschild)** | Précession énorme, ISCO à `r=6GM/c²`, plongée | Régime relativiste fort. Spectaculaire. |
| 3 | **Déviation de la lumière (lentille gravitationnelle)** | Newton prédit la moitié de l'angle observé | Confirmation de la GR (éclipse de 1919). Géodésiques nulles. |
| 4 | **Système solaire / orbites multiples** | Petites corrections relativistes par planète | Montre que Newton suffit en pratique sauf cas extrêmes. |

---

## Choix techniques

### Théorie

- **Approche GR** : métrique de Schwarzschild (statique, sphérique, à symétrie sphérique).
- **Équation maîtresse en GR (plan équatorial)** : équation de Binet-Einstein
  - Pour les particules massives : `d²u/dφ² + u = GM/h² + 3GM·u²/c²` où `u = 1/r`, `h` = moment angulaire spécifique.
  - Pour les photons : `d²u/dφ² + u = 3GM·u²/c²` (sans le terme newtonien).
- **Approche Newton** : intégration directe `F = -GMm·r̂/r²` (Runge-Kutta ou `scipy.solve_ivp`).

### Unités

- **Système SI** : mètres, secondes, kilogrammes.
- Constantes réalistes : `G = 6.674e-11`, `c = 299792458`, `M_sun = 1.989e30`, etc.
- Pour la stabilité numérique, possibilité d'utiliser des unités internes (UA, années, M_sun) puis reconversion à l'affichage.

### Visualisation — deux phases

**Phase 1 — Notebooks Jupyter avec widgets** (en cours)
- Un notebook par cas physique dans `notebooks/`.
- Cellules : théorie → code → simulation → comparaison Newton/GR.
- `ipywidgets` pour ajuster les paramètres (masse, distance, vitesse) à la volée.
- Plots `matplotlib` (trajectoires 2D, énergie/moment cinétique, écart Newton-GR).

**Phase 2 — Application Pygame** (en place pour le scénario orbital)
- `src/visualization/pygame_app.py` : visualisation temps réel, switch Newton/GR (`N`/`G`/`B`).
- Contrôles : masse (`↑↓`), périapse (`←→`), excentricité (`E`/`W`), vitesse (`+`/`-`).
- À l'état actuel : un seul scénario (orbite autour d'une masse centrale). Les cas
  *déviation de la lumière* et *N-corps système solaire* n'y sont pas (ils restent
  dans les notebooks 03 et 04). Une éventuelle sandbox plus large viendra plus tard
  si besoin.

---

## Architecture du projet

Structure adaptée du skill `thibault-ia-init`, mais **épurée** car ce n'est pas un projet ML — pas de `data/` ni `outputs/models/`.

```
simulation_gravity-general_relativity/
├── src/
│   ├── config.py               # Constantes physiques + paramètres globaux
│   ├── utils.py                # Helpers (conversions d'unités, etc.)
│   ├── newtonian/              # Intégrateurs Newton (N-corps, 2-corps)
│   ├── relativity/             # Géodésiques Schwarzschild (massif + photon)
│   └── visualization/          # Helpers matplotlib (+ pygame en phase 2)
├── notebooks/                  # 1 notebook par cas physique
├── tests/                      # Tests unitaires (conservation énergie, etc.)
├── assets/                     # Images/GIFs pour le README
├── requirements.txt
├── LICENSE                     # MIT (Thibault GAREL, 2026)
└── README.md
```

---

## Environnement

- **venv** : `pytorch_cuda_env` (il contient déjà `jupyter`, `ipywidgets`, `scipy`, `matplotlib`, `pygame`).
- Activation : `& c:\0-Code_py_temp\pytorch_cuda_env\Scripts\Activate.ps1`
- Pas besoin de PyTorch/CUDA pour ce projet (calcul purement numérique CPU).

---

## Dépôt GitHub

https://github.com/Thibault-GAREL/simulation_gravity-general_relativity

---

## Règles de travail (pour Claude)

1. **Clarté pédagogique > performance** — ce projet doit s'expliquer. Variables nommées, formules dans les docstrings, commentaires pour le *pourquoi* physique (pas le *quoi* trivial).
2. **Pas d'over-engineering** — un script qui marche et qu'on comprend > une archi parfaite.
3. **Vérifier la physique** : conservation de l'énergie et du moment cinétique sur les cas newtoniens, valeurs analytiques connues pour la GR (43"/siècle Mercure, 1.75" déviation à la surface du Soleil).
4. **Une étape à la fois** : implémenter et tester un cas avant de passer au suivant.
5. **Avant chaque changement structurel** : git commit + push (règle utilisateur).

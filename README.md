# simulation_gravity-general_relativity

Simulation comparative de la gravité selon **Newton** et selon la **Relativité Générale d'Einstein** (métrique de Schwarzschild), en Python.

Projet pédagogique : on implémente les deux théories côte à côte et on visualise les différences observables sur des cas classiques (précession de Mercure, orbite autour d'un trou noir, déviation de la lumière, système solaire).

## État du projet

- [x] Squelette du projet, modules `newtonian` & `relativity`
- [x] **Notebook 01 — Précession du périhélie de Mercure** ([notebooks/01_mercury_precession.ipynb](notebooks/01_mercury_precession.ipynb))
  - Résultat simulation GR : **42.998''/siècle** (théorie : 42.997'', observation : ~43'')
- [ ] Notebook 02 — Orbite proche d'un trou noir (Schwarzschild)
- [ ] Notebook 03 — Déviation de la lumière (lentille gravitationnelle)
- [ ] Notebook 04 — Système solaire (N-corps)
- [ ] Application Pygame interactive (phase 2)

## Installation & lancement

L'environnement Python utilisé est `pytorch_cuda_env` (contient déjà `numpy`, `scipy`, `matplotlib`, `jupyter`, `ipywidgets`, `pygame`).

```powershell
# Activer l'env
& c:\0-Code_py_temp\pytorch_cuda_env\Scripts\Activate.ps1

# Lancer Jupyter
cd D:\Loisir\Code_python\simulation_gravity-general_relativity
jupyter lab
```

Puis ouvrir `notebooks/01_mercury_precession.ipynb` et sélectionner le bon kernel.

## Structure

```
src/
├── config.py             # constantes physiques (G, c, M_sun, ...)
├── utils.py              # mesure de précession, conversions d'unités
├── newtonian/nbody.py    # intégrateur Newton (2-corps, N-corps)
├── relativity/           # géodésiques Schwarzschild
│   └── schwarzschild.py    (équation de Binet-Einstein)
└── visualization/        # helpers matplotlib (+ pygame en phase 2)

notebooks/                # 1 notebook par cas physique
```

Voir [CLAUDE.md](CLAUDE.md) pour les détails techniques et la physique sous-jacente.

## Licence

MIT — Thibault GAREL, 2026.

# Ajustement Fourier d'une RR Lyrae

Ce projet contient deux scripts Python pour travailler sur une courbe de lumiere
d'etoile variable de type RR Lyrae, a partir de fichiers CSV photometriques.

Les donnees attendues sont dans le dossier `data/`, avec un fichier CSV par nuit.
Chaque fichier doit avoir le separateur `;` et l'entete suivante :

```csv
JD;Magnitude;ErreurMagnitude
```

Le projet s'utilise avec `uv`.

## Installation / execution

Les dependances sont declarees dans `pyproject.toml` :

- `numpy`
- `matplotlib`
- `astropy`

Lancer les scripts avec :

```powershell
uv run python .\nom_du_script.py [options]
```

Les images et autres sorties generees sont ecrites par defaut dans le dossier
`output/`. Ce dossier est ignore par Git.

## Script 1 : `fourier_rrlyrae.py`

Ce script ajuste une serie de Fourier sur une courbe de lumiere lorsque la
periode est connue.

Le modele ajuste est :

```text
m(t) = A0 + somme_k ak cos(2*pi*k*(t-t0)/P) + bk sin(2*pi*k*(t-t0)/P)
```

L'ajustement est fait par moindres carres ponderes avec les incertitudes sur la
magnitude :

```text
sum(((m_obs - m_modele) / sigma_m)^2)
```

Le script affiche :

- la periode utilisee ;
- l'epoque `t0` ;
- l'ordre Fourier ;
- la magnitude moyenne `A0` ;
- le `chi2 reduit` ;
- l'amplitude pic-a-pic de la variation de magnitude estimee sur le modele ;
- le temps de montée en pourcentage du cycle et en jours ;
- les coefficients `ak`, `bk` ;
- la forme amplitude-phase equivalente `Ak`, `phik`.
- les parametres Fourier utiles pour les RR Lyrae : `R21`, `R31`, `phi21`, `phi31`.

Il genere aussi un graphique de la courbe repliee en phase, avec les barres
d'erreur et le fit Fourier superpose. Les points sont distingues par nuit avec
des couleurs et des symboles differents, et une legende indique la nuit
correspondante.

Les parametres RR Lyrae sont calcules ainsi :

```text
R21 = A2 / A1
R31 = A3 / A1
phi21 = phi2 - 2*phi1
phi31 = phi3 - 3*phi1
```

Les phases composees sont ramenees dans l'intervalle `[0, 2*pi[`.

Par convention pour les RR Lyrae, l'epoque `T0` est prise au maximum de
lumiere, c'est-a-dire au minimum de magnitude. Si `--epoch` est fourni, cette
valeur est toujours prioritaire. Sinon, `T0` est estime automatiquement avec
`--epoch-method`.

Deux methodes automatiques sont disponibles :

- `local-poly` : ajuste une parabole locale ponderee sur les mesures proches du
  maximum de lumiere. C'est la methode par defaut.
- `model` : prend le maximum de lumiere du modele Fourier ajuste.

L'amplitude est calculee sur le modele ajuste :

```text
amplitude = magnitude_max - magnitude_min
```

Le temps de montée correspond ici a la montée en luminosité, donc au trajet
entre :

```text
magnitude maximale -> magnitude minimale
luminosite minimale -> luminosite maximale
```

Elle est affichee comme fraction du cycle, en pourcentage, et comme duree en
jours :

```text
duree_montee = fraction_montee * periode
```

### Exemple avec la periode connue

```powershell
uv run python .\fourier_rrlyrae.py 0.293162 --order 4 --plot .\output\folded_light_curve.png
```

### Parametres de ligne de commande

```text
period
```

Periode connue, en jours. C'est le seul argument obligatoire.

```text
--order {3..10}
```

Ordre de la serie de Fourier. Valeur par defaut : `4`.

```text
--epoch EPOCH
```

Epoque `t0` en JD. Si fourni, cette valeur est prioritaire.

```text
--epoch-method {local-poly,model}
```

Methode automatique pour estimer `t0` si `--epoch` est absent. Valeur par
defaut : `local-poly`.

```text
--data DATA
```

Dossier contenant les fichiers CSV. Valeur par defaut : `data`.

```text
--pattern PATTERN
```

Motif des fichiers a lire. Valeur par defaut : `*.csv`.

```text
--plot PLOT
```

Chemin du graphique PNG a produire. Valeur par defaut :
`output/folded_light_curve.png`.

```text
--show
```

Affiche aussi le graphique dans une fenetre Matplotlib.

```text
--label-mode {date,file}
```

Choix des etiquettes de nuit dans la legende. `date` utilise le debut du nom de
fichier avant le premier `_`, par exemple `2026-06-18`. `file` utilise le nom
complet du fichier sans extension. Valeur par defaut : `date`.

```text
--single-style
```

Trace tous les points avec la meme couleur et le meme symbole, au lieu de
distinguer les nuits.

## Script 2 : `find_period.py`

Ce script estime la periode a partir des donnees.

Il teste une grille de periodes candidates. Pour chaque periode, il ajuste une
serie de Fourier ponderee, puis retient la periode qui minimise le `chi2 reduit`.

Par defaut, la recherche est faite dans une plage typique RRc :

```text
0.2 j <= P <= 0.5 j
```

Le script genere deux graphiques :

- `output/period_search.png` : evolution du `chi2 reduit` selon la periode candidate ;
- `output/best_period_folded_light_curve.png` : courbe repliee avec la meilleure periode.

### Exemple simple

```powershell
uv run python .\find_period.py --order 4
```

### Exemple avec plage personnalisee

```powershell
uv run python .\find_period.py --min-period 0.2 --max-period 0.5 --order 4
```

### Exemple avec plus de resolution

```powershell
uv run python .\find_period.py --order 4 --samples 20000 --refine-samples 10000
```

### Parametres de ligne de commande

```text
--min-period MIN_PERIOD
```

Periode minimale testee, en jours. Valeur par defaut : `0.2`.

```text
--max-period MAX_PERIOD
```

Periode maximale testee, en jours. Valeur par defaut : `0.5`.

```text
--samples SAMPLES
```

Nombre de periodes testees dans la grille grossiere. Valeur par defaut :
`5000`.

```text
--refine-samples REFINE_SAMPLES
```

Nombre de periodes testees autour du meilleur resultat de la grille grossiere.
Valeur par defaut : `3000`.

```text
--order {1..10}
```

Ordre Fourier utilise pour la recherche de periode. Valeur par defaut : `4`.

```text
--epoch EPOCH
```

Epoque `t0` en JD. Si fourni, cette valeur est prioritaire.

```text
--epoch-method {local-poly,model}
```

Methode automatique pour estimer `t0` du fit final si `--epoch` est absent.
Valeur par defaut : `local-poly`.

```text
--data DATA
```

Dossier contenant les fichiers CSV. Valeur par defaut : `data`.

```text
--pattern PATTERN
```

Motif des fichiers a lire. Valeur par defaut : `*.csv`.

```text
--top TOP
```

Nombre de meilleures periodes candidates affichees. Valeur par defaut : `10`.

```text
--periodogram-plot PERIODOGRAM_PLOT
```

Chemin du graphique de recherche de periode. Valeur par defaut :
`output/period_search.png`.

```text
--folded-plot FOLDED_PLOT
```

Chemin du graphique de courbe repliee avec la meilleure periode. Valeur par
defaut : `output/best_period_folded_light_curve.png`.

```text
--show
```

Affiche aussi les graphiques dans une fenetre Matplotlib.

```text
--label-mode {date,file}
```

Choix des etiquettes de nuit dans la legende de la courbe repliee. Valeur par
defaut : `date`.

```text
--single-style
```

Trace tous les points de la courbe repliee avec la meme couleur et le meme
symbole.

## Script 3 : `find_period_gls.py`

Ce script estime la periode avec la methode GLS, Generalized Lomb-Scargle, en
utilisant `astropy.timeseries.LombScargle`.

La methode GLS ajuste pour chaque frequence candidate un modele sinusoidal avec
moyenne flottante, en tenant compte des incertitudes de magnitude. Le script
retient la periode qui maximise la puissance GLS.

Comme pour les autres scripts, il lit les fichiers CSV du dossier `data/`.

Il genere deux graphiques :

- `output/gls_periodogram.png` : periodogramme GLS ;
- `output/gls_best_period_folded_light_curve.png` : courbe repliee avec la meilleure periode GLS.

La courbe repliee finale est affichee avec un fit Fourier, dont l'ordre est
choisi avec `--fit-order`.

### Exemple simple

```powershell
uv run python .\find_period_gls.py --fit-order 4
```

### Exemple avec sorties explicites

```powershell
uv run python .\find_period_gls.py --fit-order 4 --periodogram-plot .\output\gls_periodogram.png --folded-plot .\output\gls_best_period_folded_light_curve.png
```

### Parametres de ligne de commande

```text
--min-period MIN_PERIOD
```

Periode minimale testee, en jours. Valeur par defaut : `0.2`.

```text
--max-period MAX_PERIOD
```

Periode maximale testee, en jours. Valeur par defaut : `0.5`.

```text
--samples SAMPLES
```

Nombre de frequences testees dans la grille grossiere. Valeur par defaut :
`10000`.

```text
--refine-samples REFINE_SAMPLES
```

Nombre de frequences testees autour du meilleur resultat de la grille grossiere.
Valeur par defaut : `5000`.

```text
--fit-order {3..10}
```

Ordre Fourier utilise pour tracer la courbe repliee finale. Valeur par defaut :
`4`.

```text
--epoch EPOCH
```

Epoque `t0` du fit Fourier final en JD. Si fourni, cette valeur est prioritaire.

```text
--epoch-method {local-poly,model}
```

Methode automatique pour estimer `t0` du fit final si `--epoch` est absent.
Valeur par defaut : `local-poly`.

```text
--data DATA
```

Dossier contenant les fichiers CSV. Valeur par defaut : `data`.

```text
--pattern PATTERN
```

Motif des fichiers a lire. Valeur par defaut : `*.csv`.

```text
--top TOP
```

Nombre de meilleures periodes candidates affichees. Valeur par defaut : `10`.

```text
--periodogram-plot PERIODOGRAM_PLOT
```

Chemin du graphique GLS. Valeur par defaut : `output/gls_periodogram.png`.

```text
--folded-plot FOLDED_PLOT
```

Chemin du graphique de courbe repliee avec la meilleure periode GLS. Valeur par
defaut : `output/gls_best_period_folded_light_curve.png`.

```text
--show
```

Affiche aussi les graphiques dans une fenetre Matplotlib.

```text
--label-mode {date,file}
```

Choix des etiquettes de nuit dans la legende de la courbe repliee. Valeur par
defaut : `date`.

```text
--single-style
```

Trace tous les points de la courbe repliee avec la meme couleur et le meme
symbole.

## Notes

Les durees sont toujours affichees en jours. Quand la valeur est inferieure a
`1` jour, une conversion en heures est ajoutee pour faciliter la lecture.

Pour les magnitudes, l'axe vertical des graphiques est inverse : les magnitudes
plus faibles apparaissent plus haut, comme c'est l'usage en photometrie.

La recherche de periode par grille peut presenter des minima secondaires ou des
alias. Il est donc utile de comparer le resultat avec la periode attendue, la
couverture temporelle des observations et l'aspect de la courbe repliee.

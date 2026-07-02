# Ajustement Fourier d'une RR Lyrae

Ce projet contient trois scripts Python pour travailler sur une courbe de lumière
d'étoile variable de type RR Lyrae, à partir de fichiers CSV photométriques.

Les données attendues sont dans le dossier `data/`, avec un fichier CSV par nuit.
Chaque fichier doit avoir le séparateur `;` et l'entête suivante :

```csv
JD;Magnitude;ErreurMagnitude
```

Le projet s'utilise avec `uv`.

## Installation / exécution

Les dépendances sont déclarées dans `pyproject.toml` :

- `numpy`
- `matplotlib`
- `astropy`

Lancer les scripts avec :

```powershell
uv run python .\nom_du_script.py [options]
```

Les images et autres sorties générées sont écrites par défaut dans le dossier
`output/`. Ce dossier est ignoré par Git.

## Script 1 : `fourier_rrlyrae.py`

Ce script ajuste une série de Fourier sur une courbe de lumière lorsque la
période est connue.

Le modèle ajusté est :

```text
m(t) = A0 + somme_k ak cos(2*pi*k*(t-t0)/P) + bk sin(2*pi*k*(t-t0)/P)
```

L'ajustement est fait par moindres carrés pondérés avec les incertitudes sur la
magnitude :

```text
sum(((m_obs - m_modele) / sigma_m)^2)
```

Le script affiche :

- la période utilisée ;
- l'époque `t0` ;
- l'ordre Fourier ;
- la magnitude moyenne `A0` ;
- le `chi2 réduit` ;
- l'amplitude pic-à-pic de la variation de magnitude estimée sur le modèle ;
- le temps de montée en pourcentage du cycle et en jours ;
- les coefficients `ak`, `bk` ;
- la forme amplitude-phase équivalente `Ak`, `phik` ;
- les paramètres Fourier utiles pour les RR Lyrae : `R21`, `R31`, `phi21`, `phi31`.

Il génère aussi un graphique de la courbe repliée en phase, avec les barres
d'erreur et le fit Fourier superposé. Les points sont distingués par nuit avec
des couleurs et des symboles différents, et une légende indique la nuit
correspondante.

Les paramètres RR Lyrae sont calculés ainsi :

```text
R21 = A2 / A1
R31 = A3 / A1
phi21 = phi2 - 2*phi1
phi31 = phi3 - 3*phi1
```

Les phases composées sont ramenées dans l'intervalle `[0, 2*pi[`.

Par convention pour les RR Lyrae, l'époque `T0` est prise au maximum de
lumière, c'est-à-dire au minimum de magnitude. Si `--epoch` est fourni, cette
valeur est toujours prioritaire. Sinon, `T0` est estimé automatiquement avec
`--epoch-method`.

Deux méthodes automatiques sont disponibles :

- `local-poly` : ajuste une parabole locale pondérée sur les mesures proches du
  maximum de lumière. C'est la méthode par défaut.
- `model` : prend le maximum de lumière du modèle Fourier ajusté.

L'amplitude est calculée sur le modèle ajusté :

```text
amplitude = magnitude_max - magnitude_min
```

Le temps de montée correspond ici à la montée en luminosité, donc au trajet
entre :

```text
magnitude maximale -> magnitude minimale
luminosité minimale -> luminosité maximale
```

Il est affiché comme fraction du cycle, en pourcentage, et comme durée en jours :

```text
duree_montee = fraction_montee * periode
```

### Exemple avec la période connue

```powershell
uv run python .\fourier_rrlyrae.py 0.293162 --order 4 --plot .\output\folded_light_curve.png
```

### Paramètres de ligne de commande

```text
period
```

Période connue, en jours. C'est le seul argument obligatoire.

```text
--order {3..10}
```

Ordre de la série de Fourier. Valeur par défaut : `4`.

```text
--epoch EPOCH
```

Époque `t0` en JD. Si fourni, cette valeur est prioritaire.

```text
--epoch-method {local-poly,model}
```

Méthode automatique pour estimer `t0` si `--epoch` est absent. Valeur par
défaut : `local-poly`.

```text
--data DATA
```

Dossier contenant les fichiers CSV. Valeur par défaut : `data`.

```text
--pattern PATTERN
```

Motif des fichiers à lire. Valeur par défaut : `*.csv`.

```text
--plot PLOT
```

Chemin du graphique PNG à produire. Valeur par défaut :
`output/folded_light_curve.png`.

```text
--show
```

Affiche aussi le graphique dans une fenêtre Matplotlib.

```text
--label-mode {date,file}
```

Choix des étiquettes de nuit dans la légende. `date` utilise le début du nom de
fichier avant le premier `_`, par exemple `2026-06-18`. `file` utilise le nom
complet du fichier sans extension. Valeur par défaut : `date`.

```text
--single-style
```

Trace tous les points avec la même couleur et le même symbole, au lieu de
distinguer les nuits.

## Script 2 : `find_period.py`

Ce script estime la période à partir des données.

Il teste une grille de périodes candidates. Pour chaque période, il ajuste une
série de Fourier pondérée, puis retient la période qui minimise le `chi2 réduit`.

Par défaut, la recherche est faite dans une plage typique RRc :

```text
0.2 j <= P <= 0.5 j
```

Le script génère deux graphiques :

- `output/period_search.png` : évolution du `chi2 réduit` selon la période candidate ;
- `output/best_period_folded_light_curve.png` : courbe repliée avec la meilleure période.

### Exemple simple

```powershell
uv run python .\find_period.py --order 4
```

### Exemple avec plage personnalisée

```powershell
uv run python .\find_period.py --min-period 0.2 --max-period 0.5 --order 4
```

### Exemple avec plus de résolution

```powershell
uv run python .\find_period.py --order 4 --samples 20000 --refine-samples 10000
```

### Paramètres de ligne de commande

```text
--min-period MIN_PERIOD
```

Période minimale testée, en jours. Valeur par défaut : `0.2`.

```text
--max-period MAX_PERIOD
```

Période maximale testée, en jours. Valeur par défaut : `0.5`.

```text
--samples SAMPLES
```

Nombre de périodes testées dans la grille grossière. Valeur par défaut :
`5000`.

```text
--refine-samples REFINE_SAMPLES
```

Nombre de périodes testées autour du meilleur résultat de la grille grossière.
Valeur par défaut : `3000`.

```text
--order {1..10}
```

Ordre Fourier utilisé pour la recherche de période. Valeur par défaut : `4`.

```text
--epoch EPOCH
```

Époque `t0` en JD. Si fourni, cette valeur est prioritaire.

```text
--epoch-method {local-poly,model}
```

Méthode automatique pour estimer `t0` du fit final si `--epoch` est absent.
Valeur par défaut : `local-poly`.

```text
--data DATA
```

Dossier contenant les fichiers CSV. Valeur par défaut : `data`.

```text
--pattern PATTERN
```

Motif des fichiers à lire. Valeur par défaut : `*.csv`.

```text
--top TOP
```

Nombre de meilleures périodes candidates affichées. Valeur par défaut : `10`.

```text
--periodogram-plot PERIODOGRAM_PLOT
```

Chemin du graphique de recherche de période. Valeur par défaut :
`output/period_search.png`.

```text
--folded-plot FOLDED_PLOT
```

Chemin du graphique de courbe repliée avec la meilleure période. Valeur par
défaut : `output/best_period_folded_light_curve.png`.

```text
--show
```

Affiche aussi les graphiques dans une fenêtre Matplotlib.

```text
--label-mode {date,file}
```

Choix des étiquettes de nuit dans la légende de la courbe repliée. Valeur par
défaut : `date`.

```text
--single-style
```

Trace tous les points de la courbe repliée avec la même couleur et le même
symbole.

## Script 3 : `find_period_gls.py`

Ce script estime la période avec la méthode GLS, Generalized Lomb-Scargle, en
utilisant `astropy.timeseries.LombScargle`.

La méthode GLS ajuste pour chaque fréquence candidate un modèle sinusoïdal avec
moyenne flottante, en tenant compte des incertitudes de magnitude. Le script
retient la période qui maximise la puissance GLS.

Comme pour les autres scripts, il lit les fichiers CSV du dossier `data/`.

Il génère deux graphiques :

- `output/gls_periodogram.png` : périodogramme GLS ;
- `output/gls_best_period_folded_light_curve.png` : courbe repliée avec la meilleure période GLS.

La courbe repliée finale est affichée avec un fit Fourier, dont l'ordre est
choisi avec `--fit-order`.

### Exemple simple

```powershell
uv run python .\find_period_gls.py --fit-order 4
```

### Exemple avec sorties explicites

```powershell
uv run python .\find_period_gls.py --fit-order 4 --periodogram-plot .\output\gls_periodogram.png --folded-plot .\output\gls_best_period_folded_light_curve.png
```

### Paramètres de ligne de commande

```text
--min-period MIN_PERIOD
```

Période minimale testée, en jours. Valeur par défaut : `0.2`.

```text
--max-period MAX_PERIOD
```

Période maximale testée, en jours. Valeur par défaut : `0.5`.

```text
--samples SAMPLES
```

Nombre de fréquences testées dans la grille grossière. Valeur par défaut :
`10000`.

```text
--refine-samples REFINE_SAMPLES
```

Nombre de fréquences testées autour du meilleur résultat de la grille grossière.
Valeur par défaut : `5000`.

```text
--fit-order {3..10}
```

Ordre Fourier utilisé pour tracer la courbe repliée finale. Valeur par défaut :
`4`.

```text
--epoch EPOCH
```

Époque `t0` du fit Fourier final en JD. Si fourni, cette valeur est prioritaire.

```text
--epoch-method {local-poly,model}
```

Méthode automatique pour estimer `t0` du fit final si `--epoch` est absent.
Valeur par défaut : `local-poly`.

```text
--data DATA
```

Dossier contenant les fichiers CSV. Valeur par défaut : `data`.

```text
--pattern PATTERN
```

Motif des fichiers à lire. Valeur par défaut : `*.csv`.

```text
--top TOP
```

Nombre de meilleures périodes candidates affichées. Valeur par défaut : `10`.

```text
--periodogram-plot PERIODOGRAM_PLOT
```

Chemin du graphique GLS. Valeur par défaut : `output/gls_periodogram.png`.

```text
--folded-plot FOLDED_PLOT
```

Chemin du graphique de courbe repliée avec la meilleure période GLS. Valeur par
défaut : `output/gls_best_period_folded_light_curve.png`.

```text
--show
```

Affiche aussi les graphiques dans une fenêtre Matplotlib.

```text
--label-mode {date,file}
```

Choix des étiquettes de nuit dans la légende de la courbe repliée. Valeur par
défaut : `date`.

```text
--single-style
```

Trace tous les points de la courbe repliée avec la même couleur et le même
symbole.

## Notes

Les durées sont toujours affichées en jours. Quand la valeur est inférieure à
`1` jour, une conversion en heures est ajoutée pour faciliter la lecture.

Pour les magnitudes, l'axe vertical des graphiques est inversé : les magnitudes
plus faibles apparaissent plus haut, comme c'est l'usage en photométrie.

La recherche de période par grille peut présenter des minima secondaires ou des
alias. Il est donc utile de comparer le résultat avec la période attendue, la
couverture temporelle des observations et l'aspect de la courbe repliée.

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fourier_rrlyrae import (
    fit_fourier_series,
    format_duration_days,
    plot_folded_light_curve,
    print_fit_summary,
    read_light_curve_with_labels,
)


def period_grid(min_period: float, max_period: float, samples: int) -> np.ndarray:
    """Return an evenly spaced frequency grid converted to periods."""
    if min_period <= 0 or max_period <= 0:
        raise ValueError("Les periodes doivent etre positives.")
    if min_period >= max_period:
        raise ValueError("min_period doit etre inferieure a max_period.")
    if samples < 2:
        raise ValueError("samples doit etre >= 2.")

    min_frequency = 1.0 / max_period
    max_frequency = 1.0 / min_period
    frequencies = np.linspace(min_frequency, max_frequency, samples)
    return 1.0 / frequencies


def evaluate_periods(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    periods: np.ndarray,
    order: int,
    epoch: float | None,
) -> np.ndarray:
    """Fit every trial period and return reduced chi2 values."""
    reduced_chi2 = np.empty_like(periods)
    fit_epoch = float(np.min(jd)) if epoch is None else epoch
    for index, period in enumerate(periods):
        fit = fit_fourier_series(jd, mag, err, float(period), order, fit_epoch)
        reduced_chi2[index] = float(fit["reduced_chi2"])
    return reduced_chi2


def refine_around_best(
    best_period: float,
    min_period: float,
    max_period: float,
    coarse_periods: np.ndarray,
    refine_samples: int,
) -> np.ndarray:
    """Build a narrower period grid around the best coarse value."""
    if refine_samples <= 0:
        return np.array([], dtype=float)

    sorted_periods = np.sort(coarse_periods)
    local_step = float(np.median(np.diff(sorted_periods)))
    half_width = max(10.0 * local_step, best_period * 0.002)
    refined_min = max(min_period, best_period - half_width)
    refined_max = min(max_period, best_period + half_width)
    return period_grid(refined_min, refined_max, refine_samples)


def plot_period_search(
    periods: np.ndarray,
    reduced_chi2: np.ndarray,
    best_period: float,
    output_path: Path,
    show: bool = False,
) -> None:
    order = np.argsort(periods)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(periods[order], reduced_chi2[order], color="tab:blue", linewidth=1.2)
    ax.axvline(best_period, color="tab:red", linewidth=1.8, label=f"P = {best_period:.8f} j")
    ax.set_xlabel("Periode candidate (jours)")
    ax.set_ylabel("Chi2 reduit")
    ax.set_title("Recherche de periode par ajustement Fourier")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def print_best_candidates(periods: np.ndarray, reduced_chi2: np.ndarray, count: int) -> None:
    print("Meilleures periodes candidates:")
    for rank, index in enumerate(np.argsort(reduced_chi2)[:count], start=1):
        print(
            f"  {rank:2d}. P = {format_duration_days(float(periods[index]))}"
            f"    frequence = {1.0 / periods[index]:.8f} cycles/j"
            f"    chi2 reduit = {reduced_chi2[index]:.4f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recherche la periode d'une RR Lyrae par grille de periodes et fit Fourier."
    )
    parser.add_argument(
        "--min-period",
        type=float,
        default=0.2,
        help="Periode minimale testee, en jours. Defaut: 0.2.",
    )
    parser.add_argument(
        "--max-period",
        type=float,
        default=0.5,
        help="Periode maximale testee, en jours. Defaut: 0.5.",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=5000,
        help="Nombre de periodes testees dans la grille grossiere. Defaut: 5000.",
    )
    parser.add_argument(
        "--refine-samples",
        type=int,
        default=3000,
        help="Nombre de periodes testees autour du meilleur resultat. Defaut: 3000.",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=4,
        choices=range(1, 11),
        metavar="{1..10}",
        help="Ordre Fourier utilise pour la recherche. Defaut: 4.",
    )
    parser.add_argument(
        "--epoch",
        type=float,
        default=None,
        help="Epoque t0 en JD. Si fourni, cette valeur est prioritaire.",
    )
    parser.add_argument(
        "--epoch-method",
        choices=("local-poly", "model"),
        default="local-poly",
        help="Methode automatique pour t0 du fit final si --epoch est absent. Defaut: local-poly.",
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("data"),
        help="Dossier contenant les CSV. Defaut: data.",
    )
    parser.add_argument(
        "--pattern",
        default="*.csv",
        help="Motif des fichiers a lire. Defaut: *.csv.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Nombre de meilleures periodes affichees. Defaut: 10.",
    )
    parser.add_argument(
        "--periodogram-plot",
        type=Path,
        default=Path("output/period_search.png"),
        help="Chemin du graphique de recherche de periode. Defaut: output/period_search.png.",
    )
    parser.add_argument(
        "--folded-plot",
        type=Path,
        default=Path("output/best_period_folded_light_curve.png"),
        help="Chemin du graphique replie avec la meilleure periode.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Affiche aussi les graphiques dans une fenetre matplotlib.",
    )
    parser.add_argument(
        "--label-mode",
        choices=("date", "file"),
        default="date",
        help="Etiquettes des nuits dans la legende: date ou nom de fichier. Defaut: date.",
    )
    parser.add_argument(
        "--single-style",
        action="store_true",
        help="Trace tous les points de la courbe repliee avec la meme couleur et le meme symbole.",
    )
    args = parser.parse_args()

    paths = sorted(args.data.glob(args.pattern))
    if not paths:
        raise SystemExit(f"Aucun fichier trouve: {args.data / args.pattern}")

    jd, mag, err, labels = read_light_curve_with_labels(paths, args.label_mode)

    coarse_periods = period_grid(args.min_period, args.max_period, args.samples)
    coarse_chi2 = evaluate_periods(jd, mag, err, coarse_periods, args.order, args.epoch)
    coarse_best = float(coarse_periods[np.argmin(coarse_chi2)])

    refined_periods = refine_around_best(
        coarse_best,
        args.min_period,
        args.max_period,
        coarse_periods,
        args.refine_samples,
    )
    if len(refined_periods) > 0:
        refined_chi2 = evaluate_periods(jd, mag, err, refined_periods, args.order, args.epoch)
        periods = np.concatenate([coarse_periods, refined_periods])
        reduced_chi2 = np.concatenate([coarse_chi2, refined_chi2])
    else:
        periods = coarse_periods
        reduced_chi2 = coarse_chi2

    best_index = int(np.argmin(reduced_chi2))
    best_period = float(periods[best_index])
    best_fit = fit_fourier_series(
        jd,
        mag,
        err,
        best_period,
        args.order,
        args.epoch,
        args.epoch_method,
    )

    print_best_candidates(periods, reduced_chi2, args.top)
    print()
    print_fit_summary(best_fit)

    plot_period_search(periods, reduced_chi2, best_period, args.periodogram_plot, args.show)
    plot_folded_light_curve(
        jd,
        mag,
        err,
        best_fit,
        args.folded_plot,
        args.show,
        labels,
        group_by_label=not args.single_style,
    )
    print()
    print(f"Graphique recherche periode: {args.periodogram_plot}")
    print(f"Graphique courbe repliée: {args.folded_plot}")


if __name__ == "__main__":
    main()

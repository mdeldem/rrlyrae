from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle

from fourier_rrlyrae import (
    fit_fourier_series,
    format_duration_days,
    plot_folded_light_curve,
    print_fit_summary,
    read_light_curve_with_labels,
)


def frequency_grid(min_period: float, max_period: float, samples: int) -> np.ndarray:
    """Return an evenly spaced frequency grid in cycles/day."""
    if min_period <= 0 or max_period <= 0:
        raise ValueError("Les periodes doivent etre positives.")
    if min_period >= max_period:
        raise ValueError("min_period doit etre inferieure a max_period.")
    if samples < 2:
        raise ValueError("samples doit etre >= 2.")

    min_frequency = 1.0 / max_period
    max_frequency = 1.0 / min_period
    return np.linspace(min_frequency, max_frequency, samples)


def refine_frequency_grid(
    best_period: float,
    min_period: float,
    max_period: float,
    coarse_frequencies: np.ndarray,
    refine_samples: int,
) -> np.ndarray:
    """Return a narrower frequency grid around the best coarse period."""
    if refine_samples <= 0:
        return np.array([], dtype=float)

    coarse_periods = 1.0 / coarse_frequencies
    sorted_periods = np.sort(coarse_periods)
    local_period_step = float(np.median(np.diff(sorted_periods)))
    half_width = max(10.0 * local_period_step, best_period * 0.002)
    refined_min_period = max(min_period, best_period - half_width)
    refined_max_period = min(max_period, best_period + half_width)
    return frequency_grid(refined_min_period, refined_max_period, refine_samples)


def compute_gls(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    frequencies: np.ndarray,
) -> np.ndarray:
    """Compute the Generalized Lomb-Scargle power with Astropy."""
    lomb_scargle = LombScargle(
        jd,
        mag,
        dy=err,
        fit_mean=True,
        center_data=True,
        nterms=1,
    )
    return lomb_scargle.power(
        frequencies,
        normalization="standard",
        method="auto",
    )


def plot_gls_periodogram(
    periods: np.ndarray,
    power: np.ndarray,
    best_period: float,
    output_path: Path,
    show: bool = False,
) -> None:
    order = np.argsort(periods)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(periods[order], power[order], color="tab:blue", linewidth=1.2)
    ax.axvline(best_period, color="tab:red", linewidth=1.8, label=f"P = {best_period:.8f} j")
    ax.set_xlabel("Periode candidate (jours)")
    ax.set_ylabel("Puissance GLS")
    ax.set_title("Recherche de periode par Generalized Lomb-Scargle")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def print_best_candidates(periods: np.ndarray, power: np.ndarray, count: int) -> None:
    print("Meilleures periodes candidates GLS:")
    for rank, index in enumerate(np.argsort(power)[::-1][:count], start=1):
        print(
            f"  {rank:2d}. P = {format_duration_days(float(periods[index]))}"
            f"    frequence = {1.0 / periods[index]:.8f} cycles/j"
            f"    puissance = {power[index]:.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recherche la periode d'une RR Lyrae avec Astropy LombScargle GLS."
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
        default=10000,
        help="Nombre de frequences testees dans la grille grossiere. Defaut: 10000.",
    )
    parser.add_argument(
        "--refine-samples",
        type=int,
        default=5000,
        help="Nombre de frequences testees autour du meilleur resultat. Defaut: 5000.",
    )
    parser.add_argument(
        "--fit-order",
        type=int,
        default=4,
        choices=range(3, 11),
        metavar="{3..10}",
        help="Ordre Fourier du graphique replie final. Defaut: 4.",
    )
    parser.add_argument(
        "--epoch",
        type=float,
        default=None,
        help="Epoque t0 du fit Fourier final en JD. Defaut: maximum de lumiere du modele.",
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
        help="Nombre de meilleures periodes candidates affichees. Defaut: 10.",
    )
    parser.add_argument(
        "--periodogram-plot",
        type=Path,
        default=Path("output/gls_periodogram.png"),
        help="Chemin du graphique GLS. Defaut: output/gls_periodogram.png.",
    )
    parser.add_argument(
        "--folded-plot",
        type=Path,
        default=Path("output/gls_best_period_folded_light_curve.png"),
        help="Chemin du graphique replie avec la meilleure periode GLS.",
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

    coarse_frequencies = frequency_grid(args.min_period, args.max_period, args.samples)
    coarse_power = compute_gls(jd, mag, err, coarse_frequencies)
    coarse_periods = 1.0 / coarse_frequencies
    coarse_best_period = float(coarse_periods[np.argmax(coarse_power)])

    refined_frequencies = refine_frequency_grid(
        coarse_best_period,
        args.min_period,
        args.max_period,
        coarse_frequencies,
        args.refine_samples,
    )
    if len(refined_frequencies) > 0:
        refined_power = compute_gls(jd, mag, err, refined_frequencies)
        frequencies = np.concatenate([coarse_frequencies, refined_frequencies])
        power = np.concatenate([coarse_power, refined_power])
    else:
        frequencies = coarse_frequencies
        power = coarse_power

    periods = 1.0 / frequencies
    best_index = int(np.argmax(power))
    best_period = float(periods[best_index])
    best_fit = fit_fourier_series(jd, mag, err, best_period, args.fit_order, args.epoch)

    print_best_candidates(periods, power, args.top)
    print()
    print_fit_summary(best_fit)

    plot_gls_periodogram(periods, power, best_period, args.periodogram_plot, args.show)
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
    print(f"Graphique GLS: {args.periodogram_plot}")
    print(f"Graphique courbe repliée: {args.folded_plot}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_light_curve_with_labels(
    paths: list[Path],
    label_mode: str = "date",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read one or more CSV files with columns JD;Magnitude;ErreurMagnitude."""
    jd_parts = []
    mag_parts = []
    err_parts = []
    label_parts = []

    for path in paths:
        data = np.genfromtxt(
            path,
            delimiter=";",
            names=True,
            dtype=float,
            encoding="utf-8",
        )
        jd_parts.append(np.asarray(data["JD"], dtype=float))
        mag_parts.append(np.asarray(data["Magnitude"], dtype=float))
        err_parts.append(np.asarray(data["ErreurMagnitude"], dtype=float))
        label = path.stem if label_mode == "file" else path.stem.split("_")[0]
        labels = np.full(np.asarray(data["JD"]).shape, label, dtype=object)
        label_parts.append(labels)

    jd = np.concatenate(jd_parts)
    mag = np.concatenate(mag_parts)
    err = np.concatenate(err_parts)
    labels = np.concatenate(label_parts)

    valid = np.isfinite(jd) & np.isfinite(mag) & np.isfinite(err) & (err > 0)
    order = np.argsort(jd[valid])
    return jd[valid][order], mag[valid][order], err[valid][order], labels[valid][order]


def read_light_curve(paths: list[Path]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    jd, mag, err, _labels = read_light_curve_with_labels(paths)
    return jd, mag, err


def phase_from_jd(jd: np.ndarray, period: float, epoch: float | None = None) -> np.ndarray:
    """Return phases in [0, 1). If epoch is omitted, the first JD is used."""
    if epoch is None:
        epoch = float(np.min(jd))
    return ((jd - epoch) / period) % 1.0


def fourier_design_matrix(
    jd: np.ndarray,
    period: float,
    order: int,
    epoch: float | None = None,
) -> np.ndarray:
    """Build [1, cos(wt), sin(wt), ..., cos(nwt), sin(nwt)]."""
    if not 1 <= order:
        raise ValueError("order must be >= 1")
    if period <= 0:
        raise ValueError("period must be > 0")

    if epoch is None:
        epoch = float(np.min(jd))

    x = 2.0 * np.pi * (jd - epoch) / period
    columns = [np.ones_like(jd)]
    for k in range(1, order + 1):
        columns.append(np.cos(k * x))
        columns.append(np.sin(k * x))
    return np.column_stack(columns)


def solve_weighted_fourier(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    period: float,
    order: int,
    epoch: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    design = fourier_design_matrix(jd, period, order, epoch)
    weighted_design = design / err[:, None]
    weighted_mag = mag / err
    coeffs, *_ = np.linalg.lstsq(weighted_design, weighted_mag, rcond=None)

    model = design @ coeffs
    residuals = mag - model
    chi2 = float(np.sum((residuals / err) ** 2))
    dof = int(max(len(mag) - len(coeffs), 0))
    reduced_chi2 = chi2 / dof if dof > 0 else np.nan
    return coeffs, model, residuals, chi2, reduced_chi2


def epoch_of_maximum_light(
    coeffs: np.ndarray,
    period: float,
    epoch: float,
    samples: int = 20000,
) -> float:
    """Estimate T0 as the model maximum-light time, i.e. minimum magnitude."""
    phase = np.linspace(0.0, 1.0, samples, endpoint=False)
    jd = epoch + phase * period
    model_mag = predict_magnitude(jd, coeffs, period, epoch)
    brightest_index = int(np.argmin(model_mag))
    return float(jd[brightest_index])


def epoch_from_local_polynomial(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    period: float,
    reference_epoch: float,
    window_phase: float = 0.08,
) -> float:
    """Estimate maximum-light T0 from measurements with a local quadratic fit."""
    dt = ((jd - reference_epoch + 0.5 * period) % period) - 0.5 * period
    half_window_days = window_phase * period
    in_window = np.abs(dt) <= half_window_days

    if np.count_nonzero(in_window) < 5:
        brightest_index = int(np.argmin(mag))
        observed_epoch = float(jd[brightest_index])
        cycles_to_reference = np.rint((reference_epoch - observed_epoch) / period)
        return float(observed_epoch + cycles_to_reference * period)

    x = dt[in_window]
    y = mag[in_window]
    sigma = err[in_window]

    x_scale = half_window_days if half_window_days > 0 else period
    x_scaled = x / x_scale
    design = np.column_stack([np.ones_like(x_scaled), x_scaled, x_scaled**2])
    weighted_design = design / sigma[:, None]
    weighted_mag = y / sigma
    coeffs, *_ = np.linalg.lstsq(weighted_design, weighted_mag, rcond=None)

    curvature = coeffs[2]
    if curvature <= 0:
        return reference_epoch

    vertex_scaled = -coeffs[1] / (2.0 * curvature)
    if abs(vertex_scaled) > 1.0:
        return reference_epoch

    return float(reference_epoch + vertex_scaled * x_scale)


def automatic_epoch(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    period: float,
    order: int,
    epoch_method: str,
) -> float:
    initial_epoch = float(np.min(jd))
    initial_coeffs, *_ = solve_weighted_fourier(jd, mag, err, period, order, initial_epoch)

    if epoch_method == "model":
        return epoch_of_maximum_light(initial_coeffs, period, initial_epoch)
    if epoch_method == "local-poly":
        model_epoch = epoch_of_maximum_light(initial_coeffs, period, initial_epoch)
        return epoch_from_local_polynomial(jd, mag, err, period, model_epoch)
    raise ValueError(f"epoch_method inconnu: {epoch_method}")


def fit_fourier_series(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    period: float,
    order: int,
    epoch: float | None = None,
    epoch_method: str = "local-poly",
) -> dict[str, object]:
    """Weighted least-squares Fourier fit for a known period."""
    epoch_was_manual = epoch is not None
    if epoch is None:
        epoch = automatic_epoch(jd, mag, err, period, order, epoch_method)

    coeffs, model, residuals, chi2, reduced_chi2 = solve_weighted_fourier(
        jd, mag, err, period, order, float(epoch)
    )

    harmonic_amplitudes = np.empty(order)
    harmonic_phases = np.empty(order)
    for k in range(order):
        ak = coeffs[1 + 2 * k]
        bk = coeffs[2 + 2 * k]
        harmonic_amplitudes[k] = np.hypot(ak, bk)
        harmonic_phases[k] = np.arctan2(-bk, ak) % (2.0 * np.pi)

    return {
        "period": period,
        "order": order,
        "epoch": float(epoch),
        "epoch_method": "manual" if epoch_was_manual else epoch_method,
        "coeffs": coeffs,
        "model": model,
        "residuals": residuals,
        "chi2": chi2,
        "reduced_chi2": reduced_chi2,
        "harmonic_amplitudes": harmonic_amplitudes,
        "harmonic_phases": harmonic_phases,
    }


def predict_magnitude(
    jd: np.ndarray,
    coeffs: np.ndarray,
    period: float,
    epoch: float,
) -> np.ndarray:
    order = (len(coeffs) - 1) // 2
    design = fourier_design_matrix(jd, period, order, epoch)
    return design @ coeffs


def format_duration_days(days: float, decimals: int = 10) -> str:
    text = f"{days:.{decimals}f} jours"
    if abs(days) < 1.0:
        text += f" ({days * 24.0:.4f} heures)"
    return text


def light_curve_shape_parameters(
    fit: dict[str, object],
    samples: int = 20000,
) -> dict[str, float]:
    """Estimate amplitude and rising-branch duration from the fitted curve."""
    period = float(fit["period"])
    epoch = float(fit["epoch"])
    coeffs = np.asarray(fit["coeffs"])

    phase = np.linspace(0.0, 1.0, samples, endpoint=False)
    jd = epoch + phase * period
    model_mag = predict_magnitude(jd, coeffs, period, epoch)

    faintest_index = int(np.argmax(model_mag))
    brightest_index = int(np.argmin(model_mag))
    faintest_phase = float(phase[faintest_index])
    brightest_phase = float(phase[brightest_index])
    rise_phase_fraction = (brightest_phase - faintest_phase) % 1.0

    amplitude_mag = float(np.max(model_mag) - np.min(model_mag))
    return {
        "amplitude_mag": amplitude_mag,
        "faintest_phase": faintest_phase,
        "brightest_phase": brightest_phase,
        "rise_phase_fraction": rise_phase_fraction,
        "rise_duration_days": rise_phase_fraction * period,
        "rise_percentage": 100.0 * rise_phase_fraction,
    }


def print_fit_summary(fit: dict[str, object]) -> None:
    coeffs = np.asarray(fit["coeffs"])
    amplitudes = np.asarray(fit["harmonic_amplitudes"])
    phases = np.asarray(fit["harmonic_phases"])
    order = int(fit["order"])

    print(f"Periode: {format_duration_days(float(fit['period']))}")
    print(f"Epoque t0: {fit['epoch']:.8f} JD")
    print(f"Methode t0: {fit['epoch_method']}")
    print(f"Ordre Fourier: {order}")
    print(f"Magnitude moyenne A0: {coeffs[0]:.6f}")
    print(f"Chi2 reduit: {fit['reduced_chi2']:.4f}")
    shape = light_curve_shape_parameters(fit)
    print(f"Amplitude modele pic-a-pic: {shape['amplitude_mag']:.6f} mag")
    print(
        "Temps de montée: "
        f"{shape['rise_percentage']:.2f}% du cycle, "
        f"soit {format_duration_days(shape['rise_duration_days'], decimals=6)}"
    )
    print(
        "  de la phase "
        f"{shape['faintest_phase']:.4f} (magnitude max, luminosite min) "
        f"a {shape['brightest_phase']:.4f} (magnitude min, luminosite max)"
    )
    print()
    print("Coefficients du modele:")
    print("  m(t) = A0 + somme_k ak cos(2*pi*k*(t-t0)/P) + bk sin(...)")
    print(f"  A0 = {coeffs[0]: .8f}")
    for k in range(1, order + 1):
        ak = coeffs[1 + 2 * (k - 1)]
        bk = coeffs[2 + 2 * (k - 1)]
        print(f"  a{k} = {ak: .8f}    b{k} = {bk: .8f}")

    print()
    print("Forme amplitude-phase equivalente:")
    print("  m(t) = A0 + somme_k Ak cos(2*pi*k*(t-t0)/P + phik)")
    for k in range(1, order + 1):
        print(f"  A{k} = {amplitudes[k - 1]: .8f}    phi{k} = {phases[k - 1]: .8f} rad")

    if order >= 2 and amplitudes[0] > 0:
        print()
        print("Parametres Fourier RR Lyrae:")
        print(f"  R21 = {amplitudes[1] / amplitudes[0]: .8f}")
        phi21 = (phases[1] - 2.0 * phases[0]) % (2.0 * np.pi)
        print(f"  phi21 = {phi21: .8f} rad")
        if order >= 3:
            print(f"  R31 = {amplitudes[2] / amplitudes[0]: .8f}")
            phi31 = (phases[2] - 3.0 * phases[0]) % (2.0 * np.pi)
            print(f"  phi31 = {phi31: .8f} rad")


def plot_folded_light_curve(
    jd: np.ndarray,
    mag: np.ndarray,
    err: np.ndarray,
    fit: dict[str, object],
    output_path: Path,
    show: bool = False,
    labels: np.ndarray | None = None,
    group_by_label: bool = True,
) -> None:
    """Plot the folded light curve with error bars and the Fourier fit."""
    period = float(fit["period"])
    epoch = float(fit["epoch"])
    coeffs = np.asarray(fit["coeffs"])
    order = int(fit["order"])

    phase = phase_from_jd(jd, period, epoch)
    if labels is None or not group_by_label:
        labels = np.full(jd.shape, "Mesures", dtype=object)
    labels = np.asarray(labels)
    unique_labels = list(dict.fromkeys(labels.tolist()))
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fit_phase = np.linspace(0.0, 1.0, 800, endpoint=True)
    fit_jd = epoch + fit_phase * period
    fit_mag = predict_magnitude(fit_jd, coeffs, period, epoch)

    fig, ax = plt.subplots(figsize=(10, 6))
    for label_index, label in enumerate(unique_labels):
        mask = labels == label
        order_by_phase = np.argsort(phase[mask])
        label_phase = phase[mask][order_by_phase]
        label_mag = mag[mask][order_by_phase]
        label_err = err[mask][order_by_phase]
        color = colors[label_index % len(colors)]
        marker = markers[label_index % len(markers)]

        for offset in (0.0, 1.0):
            ax.errorbar(
                label_phase + offset,
                label_mag,
                yerr=label_err,
                fmt=marker,
                linestyle="none",
                markersize=4,
                elinewidth=0.6,
                capsize=1.5,
                alpha=0.42,
                color=color,
                label=str(label) if offset == 0.0 else None,
                zorder=2,
            )

    for offset in (0.0, 1.0):
        ax.plot(
            fit_phase + offset,
            fit_mag,
            color="tab:red",
            linewidth=1.8,
            label=f"Fit Fourier ordre {order}" if offset == 0.0 else None,
            zorder=5,
        )

    ax.set_xlabel("Phase")
    ax.set_ylabel("Magnitude")
    ax.set_title(f"Courbe de lumière repliée - P = {period:.6f} j")
    ax.set_xlim(0.0, 2.0)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.25)
    ax.legend(
        fontsize=8,
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
    )
    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ajuste une serie de Fourier ponderee sur une courbe de lumiere RR Lyrae."
    )
    parser.add_argument("period", type=float, help="Periode connue, en jours.")
    parser.add_argument(
        "--order",
        type=int,
        default=4,
        choices=range(3, 11),
        metavar="{3..10}",
        help="Ordre Fourier, entre 3 et 10. Defaut: 4.",
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
        help="Methode automatique pour t0 si --epoch est absent. Defaut: local-poly.",
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
        "--plot",
        type=Path,
        default=Path("output/folded_light_curve.png"),
        help="Chemin du graphique PNG a produire. Defaut: output/folded_light_curve.png.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Affiche aussi le graphique dans une fenetre matplotlib.",
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
        help="Trace tous les points avec la meme couleur et le meme symbole.",
    )
    args = parser.parse_args()

    paths = sorted(args.data.glob(args.pattern))
    if not paths:
        raise SystemExit(f"Aucun fichier trouve: {args.data / args.pattern}")

    jd, mag, err, labels = read_light_curve_with_labels(paths, args.label_mode)
    fit = fit_fourier_series(
        jd,
        mag,
        err,
        args.period,
        args.order,
        args.epoch,
        args.epoch_method,
    )
    print_fit_summary(fit)
    plot_folded_light_curve(
        jd,
        mag,
        err,
        fit,
        args.plot,
        args.show,
        labels,
        group_by_label=not args.single_style,
    )
    print()
    print(f"Graphique enregistre: {args.plot}")


if __name__ == "__main__":
    main()

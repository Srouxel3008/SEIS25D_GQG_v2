"""
Simple Seis2D plotting script.

Edit the USER SETTINGS section first, then run:

    python seis2d_plot.py


"""

from __future__ import annotations

import math  # Used only to choose a neat row/column layout for subplots.
import re  # Used to read parameter/frequency/iteration from file names.
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ModuleNotFoundError as exc:
    missing_name = exc.name
    raise SystemExit(
        f"Missing Python package: {missing_name}\n"
        "Install the plotting requirements with:\n"
        "    python -m pip install -r requirements.txt"
    ) from exc


# =============================================================================
# USER SETTINGS
# =============================================================================

INPUT_FOLDER = Path(
    r"C:\Users\sedar\Desktop\SEIS25D_GQG_github\Github_version\runs\Small_TTI\output_20260807_124004"
)

# Figures are saved here, inside INPUT_FOLDER when this is a relative path.
OUTPUT_FOLDER = Path("python_plots")

# Case-specific plotting controls.
WELL_X = 100
DEPTH_LIMITS = (-210, 10.0)  # (z_min, z_max), in metres
X_LIMITS = None  # example: (0.0, 4500.0); use None for automatic limits
TOP_BOUNDARY_ZONE_THICKNESS = 100.0  # to adjust receiver depth

# Export quality.
DPI = 150

# Spyder note:
# False is best for making many figures because each figure is saved then closed.
# Set True when you want figures to remain visible in Spyder's Plots pane.
SHOW_FIGURES = True

# Select which plot families to make.
PLOT_MODEL_MAPS = True
PLOT_VERTICAL_PROFILES = True
PLOT_GRADIENTS = False
PLOT_DIAGNOSTICS = False
PLOT_SPECTRAL_DATA = True
PLOT_RESIDUAL_HEATMAPS = False

# Spectral data plot controls. These read files such as GT0_1_0.txt,
# GTO_1_0.txt, and GO_1_1.txt, then plot real/imag versus receiver ID.
SPECTRAL_SOURCE_ID = 1
SPECTRAL_FREQUENCY = 1

# Limit plotting to a few parameters if wanted. Use [] for all discovered.
PARAMETERS_TO_PLOT = []  # examples: ["C33", "Q33"]

# Optional measured/true well profiles when plotting real data inversion results. Each file should have at least two
# columns: depth, value. Leave the dictionary empty if not needed.
TRUE_PROFILE_FILES = {
    # "C33": r"C:\path\to\S13_C33_QC_well.csv",
    # "Q33": r"C:\path\to\Qp_QC.txt",
}


# =============================================================================
# PARAMETER NAMES AND UNITS
# =============================================================================

PARAMETER_LABELS = {
    "rho": "rho [kg/m3]",
    "C11": "C11 [Pa]",
    "C13": "C13 [Pa]",
    "C33": "C33 [Pa]",
    "C44": "C44 [Pa]",
    "C66": "C66 [Pa]",
    "Q11": "Q11",
    "Q13": "Q13",
    "Q33": "Q33",
    "Q44": "Q44",
    "Q66": "Q66",
    "Vp": "Vp [m/s]",
    "Vs": "Vs [m/s]",
    "Qp": "Qp",
    "Qs": "Qs",
}


# =============================================================================
# FILE READING FUNCTIONS
# =============================================================================

def read_numeric_table(path: Path) -> np.ndarray:
    """Read a whitespace or comma separated numeric file."""
    path = Path(path)
    try:
        data = np.genfromtxt(path, comments="#")
    except ValueError:
        data = np.genfromtxt(path, delimiter=",", comments="#")
    if data.ndim == 1 and data.size <= 1:
        data = np.genfromtxt(path, delimiter=",", comments="#")
    if np.size(data) > 0 and np.all(~np.isfinite(data)):
        data = np.genfromtxt(path, delimiter=",", comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def read_model_grid(path: Path, missing_value: float = 1.0e20) -> dict:
    """Read mI/mT/PCBB/GRAD style files."""
    data = read_numeric_table(path)
    x = data[1:, 0]
    z = data[0, 1:]
    values = data[1:, 1:].astype(float)
    values[values >= missing_value] = np.nan
    return {"x": x, "z": z, "values": values}


def make_grid(x: np.ndarray, z: np.ndarray, values: np.ndarray) -> dict:
    """Create the same simple grid dictionary used by read_model_grid()."""
    return {"x": x, "z": z, "values": values}


def read_topography(folder: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """Read m_Top.dat when available."""
    path = folder / "m_Top.dat"
    if not path.exists():
        return None
    data = read_numeric_table(path)
    return data[:, -2], data[:, -1]


def read_sr_geometry(folder: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Read source and receiver coordinates from m_SR.dat."""
    path = folder / "m_SR.dat"
    if not path.exists():
        path = folder / "m_S-R.dat"
    if not path.exists():
        return None
    data = read_numeric_table(path)
    is_source = data[:, 0] == 1
    x = data[:, -2]
    z = data[:, -1]
    return x[is_source], z[is_source], x[~is_source], z[~is_source]


def read_profile_at_x(grid_path: Path, x_target: float) -> tuple[np.ndarray, np.ndarray]:
    """Extract the nearest vertical profile from a model grid."""
    grid = read_model_grid(grid_path)
    ix = int(np.nanargmin(np.abs(grid["x"] - x_target)))
    return grid["z"], grid["values"][ix, :]


def read_true_profile(path: Path, depth_limits: tuple[float, float] | None) -> tuple[np.ndarray, np.ndarray]:
    """Read optional well/QC profile files with columns depth, value."""
    data = read_numeric_table(path)
    z = data[:, 0]
    values = data[:, 1]
    if depth_limits is not None:
        zmin, zmax = depth_limits
        keep = (z >= zmin) & (z <= zmax)
        z = z[keep]
        values = values[keep]
    return z, values


def read_go_or_residual(path: Path, max_rows: int | None = None) -> dict[str, np.ndarray]:
    """Read GO/GT0/GTO or residual file columns into named arrays."""
    data = read_numeric_table(path)
    if max_rows is not None:
        data = data[:max_rows, :]

    if data.shape[1] >= 14:
        # GO/GT format: 4=sid, 5=scomp, 6=rid, 7=rcomp, 12=real, 13=imag, 14=amp.
        sid, scomp, rid, rcomp = data[:, 3], data[:, 4], data[:, 5], data[:, 6]
        real, imag, amp = data[:, 11], data[:, 12], data[:, 13]
    else:
        # Residual format: 3=sid, 4=scomp, 5=rid, 6=rcomp, 11=real, 12=imag.
        sid, scomp, rid, rcomp = data[:, 2], data[:, 3], data[:, 4], data[:, 5]
        real, imag = data[:, 10], data[:, 11]
        amp = np.hypot(real, imag)

    return {
        "sid": sid,
        "scomp": scomp,
        "rid": rid,
        "rcomp": rcomp,
        "real": real,
        "imag": imag,
        "amp": amp,
    }


# =============================================================================
# FILE DISCOVERY FUNCTIONS
# =============================================================================

def clean_parameter_name(text: str) -> str:
    """Turn '{C33}', '{{C33}', or 'C33}' into 'C33'."""
    return text.replace("{", "").replace("}", "")


def parse_iteration_number(label: str) -> int:
    """Turn 'IT03' into 3; return 0 when no number is found."""
    match = re.search(r"(\d+)$", label)
    return int(match.group(1)) if match else 0


def parse_pcbb_file(path: Path) -> dict | None:
    """Parse names like PCBB_{{C33}_05.65_IT03}.dat."""
    match = re.match(
        r"^(PCBB|PCBA)_+\{*(?P<param>[^}_]+)\}*_(?P<freq>[^_]+)_(?P<it>IT\d+)\}*\.dat$", path.name)
    if not match:
        return None
    return {
        "path": path,
        "family": match.group(1),
        "parameter": clean_parameter_name(match.group("param")),
        "frequency": match.group("freq"),
        "iteration_label": match.group("it"),
        "iteration_number": parse_iteration_number(match.group("it")),
    }


def parse_gradient_file(path: Path) -> dict | None:
    """Parse names like GRADr_{C33}_05.65_IT02.dat."""
    match = re.match(
        r"^(?P<family>GRAD[^_]*?)_+\{*(?P<param>[^}_]+)\}*_(?P<freq>[^_]+)_(?P<it>IT\d+)\.dat$", path.name)
    if not match:
        return None
    return {
        "path": path,
        "family": match.group("family"),
        "parameter": clean_parameter_name(match.group("param")),
        "frequency": match.group("freq"),
        "iteration_label": match.group("it"),
        "iteration_number": parse_iteration_number(match.group("it")),
    }


def find_parameter_files(folder: Path, prefix: str) -> dict[str, Path]:
    """Find starting or true model files, for example prefix='mI' or 'mT'."""
    found: dict[str, Path] = {}
    for path in folder.glob(f"{prefix}_*.dat"):
        match = re.search(r"\{([^}]+)\}", path.name)
        if match:
            found[clean_parameter_name(match.group(1))] = path
    return found


def find_final_models(folder: Path) -> list[dict]:
    """Find all PCBB files and parse parameter/frequency/iteration from names."""
    out = []
    for path in folder.glob("PCBB_*.dat"):
        parsed = parse_pcbb_file(path)
        if parsed is not None:
            out.append(parsed)
    return sorted(out, key=lambda f: (f["parameter"], numeric_token(f["frequency"]), f["iteration_number"]))


def find_gradient_files(folder: Path) -> list[dict]:
    """Find all GRAD*.dat files."""
    out = []
    for path in folder.glob("GRAD*.dat"):
        parsed = parse_gradient_file(path)
        if parsed is not None:
            out.append(parsed)
    return sorted(out, key=lambda f: (f["family"], f["parameter"], numeric_token(f["frequency"]), f["iteration_number"]))


def find_data_files(folder: Path) -> dict[str, list[Path]]:
    """Find common Seis2D data/residual files."""
    return {
        "go": sorted(folder.glob("GO_*.txt")),
        "gt0": sorted(folder.glob("GT0_*.txt")),
        "gto": sorted(folder.glob("GTO_*.txt")),
        "residual": sorted(folder.glob("out_resid*.txt")),
    }


def find_spectral_files(folder: Path, frequency: float | int | str) -> list[Path]:
    """Find GT0/GTO/GO files for one frequency number."""
    files = find_data_files(folder)
    candidates = files["gt0"] + files["gto"] + files["go"]
    return [path for path in candidates if data_file_matches_frequency(path, frequency)]


def find_diagnostic_files(folder: Path) -> dict[str, list[Path]]:
    """Find inversion diagnostic text files."""
    return {
        "diag": sorted(folder.glob("out_diag*.txt")) + sorted(folder.glob("out_diagt*.txt")),
        "linesearch": sorted(folder.glob("out_linesearch*.txt")),
        "lbfgs": sorted(folder.glob("out_lbfgs_hist*.txt")),
    }


def numeric_token(text: str) -> float:
    """Sort strings by their first number when possible."""
    match = re.search(r"[-+]?\d*\.?\d+", str(text))
    return float(match.group(0)) if match else math.inf


def data_file_matches_frequency(path: Path, frequency: float | int | str) -> bool:
    """Check whether a GO/GT0/GTO file name starts with the selected frequency."""
    wanted = float(frequency)
    match = re.match(r"^(?:GO|GT0|GTO)_(?P<freq>[-+]?\d*\.?\d+)_", path.name)
    if not match:
        return False
    return abs(float(match.group("freq")) - wanted) <= max(1e-8, abs(wanted) * 1e-8)


def group_iteration_files(files: list[dict]) -> dict[tuple[str, str, str], list[dict]]:
    """Group files by family, parameter, and frequency."""
    groups: dict[tuple[str, str, str], list[dict]] = {}
    for item in files:
        key = (item["family"], item["parameter"], item["frequency"])
        groups.setdefault(key, []).append(item)
    for key in groups:
        groups[key].sort(key=lambda f: f["iteration_number"])
    return groups


def first_existing(candidates: list[Path]) -> Path | None:
    """Return the first path in candidates that exists."""
    for path in candidates:
        if path.exists():
            return path
    return None


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def setup_matplotlib() -> None:
    """Use publication-style defaults without depending on MATLAB."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.grid": False,
            "figure.dpi": 120,
            "savefig.dpi": DPI,
        }
    )


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    """Save and close a figure."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    if SHOW_FIGURES:
        fig.show()
    else:
        plt.close(fig)
    print(f"Saved: {output_path}")


def parameter_label(parameter: str) -> str:
    """Display label for a parameter."""
    return PARAMETER_LABELS.get(parameter, parameter)


def apply_common_axes(ax: plt.Axes) -> None:
    """Apply case-specific axes limits and depth orientation."""
    if X_LIMITS is not None:
        ax.set_xlim(*X_LIMITS)
    if DEPTH_LIMITS is not None:
        ax.set_ylim(*DEPTH_LIMITS)
    ax.set_xlabel("X-distance [m]")
    ax.set_ylabel("Depth [m]")


def plot_grid_image(
    ax: plt.Axes,
    grid: dict,
    label: str,
    cmap: str = "viridis",
    clim: tuple[float, float] | None = None,
) -> None:
    """Draw a model grid with a colorbar."""
    mesh = ax.pcolormesh(grid["x"], grid["z"],
                         grid["values"].T, shading="auto", cmap=cmap)
    if clim is not None:
        mesh.set_clim(*clim)
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label(label)
    apply_common_axes(ax)


def overlay_topography_and_acquisition(ax: plt.Axes, folder: Path) -> None:
    """Overlay m_Top and m_SR/m_S-R when available."""
    topo = read_topography(folder)
    if topo is not None:
        ax.plot(topo[0], topo[1], color="cyan", linewidth=1.4)
        # if TOP_BOUNDARY_ZONE_THICKNESS > 0:
        #     ax.fill_between(
        #         topo[0],
        #         topo[1],
        #         topo[1] - TOP_BOUNDARY_ZONE_THICKNESS,
        #         color="cyan",
        #         alpha=0.12,
        #         linewidth=0,
        #     )

    geom = read_sr_geometry(folder)
    if geom is not None:
        x_src, z_src, x_rec, z_rec = geom
        z_src= z_src + TOP_BOUNDARY_ZONE_THICKNESS
        z_rec= z_rec + TOP_BOUNDARY_ZONE_THICKNESS


        ax.scatter(x_rec, z_rec, s=12, c="limegreen",
                   edgecolors="black", linewidths=0.2, label="Receivers")
        ax.scatter(x_src, z_src, s=45, c="red", marker="*",
                   edgecolors="black", linewidths=0.3, label="Sources")
        ax.legend(loc="best", fontsize=8)


def plot_model_maps(folder: Path, output_folder: Path, final_models: list[dict]) -> None:
    """Plot true/start/final model maps and simple differences."""
    true_models = find_parameter_files(folder, "mT")
    start_models = find_parameter_files(folder, "mI")
    groups = group_iteration_files(final_models)

    for (_family, parameter, frequency), files in groups.items():
        if PARAMETERS_TO_PLOT and parameter not in PARAMETERS_TO_PLOT:
            continue

        final_file = files[-1]["path"]
        true_file = true_models.get(parameter)
        start_file = start_models.get(parameter)
        reference_file = true_file or start_file
        if reference_file is None:
            print(f"Skipping model comparison for {
                  parameter}: no mT/mI reference found.")
            continue

        ref_grid = read_model_grid(reference_file)
        final_grid = read_model_grid(final_file)
        diff = final_grid["values"] - ref_grid["values"]
        percent = 100.0 * diff / ref_grid["values"]

        fig, axes = plt.subplots(2, 2, figsize=(
            10, 7), constrained_layout=True)
        label = parameter_label(parameter)
        finite_ref = ref_grid["values"][np.isfinite(ref_grid["values"])]
        clim = (float(np.nanmin(finite_ref)), float(
            np.nanmax(finite_ref))) if finite_ref.size else None

        plot_grid_image(axes[0, 0], ref_grid, label, cmap="jet", clim=clim)
        axes[0, 0].set_title("Reference")
        plot_grid_image(axes[0, 1], final_grid, label, cmap="jet", clim=clim)
        axes[0, 1].set_title(f"Final {frequency} Hz, {
                             files[-1]['iteration_label']}")

        plot_grid_image(axes[1, 0], make_grid(final_grid["x"], final_grid["z"],
                        percent), "Difference [%]", cmap="seismic", clim=(-10, 10))
        axes[1, 0].set_title("Percent difference")
        vmax = robust_abs_limit(diff)
        plot_grid_image(axes[1, 1], make_grid(
            final_grid["x"], final_grid["z"], diff), label, cmap="seismic", clim=(-vmax, vmax))
        axes[1, 1].set_title("Absolute difference")

        for ax in axes.ravel():
            overlay_topography_and_acquisition(ax, folder)

        out = output_folder / f"model_compare_{parameter}_{frequency}.png"
        save_figure(fig, out)


def plot_vertical_profiles(folder: Path, output_folder: Path, final_models: list[dict]) -> None:
    """Plot well-location vertical profiles for each parameter/frequency."""
    true_models = find_parameter_files(folder, "mT")
    start_models = find_parameter_files(folder, "mI")
    groups = group_iteration_files(final_models)

    for (_family, parameter, frequency), files in groups.items():
        if PARAMETERS_TO_PLOT and parameter not in PARAMETERS_TO_PLOT:
            continue

        fig, ax = plt.subplots(figsize=(6, 8), constrained_layout=True)
        profiles_for_xlim = []

        true_file = true_models.get(parameter)
        if true_file is not None:
            z, values = read_profile_at_x(true_file, WELL_X)
            ax.plot(values, z, "k-", linewidth=1.6, label="True")
            profiles_for_xlim.append(values)

        true_profile_path = TRUE_PROFILE_FILES.get(parameter)
        if true_profile_path:
            z, values = read_true_profile(
                Path(true_profile_path), DEPTH_LIMITS)
            ax.plot(values, z, "k:", linewidth=1.6, label="Well/QC")
            profiles_for_xlim.append(values)

        start_file = start_models.get(parameter)
        if start_file is not None:
            z, values = read_profile_at_x(start_file, WELL_X)
            ax.plot(values, z, "k--", linewidth=1.4, label="Start")
            profiles_for_xlim.append(values)

        colors = plt.cm.hsv(np.linspace(
            0, 1, max(len(files), 1), endpoint=False))
        for color, item in zip(colors, files):
            z, values = read_profile_at_x(item["path"], WELL_X)
            ax.plot(values, z, color=color, linewidth=1.1,
                    label=f"Iter {item['iteration_number']}")
            profiles_for_xlim.append(values)

        ax.set_xlabel(parameter_label(parameter))
        ax.set_ylabel("Depth [m]")
        if DEPTH_LIMITS is not None:
            ax.set_ylim(*DEPTH_LIMITS)
        set_sensible_profile_xlim(ax, profiles_for_xlim)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=8)
        ax.set_title(f"{parameter} at x={WELL_X:g} m, {frequency} Hz")

        out = output_folder / \
            f"profile_{parameter}_{frequency}_x{WELL_X:g}.png"
        save_figure(fig, out)


def plot_gradient_series(folder: Path, output_folder: Path, gradient_files: list[dict]) -> None:
    """Plot all iterations for each gradient family/parameter/frequency."""
    groups = group_iteration_files(gradient_files)
    for (family, parameter, frequency), files in groups.items():
        if PARAMETERS_TO_PLOT and parameter not in PARAMETERS_TO_PLOT:
            continue

        n = len(files)
        ncols = int(math.ceil(math.sqrt(n)))
        nrows = int(math.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(
            4.0 * ncols, 3.2 * nrows), constrained_layout=True)
        axes = np.atleast_1d(axes).ravel()

        all_values = []
        for item in files:
            vals = read_model_grid(item["path"])["values"]
            all_values.append(vals[np.isfinite(vals)])
        all_values = np.concatenate([v for v in all_values if v.size]) if any(
            v.size for v in all_values) else np.array([1.0])
        vmax = float(np.nanpercentile(np.abs(all_values), 99))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0

        for ax, item in zip(axes, files):
            grid = read_model_grid(item["path"])
            plot_grid_image(ax, grid, parameter_label(
                parameter), cmap="seismic", clim=(-vmax, vmax))
            overlay_topography_and_acquisition(ax, folder)
            ax.set_title(item["iteration_label"])

        for ax in axes[n:]:
            ax.axis("off")

        fig.suptitle(f"{family} gradient for {parameter}, {frequency} Hz")
        out = output_folder / "Grad" / f"{family}_{parameter}_{frequency}.png"
        save_figure(fig, out)


def plot_diagnostics(folder: Path, output_folder: Path) -> None:
    """Plot simple diagnostics from out_diag.txt."""
    files = find_diagnostic_files(folder)
    if not files["diag"]:
        print("No out_diag*.txt files found.")
        return

    for diag_file in files["diag"]:
        data = read_numeric_table(diag_file)
        if data.shape[1] < 3:
            continue
        freq = data[:, 0]
        iteration = data[:, 1]
        cost = data[:, 2]
        grad = data[:, 3] if data.shape[1] >= 4 else None
        residual = data[:, 4] if data.shape[1] >= 5 else None

        rows = 3 if residual is not None else 2
        fig, axes = plt.subplots(rows, 1, figsize=(
            7, 3.0 * rows), constrained_layout=True, sharex=True)
        axes = np.atleast_1d(axes)

        for f in unique_stable(freq):
            mask = freq == f
            order = np.argsort(iteration[mask])
            it = iteration[mask][order]
            c = cost[mask][order]
            base = c[0] if c[0] != 0 else 1.0
            axes[0].semilogy(it, c / base, "-o",
                             markersize=3, label=f"{f:g} Hz")
            if grad is not None:
                axes[1].semilogy(it, grad[mask][order], "-o", markersize=3)
            if residual is not None:
                axes[2].plot(it, residual[mask][order], "-o", markersize=3)

        axes[0].set_ylabel("Normalized misfit")
        axes[0].legend(loc="best", fontsize=8)
        if grad is not None:
            axes[1].set_ylabel("Gradient norm")
        if residual is not None:
            axes[2].set_ylabel("Residual RMS")
        axes[-1].set_xlabel("Iteration")
        for ax in axes:
            ax.grid(True, alpha=0.25)

        out = output_folder / f"diagnostic_{diag_file.stem}.png"
        save_figure(fig, out)


def plot_residual_heatmaps(folder: Path, output_folder: Path) -> None:
    """Plot real/imag/amplitude heatmaps from the first residual or GT file."""
    files = find_data_files(folder)
    candidates = files["residual"] or files["gt0"] or files["gto"] or files["go"]
    if not candidates:
        print("No residual, GO, GT0, or GTO files found for heatmaps.")
        return

    data_file = candidates[0]
    data = read_go_or_residual(data_file)
    sids = unique_stable(data["sid"])
    rids = unique_stable(data["rid"])
    scomps = unique_stable(data["scomp"])
    rcomps = unique_stable(data["rcomp"])

    ncols = max(1, len(scomps) * len(rcomps))
    fig, axes = plt.subplots(3, ncols, figsize=(
        3.4 * ncols, 8.0), constrained_layout=True)
    axes = np.asarray(axes).reshape(3, ncols)

    col = 0
    for scomp in scomps:
        for rcomp in rcomps:
            mask = (data["scomp"] == scomp) & (data["rcomp"] == rcomp)
            matrices = [
                make_sid_rid_matrix(
                    data["sid"][mask], data["rid"][mask], data["real"][mask], sids, rids),
                make_sid_rid_matrix(
                    data["sid"][mask], data["rid"][mask], data["imag"][mask], sids, rids),
                make_sid_rid_matrix(
                    data["sid"][mask], data["rid"][mask], data["amp"][mask], sids, rids),
            ]
            labels = ["real", "imag", "amplitude"]
            cmaps = ["seismic", "seismic", "viridis"]
            for row in range(3):
                ax = axes[row, col]
                image = ax.imshow(
                    matrices[row].T, aspect="auto", origin="upper", cmap=cmaps[row])
                fig.colorbar(image, ax=ax, label=labels[row])
                ax.set_xlabel("Source ID")
                ax.set_ylabel("Receiver ID")
                ax.set_title(f"s{scomp:g} r{rcomp:g}" if row == 0 else "")
            col += 1

    out = output_folder / f"heatmap_{data_file.stem}.png"
    save_figure(fig, out)


def plot_spectral_data(folder: Path, output_folder: Path, source_id: int, frequency: float | int | str) -> None:
    """Plot real and imaginary traces for one source and one frequency."""
    spectral_files = find_spectral_files(folder, frequency)
    if not spectral_files:
        print(f"No GT0/GTO/GO files found for frequency {frequency}.")
        return

    loaded = []
    receiver_components = []
    for path in spectral_files:
        data = read_go_or_residual(path)
        source_mask = data["sid"] == source_id
        if not np.any(source_mask):
            continue
        loaded.append((path, data))
        receiver_components.extend(list(unique_stable(data["rcomp"][source_mask])))

    if not loaded:
        print(f"No spectral rows found for source ID {source_id} at frequency {frequency}.")
        return

    receiver_components = unique_stable(np.asarray(receiver_components, dtype=float))
    ncols = len(receiver_components)
    fig, axes = plt.subplots(2, ncols, figsize=(4.2 * ncols, 7.0), constrained_layout=True, squeeze=False)

    for col, receiver_comp in enumerate(receiver_components):
        ax_real = axes[0, col]
        ax_imag = axes[1, col]

        for path, data in loaded:
            mask = (data["sid"] == source_id) & (data["rcomp"] == receiver_comp)
            if not np.any(mask):
                continue
            receiver_id = data["rid"][mask]
            order = np.argsort(receiver_id)
            label = path.stem

            ax_real.plot(receiver_id[order], data["real"][mask][order], "-", linewidth=1.2, label=label)
            ax_imag.plot(receiver_id[order], data["imag"][mask][order], "-", linewidth=1.2, label=label)

        ax_real.set_title(f"Receiver component {receiver_comp:g}")
        ax_real.set_ylabel("Real")
        ax_imag.set_ylabel("Imaginary")
        ax_imag.set_xlabel("Receiver ID")
        ax_real.grid(True, alpha=0.25)
        ax_imag.grid(True, alpha=0.25)
        ax_real.legend(loc="best", fontsize=8)

    fig.suptitle(f"Spectral data, source {source_id}, frequency {frequency}")
    out = output_folder / f"spectral_source{source_id}_freq{frequency}.png"
    save_figure(fig, out)


def robust_abs_limit(values: np.ndarray, percentile: float = 95.0, default: float = 1.0) -> float:
    """Symmetric color limit based on a percentile of absolute finite values."""
    finite = np.abs(values[np.isfinite(values)])
    if finite.size == 0:
        return default
    limit = float(np.nanpercentile(finite, percentile))
    return limit if np.isfinite(limit) and limit > 0 else default


def set_sensible_profile_xlim(ax: plt.Axes, profiles: list[np.ndarray]) -> None:
    """Set a padded x-limit from all plotted profiles."""
    values = [p[np.isfinite(p)]
              for p in profiles if p is not None and np.any(np.isfinite(p))]
    if not values:
        return
    values = np.concatenate(values)
    xmin = float(np.nanmin(values))
    xmax = float(np.nanmax(values))
    if xmin == xmax:
        pad = 0.15 * max(1.0, abs(xmin))
    else:
        pad = 0.05 * (xmax - xmin)
    ax.set_xlim(xmin - pad, xmax + pad)


def unique_stable(values: np.ndarray) -> np.ndarray:
    """Unique values while preserving first-seen order."""
    _, index = np.unique(values, return_index=True)
    return values[np.sort(index)]


def make_sid_rid_matrix(
    sid: np.ndarray,
    rid: np.ndarray,
    values: np.ndarray,
    sids: np.ndarray,
    rids: np.ndarray,
) -> np.ndarray:
    """Create source-id by receiver-id matrix using mean where duplicates exist."""
    matrix = np.full((len(sids), len(rids)), np.nan)
    for i, source_id in enumerate(sids):
        for j, receiver_id in enumerate(rids):
            keep = (sid == source_id) & (rid == receiver_id)
            if np.any(keep):
                matrix[i, j] = np.nanmean(values[keep])
    return matrix


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    """Discover files and make selected plots."""
    setup_matplotlib()

    folder = INPUT_FOLDER.expanduser()
    output_folder = OUTPUT_FOLDER
    if not output_folder.is_absolute():
        output_folder = folder / output_folder
    output_folder.mkdir(parents=True, exist_ok=True)

    if not folder.exists():
        raise FileNotFoundError(f"INPUT_FOLDER does not exist: {folder}")

    final_models = find_final_models(folder)
    gradient_files = find_gradient_files(folder)

    print(f"Input folder: {folder}")
    print(f"Found {len(final_models)} PCBB final/inversion files.")
    print(f"Found {len(gradient_files)} gradient files.")

    if PLOT_MODEL_MAPS and final_models:
        plot_model_maps(folder, output_folder, final_models)
    if PLOT_VERTICAL_PROFILES and final_models:
        plot_vertical_profiles(folder, output_folder, final_models)
    if PLOT_GRADIENTS and gradient_files:
        plot_gradient_series(folder, output_folder, gradient_files)
    if PLOT_DIAGNOSTICS:
        plot_diagnostics(folder, output_folder)
    if PLOT_SPECTRAL_DATA:
        plot_spectral_data(folder, output_folder, SPECTRAL_SOURCE_ID, SPECTRAL_FREQUENCY)
    if PLOT_RESIDUAL_HEATMAPS:
        plot_residual_heatmaps(folder, output_folder)

    print("Done.")


if __name__ == "__main__":
    main()

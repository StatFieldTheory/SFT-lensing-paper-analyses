"""Select explicitly pinned publication displays without changing computation."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tempfile
from typing import TYPE_CHECKING, Callable, TypeAlias

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import ArrayLike, NDArray

    FloatArray: TypeAlias = NDArray[np.float64]
    Pair: TypeAlias = tuple[int, int]
    Sweep: TypeAlias = tuple[FloatArray, dict[Pair, FloatArray]]
    MainPaths: TypeAlias = tuple[Path, Path, Path]
    MainSweeps: TypeAlias = tuple[Sweep, Sweep, Sweep]
    MainDisplay: TypeAlias = tuple[
        FloatArray, dict[Pair, FloatArray], FloatArray, dict[Pair, FloatArray],
        FloatArray, dict[Pair, FloatArray]]

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
MANIFEST = HERE / "publication_displays.json"


def _manifest() -> dict:
    data = json.loads(MANIFEST.read_text())
    if (data.get("schema_version") != 1
            or data.get("policy") != "AUTHOR_RETAINED_APPROXIMATE_DISPLAY"
            or set(data.get("retained_pdf_figures", {})) != {"4", "5", "6", "17"}):
        raise ValueError("Explicit author publication-display selection required")
    active = HERE / "active_products.json"
    if (not active.is_file()
            or _digest(active) != data.get("computational_manifest_sha256")):
        raise ValueError("Publication displays do not bind the selected computational manifest")
    return data


def _digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _local_source(record: dict) -> Path:
    path = (HERE / record["path"]).resolve()
    if (not path.is_relative_to(PKG) or not path.is_file()
            or _digest(path) != record["sha256"]):
        raise ValueError(f"Publication-display source changed or is missing: {path}")
    return path


def stage_retained_pdf(figure: int, outputs: list[str]) -> str | None:
    """Archive a previous generated PDF, then stage the verified reference bytes."""
    if figure not in (4, 5, 6, 17):
        return None
    record = _manifest()["retained_pdf_figures"][str(figure)]
    if outputs != [record["output"]]:
        raise ValueError("Retained display does not match the figure's output slot")
    source = _local_source(record)
    payload = source.read_bytes()
    if hashlib.sha256(payload).hexdigest() != record["sha256"]:
        raise ValueError("Publication-display source changed while reading")
    output = (PKG / outputs[0]).resolve()
    if not output.is_relative_to(PKG) or output.suffix != ".pdf":
        raise ValueError("Publication-display output must be a package-local PDF")
    output.parent.mkdir(parents=True, exist_ok=True)
    archive = None
    # Write replacement bytes before moving the previous output into its archive.
    with tempfile.NamedTemporaryFile(dir=output.parent, suffix=".pdf", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(payload)
    if output.exists():
        date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
        archive = Path(tempfile.mkdtemp(
            prefix=f"{output.parent.name}_archive_{date}_author_display_fig{figure}_",
            dir=output.parent.parent)) / output.name
        output.rename(archive)
    temporary.replace(output)
    return (f"{record['scope']}: {source} sha256={record['sha256']} -> {output}; "
            f"previous_output_archive={archive}")


def select_main_display(
    figure: int,
    loader: Callable[[Path, int], Sweep],
    angle: Callable[[ArrayLike, ArrayLike], float],
    paths: MainPaths,
    sweeps: MainSweeps,
) -> MainDisplay:
    """Use historical diagonal pairs and the unchanged current FF/FK cross pair."""
    import numpy as np

    spec = _manifest()["main_panels"]
    if (figure not in spec["figures"] or spec["figures"] != [2, 3]
            or spec["retained_pairs"] != [[0, 0], [1, 1], [2, 2]]
            or spec["current_cross_pair"] != [0, 1]
            or spec["grid_policy"] != "EXACT_DERIVED_ANGLE_EQUALITY_NO_INTERPOLATION"):
        raise ValueError("Exact author main-panel selection required")
    keys, orders = ("order0", "ff", "fk"), (0, 2, 2)
    historical = [_local_source(spec["historical_sources"][key]) for key in keys]
    for key, path in zip(keys, paths):
        if _digest(path) != spec["current_source_sha256"][key]:
            raise ValueError(f"Current cross-panel source changed: {key}: {path}")
    old = [loader(path, order) for path, order in zip(historical, orders)]
    reference = old[0][0]

    def require_grid(grid: FloatArray, label: str) -> None:
        if grid.shape != reference.shape or not np.array_equal(grid, reference):
            maximum = (float(np.max(np.abs(grid - reference)))
                       if grid.shape == reference.shape else None)
            raise ValueError(f"Publication-display angle grid differs: {label}; "
                             f"max_abs_difference_arcmin={maximum}")

    def require_pairs(
        path: Path, order: int, grid: FloatArray,
        grouped: dict[Pair, FloatArray], pairs: tuple[Pair, ...],
    ) -> None:
        require_grid(grid, str(path))
        with np.load(path, allow_pickle=True) as data:
            a, b, o = data["a"], data["b"], data["order"]
            for pair in pairs:
                if pair not in grouped or grouped[pair].shape != reference.shape:
                    raise ValueError(f"Publication-display pair missing or malformed: {path}: {pair}")
                if not np.all(np.isfinite(grouped[pair])):
                    raise ValueError(f"Nonfinite publication-display pair: {path}: {pair}")
                indices = np.flatnonzero((a == pair[0]) & (b == pair[1]) & (o == order))
                pair_grid = np.sort(np.array([angle(data["x"][i], data["y"][i]) for i in indices]))
                require_grid(pair_grid, f"{path}: {pair}")

    diagonal = ((0, 0), (1, 1), (2, 2))
    for path, order, (grid, grouped) in zip(historical, orders, old):
        require_pairs(path, order, grid, grouped, diagonal + (((0, 1),) if order == 0 else ()))
    for path, order, (grid, grouped) in zip(paths, orders, sweeps):
        require_pairs(path, order, grid, grouped, ((0, 1),))
    display = [dict(grouped) for _, grouped in old]
    display[1][(0, 1)] = sweeps[1][1][(0, 1)]
    display[2][(0, 1)] = sweeps[2][1][(0, 1)]
    print("[publication display: historical R1 diagonal panels; current C23 FF/FK cross panel]")
    return reference, display[0], reference, display[1], reference, display[2]

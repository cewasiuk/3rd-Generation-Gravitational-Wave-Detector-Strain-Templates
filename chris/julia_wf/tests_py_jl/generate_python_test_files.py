from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT
TEST_DIR = Path(__file__).resolve().parent / "test_files" / "python"

sys.path.insert(0, str(PYTHON_DIR))

from annihilation import iso_gatom_ann_strain
from level_transition import iso_gatom_level_tr_strain


def save_real(path: Path, values) -> None:
    array = np.asarray(values, dtype=float)
    np.savetxt(path, array.reshape(-1, 1), fmt="%.18e")


def save_complex(path: Path, values) -> None:
    array = np.asarray(values, dtype=complex)
    stacked = np.column_stack([array.real, array.imag])
    np.savetxt(path, stacked, fmt="%.18e")


def save_row(path: Path, values) -> None:
    array = np.asarray(values, dtype=float).reshape(1, -1)
    np.savetxt(path, array, fmt="%.18e")


def main() -> None:
    TEST_DIR.mkdir(parents=True, exist_ok=True)

    level = iso_gatom_level_tr_strain(
        M_solar=1e-6,
        a_spin=0.999999,
        alpha=0.75,
        ne=6,
        ng=5,
        distance_kpc=1.0,
        N_g0=1e-6,
        n_time=512,
        n_fft=4096,
        n_top=120,
        verbose=False,
    )

    annihilation = iso_gatom_ann_strain(
        M_solar=3.1e-4,
        mua=2e-16,
        n=4,
        l=None,
        alpha=None,
        distance_kpc=1.0,
        iota=0.0,
        phase=0.0,
        f_min_Hz=1e7,
        f_max_Hz=1e9,
        n_f=1024,
        verbose=False,
    )

    for name, values in level.items():
        path = TEST_DIR / f"level_{name}.txt"
        if name == "lorentz_params":
            save_row(path, values)
        elif np.iscomplexobj(values):
            save_complex(path, values)
        else:
            save_real(path, values)

    for name, values in annihilation.items():
        path = TEST_DIR / f"ann_{name}.txt"
        if np.isscalar(values):
            save_real(path, [values])
        elif np.iscomplexobj(values):
            save_complex(path, values)
        else:
            save_real(path, values)


if __name__ == "__main__":
    main()

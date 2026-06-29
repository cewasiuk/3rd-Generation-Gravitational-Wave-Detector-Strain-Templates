from __future__ import annotations

from pathlib import Path
import sys
import types

import numpy as np
import scipy.special as scipy_special

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python_wf"
TEST_DIR = Path(__file__).resolve().parent / "test_files" / "python"


class _DummyMinuit:
    """Enough to import Henry's source when iminuit is unavailable."""

    def __init__(self, func, **values):
        self.func = func
        self.values = dict(values)
        self.limits = {}
        self.tol = None

    def migrad(self):
        return self


sys.modules.setdefault("iminuit", types.SimpleNamespace(Minuit=_DummyMinuit))
if not hasattr(scipy_special, "sph_harm") and hasattr(scipy_special, "sph_harm_y"):
    scipy_special.sph_harm = (
        lambda m, n, theta, phi: scipy_special.sph_harm_y(n, m, phi, theta)
    )
sys.path.insert(0, str(PYTHON_DIR))

import binary_level_transition_functions as wf  # noqa: E402


def save_real(path: Path, values) -> None:
    array = np.asarray(values, dtype=float)
    np.savetxt(path, array.reshape(-1, 1), fmt="%.18e")


def save_complex(path: Path, values) -> None:
    array = np.asarray(values, dtype=complex)
    np.savetxt(path, np.column_stack([array.real, array.imag]), fmt="%.18e")


def save_scalar(path: Path, value) -> None:
    save_real(path, [value])


def main() -> None:
    TEST_DIR.mkdir(parents=True, exist_ok=True)

    alpha = 0.42
    q = 0.4
    M_sol = 10.0
    M = M_sol * wf.Msol_in_eV
    r = 1.0e28
    m_i = 1
    m_f = -1
    n = 2
    l_i = 1
    Omega0 = wf.Omega0_binary_natural_unit(m_i, alpha, wf.G, M, n, l_i)
    eta = wf.eta_parameter(alpha, q, M_sol)
    Gamma_abs = 1.0e-20
    f = np.linspace(0.5, 1.5, 32) * wf.fc_from_Omega0(Omega0)
    h_of_iota = wf.htilde_plus(
        f, M, r, alpha, Omega0, q, m_i, m_f, eta, Gamma_abs,
        use_z_scaling=False, numerical_qc=False,
    )

    save_scalar(TEST_DIR / "eta.txt", eta)
    save_scalar(TEST_DIR / "omega0_binary.txt", Omega0)
    save_scalar(TEST_DIR / "gamma_rate.txt", wf.gamma_rate(q, M, Omega0))
    save_scalar(TEST_DIR / "q_c.txt", wf.q_c(alpha, m_i))
    save_scalar(TEST_DIR / "z_scaling.txt", wf.z_scaling_211_to_21m1(alpha, q))
    save_real(TEST_DIR / "f_grid.txt", f)
    save_real(TEST_DIR / "htilde_iota0.txt", h_of_iota(0.0))
    save_real(TEST_DIR / "htilde_iota1.txt", h_of_iota(1.0))

    mu = 1.0e-12
    mbh = 1.0e-6
    astar = 0.7
    save_scalar(TEST_DIR / "rg.txt", wf.rg(mbh))
    save_scalar(TEST_DIR / "alpha_bhsr.txt", wf.alpha(mu, mbh))
    save_scalar(TEST_DIR / "r_plus.txt", wf.r_plus(mbh, astar))
    save_scalar(TEST_DIR / "omega_hyperfine.txt", wf.omegaHyperfine(mu, mbh, astar, 2, 1, 1))
    save_scalar(TEST_DIR / "omega0_bxzh.txt", wf.omega0_bxzh(mu, mbh, 2))
    save_scalar(TEST_DIR / "omega1_bxzh.txt", wf.omega1_bxzh(mu, mbh, 2))
    save_complex(TEST_DIR / "omega_nlm_bxzh.txt", [complex(*wf.omega_nlm_bxzh(mu, mbh, astar, 2, 1, 1))])
    save_complex(TEST_DIR / "angular_ev.txt", [wf.angular_ev(0.99 * mu + 1j * 1.0e-20, mbh, astar, mu, 1, 1)])
    save_complex(TEST_DIR / "continued_fraction.txt", [
        wf.continued_fraction(0.99 * mu + 1j * 1.0e-20, mbh, astar, mu, wf.angular_ev(0.99 * mu + 1j * 1.0e-20, mbh, astar, mu, 1, 1), 1, nmax=80)
    ])


if __name__ == "__main__":
    main()

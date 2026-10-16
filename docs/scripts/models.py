"""Independent calculations used by the archive audit (raw lattice parameters)."""

from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]


def aicc(rss, n, k):
    if n <= k + 1 or rss <= 0:
        raise ValueError("AICc requires positive RSS and n > K + 1")
    return float(n * np.log(rss / n) + 2 * k * n / (n - k - 1))


def weights(values):
    values = np.asarray(values)
    w = np.exp(-0.5 * (values - values.min()))
    return w / w.sum()


def datasets():
    for source, variables in (
        ("swaminathan", ("a", "c")),
        ("julian", ("a", "c")),
        ("anthracene", ("a", "b", "c", "beta")),
    ):
        data = np.loadtxt(ROOT / "data" / f"{source}.txt")
        if source == "swaminathan":
            data = data.T
        for column, variable in enumerate(variables, 1):
            yield (
                f"{source}_{variable}",
                data[:, 0],
                data[:, column],
                10 if source == "anthracene" else 3,
            )


def polynomial(x, y, degree):
    # Scale the abscissa before least squares; avoid the archived normal equations.
    fit = Polynomial.fit(x, y, degree)
    residual = y - fit(x)
    rss = float(residual @ residual)
    n = len(y)
    return {
        "model": f"Polynomial {degree}",
        "degree": degree,
        "parameters": fit.convert().coef.tolist(),
        "rss": rss,
        "aicc": aicc(rss, n, degree + 2),
        "adjusted_r2": float(
            1 - rss / np.sum((y - y.mean()) ** 2) * (n - 1) / (n - degree - 1)
        ),
    }


def einstein(p, temperature):
    return p[0] + p[1] / np.expm1(p[2] / temperature)


def debye(p, temperature):
    z = p[2] / temperature
    numerator = (
        3953.632 - 800.6087 * z + 85.07724 * z**2 - 4.43282 * z**3 + 0.0946173 * z**4
    )
    denominator = 3953.632 + 682.0012 * z + 143.155337 * z**2 + 15.121491 * z**3 + z**4
    return p[0] + p[1] * temperature * numerator / denominator


def nonlinear(name, temperature, y):
    function = {"Einstein": einstein, "Debye": debye}[name]
    fits = []
    for theta in (80.0, 150.0, 280.0, 500.0):
        amplitude = 0.04 if name == "Einstein" else 0.0004
        result = least_squares(
            lambda p: function(p, temperature) - y,
            [float(y.min()), amplitude, theta],
            bounds=([-np.inf, 0.0, 1.0], [np.inf, np.inf, 2000.0]),
            x_scale="jac",
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            max_nfev=10000,
        )
        if not result.success:
            raise AssertionError(f"{name} optimizer failed: {result.message}")
        fits.append(result)
    scores = [aicc(2 * f.cost, len(y), 4) for f in fits]
    if np.ptp(scores) > 1e-5:
        raise AssertionError(f"{name} starts disagree: {scores}")
    fit = min(fits, key=lambda f: f.cost)
    return {
        "model": name,
        "parameters": fit.x.tolist(),
        "rss": float(2 * fit.cost),
        "aicc": aicc(2 * fit.cost, len(y), 4),
        "start_aicc_range": float(np.ptp(scores)),
    }

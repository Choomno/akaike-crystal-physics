"""Independent mathematical checks for the maintained fitting calculations."""

import numpy as np
import pytest
from models import ROOT, aicc, debye, einstein, nonlinear, polynomial, weights


@pytest.mark.parametrize("rss,n,k", [(0.02, 7, 4), (12.5, 18, 5), (3.0, 40, 2)])
def test_aicc_is_aic_plus_finite_sample_correction(rss, n, k):
    expected = n * np.log(rss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
    assert aicc(rss, n, k) == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize(
    "rss,n,k", [(0.0, 10, 3), (-1.0, 10, 3), (1.0, 4, 3), (1.0, 3, 3)]
)
def test_aicc_rejects_nonpositive_rss_or_insufficient_observations(rss, n, k):
    with pytest.raises(ValueError, match="positive RSS and n > K"):
        aicc(rss, n, k)


def test_akaike_weights_are_normalized_likelihood_ratios():
    scores = np.array([0.0, 2 * np.log(2), 2 * np.log(4)])
    actual = weights(scores)
    np.testing.assert_allclose(actual, [4 / 7, 2 / 7, 1 / 7], rtol=1e-14)
    assert actual.sum() == pytest.approx(1.0)
    assert np.all(actual >= 0)
    assert np.log(actual[0] / actual[2]) == pytest.approx(-(scores[0] - scores[2]) / 2)


def test_akaike_weights_are_invariant_under_translation_and_permutation():
    scores = np.array([-83.2, -81.7, -72.4])
    permutation = np.array([2, 0, 1])
    np.testing.assert_allclose(weights(scores + 10000), weights(scores), rtol=1e-12)
    np.testing.assert_allclose(
        weights(scores[permutation]), weights(scores)[permutation]
    )


def test_akaike_weights_handle_ties_and_large_score_differences():
    np.testing.assert_allclose(weights([10000.0] * 3), np.full(3, 1 / 3))
    actual = weights([10000.0, 12000.0, 14000.0])
    assert np.all(np.isfinite(actual))
    np.testing.assert_array_equal(actual, [1.0, 0.0, 0.0])


@pytest.mark.parametrize(
    "degree,expected_parameters,signal_sum_squares",
    [(1, [3.5, -0.015], 2.52), (2, [16.0, -0.265, 0.00125], 23.52)],
)
def test_polynomial_refit_recovers_physical_coefficients_and_fit_statistics(
    degree, expected_parameters, signal_sum_squares
):
    z = np.arange(-3.0, 4.0)
    x = 100 + 20 * z
    signal = 2 - 0.3 * z + (0.5 * z**2 if degree == 2 else 0)
    # Fourth differences annihilate polynomials through degree three, so this
    # residual is orthogonal to both candidate design spaces.
    residual = 0.01 * np.array([1.0, -4.0, 6.0, -4.0, 1.0, 0.0, 0.0])
    y = signal + residual
    fit = polynomial(x, y, degree)

    np.testing.assert_allclose(
        fit["parameters"], expected_parameters, rtol=0, atol=5e-12
    )
    fitted = np.polynomial.polynomial.polyval(x, fit["parameters"])
    np.testing.assert_allclose(fitted, signal, rtol=0, atol=1e-11)
    rss = 0.007  # 0.01**2 * (1 + 16 + 36 + 16 + 1)
    assert fit["rss"] == pytest.approx(rss, rel=1e-11)
    assert np.sum((y - fitted) ** 2) == pytest.approx(rss, rel=1e-11)
    expected_adjusted_r2 = 1 - rss / (signal_sum_squares + rss) * 6 / (6 - degree)
    assert fit["adjusted_r2"] == pytest.approx(expected_adjusted_r2, abs=1e-13)

    # K counts degree+1 regression coefficients and the fitted error variance.
    k = degree + 2
    expected_aicc = 7 * np.log(rss / 7) + 2 * k + 2 * k * (k + 1) / (6 - k)
    assert fit["aicc"] == pytest.approx(expected_aicc, abs=1e-10)


def test_einstein_model_at_exact_occupation_numbers():
    parameters = [2.0, 3.0, 100 * np.log(2)]
    # exp(theta/T)-1 is respectively 1 and 3 at these temperatures.
    np.testing.assert_allclose(
        einstein(parameters, np.array([100.0, 50.0])), [5.0, 3.0]
    )
    assert einstein(parameters, 100.0) == pytest.approx(5.0)


def test_debye_rational_model_at_unit_and_double_scaled_temperature():
    parameters = [5.0, 0.001, 100.0]
    # Literal numerator/denominator values at theta/T=1 and 2 check the
    # rational approximation independently of its polynomial evaluation.
    expected = [
        5 + 0.1 * 3233.7623373 / 4794.910028,
        5 + 0.05 * 2658.7748768 / 6027.227676,
    ]
    np.testing.assert_allclose(
        debye(parameters, np.array([100.0, 50.0])), expected, rtol=1e-14
    )
    assert debye(parameters, 100.0) == pytest.approx(expected[0], abs=1e-13)


def test_urea_debye_refit_retains_corrected_minimum_and_counts_error_variance():
    temperature, observed, _ = np.loadtxt(ROOT / "data" / "swaminathan.txt")
    fit = nonlinear("Debye", temperature, observed)
    # Guard the corrected fit, rather than merely requiring improvement over
    # the archive's early stop (AICc=-68.3644). Allow optimizer roundoff.
    assert fit["aicc"] == pytest.approx(-76.7703893575, abs=2e-4)
    residual = observed - debye(fit["parameters"], temperature)
    rss = float(residual @ residual)
    assert fit["rss"] == pytest.approx(rss, rel=1e-12)
    # Three model parameters plus the error variance: n=7, K=4.
    assert fit["aicc"] == pytest.approx(7 * np.log(rss / 7) + 28, abs=1e-10)

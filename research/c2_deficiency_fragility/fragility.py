"""Fragility and reserve-ratio utilities for the C2 deficiency pilot.

The fitted lognormal dispersion from record-to-record collapse observations is
named ``beta_rtr``.  The deficiency framework's reserve-ratio index is named
``beta_eff``.  The two quantities are deliberately not treated as identical.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, exp, isfinite, log, sqrt
from statistics import NormalDist
from typing import Iterable


_STANDARD_NORMAL = NormalDist()


def _require_positive(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _require_probability(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value) or not 0.0 < value < 1.0:
        raise ValueError(f"{name} must be strictly between zero and one")
    return value


def standard_normal_cdf(value: float) -> float:
    """Return the standard-normal cumulative distribution function."""

    return 0.5 * (1.0 + erf(float(value) / sqrt(2.0)))


def fragility_probability(intensity: float, median: float, beta_rtr: float) -> float:
    """Return lognormal collapse probability at ``intensity``.

    ``beta_rtr`` is the logarithmic dispersion used in the fitted fragility.
    """

    intensity = _require_positive("intensity", intensity)
    median = _require_positive("median", median)
    beta_rtr = _require_positive("beta_rtr", beta_rtr)
    return standard_normal_cdf(log(intensity / median) / beta_rtr)


def reserve_ratio(beta_eff: float, probability: float = 0.01) -> float:
    """Return ``IM50 / IM_p`` for a reserve-ratio index ``beta_eff``."""

    beta_eff = _require_positive("beta_eff", beta_eff)
    probability = _require_probability("probability", probability)
    if probability >= 0.5:
        raise ValueError("reserve-ratio probability must be below 0.5")
    z_value = abs(_STANDARD_NORMAL.inv_cdf(probability))
    return exp(z_value * beta_eff)


def effective_beta_from_reserve(
    reserve: float, probability: float = 0.01
) -> float:
    """Return the framework's ``beta_eff`` implied by a reserve ratio."""

    reserve = _require_positive("reserve", reserve)
    probability = _require_probability("probability", probability)
    if reserve < 1.0:
        raise ValueError("reserve must be at least one")
    if probability >= 0.5:
        raise ValueError("reserve-ratio probability must be below 0.5")
    return log(reserve) / abs(_STANDARD_NORMAL.inv_cdf(probability))


def effective_beta_after_amplification(
    reference_beta_eff: float,
    amplification: float,
    probability: float = 0.01,
) -> float:
    """Apply the framework's demand-amplification penalty to ``beta_eff``."""

    reference_beta_eff = _require_positive(
        "reference_beta_eff", reference_beta_eff
    )
    amplification = _require_positive("amplification", amplification)
    probability = _require_probability("probability", probability)
    if probability >= 0.5:
        raise ValueError("reserve-ratio probability must be below 0.5")
    return reference_beta_eff - log(amplification) / abs(
        _STANDARD_NORMAL.inv_cdf(probability)
    )


def observed_amplification(
    reference_median: float,
    case_median: float,
    *,
    reference_anchor: float | None = None,
    case_anchor: float | None = None,
) -> float:
    """Return observed reserve loss relative to a reference case.

    If anchors are omitted, this reduces to the median-capacity ratio and is
    valid only for a study that enforces a common strength/low-probability
    anchor.  Supplying both anchors computes the ratio of reserve ratios.
    """

    reference_median = _require_positive("reference_median", reference_median)
    case_median = _require_positive("case_median", case_median)
    if (reference_anchor is None) != (case_anchor is None):
        raise ValueError("supply both anchors or neither anchor")
    if reference_anchor is None:
        return reference_median / case_median
    reference_anchor = _require_positive("reference_anchor", reference_anchor)
    case_anchor = _require_positive("case_anchor", case_anchor)
    return (reference_median / reference_anchor) / (case_median / case_anchor)


@dataclass(frozen=True)
class CollapseObservation:
    """One record's collapse-intensity information.

    ``lower_im`` is a demonstrated noncollapse intensity and ``upper_im`` is a
    demonstrated collapse intensity.  Both values produce an interval-censored
    observation.  A missing upper value is right-censored; a missing lower
    value is left-censored.
    """

    record_id: str
    lower_im: float | None
    upper_im: float | None

    def __post_init__(self) -> None:
        if not self.record_id:
            raise ValueError("record_id must not be empty")
        if self.lower_im is None and self.upper_im is None:
            raise ValueError("at least one censoring bound is required")
        if self.lower_im is not None:
            _require_positive("lower_im", self.lower_im)
        if self.upper_im is not None:
            _require_positive("upper_im", self.upper_im)
        if (
            self.lower_im is not None
            and self.upper_im is not None
            and self.lower_im >= self.upper_im
        ):
            raise ValueError("lower_im must be less than upper_im")


@dataclass(frozen=True)
class CensoredFit:
    """Maximum-likelihood result for a censored lognormal fragility."""

    median: float
    beta_rtr: float
    negative_log_likelihood: float
    converged: bool
    message: str
    observation_count: int


def _log_difference(log_a: float, log_b: float) -> float:
    """Return log(exp(log_a) - exp(log_b)) for log_a > log_b."""

    if log_b >= log_a:
        return float("-inf")
    from math import expm1

    return log_a + log(-expm1(log_b - log_a))


def censored_negative_log_likelihood(
    log_median: float,
    log_beta_rtr: float,
    observations: Iterable[CollapseObservation],
) -> float:
    """Return the censored lognormal negative log likelihood.

    SciPy supplies stable normal log-CDF calculations.  It is an analysis-only
    dependency and is not required by the QuakeCore C++ runtime.
    """

    try:
        from scipy.special import log_ndtr
    except ImportError as exc:  # pragma: no cover - exercised by environment
        raise RuntimeError("SciPy is required for censored likelihood fitting") from exc

    median = exp(float(log_median))
    beta_rtr = exp(float(log_beta_rtr))
    total_log_likelihood = 0.0
    count = 0
    for observation in observations:
        count += 1
        lower_z = (
            None
            if observation.lower_im is None
            else log(observation.lower_im / median) / beta_rtr
        )
        upper_z = (
            None
            if observation.upper_im is None
            else log(observation.upper_im / median) / beta_rtr
        )

        if lower_z is None:
            term = float(log_ndtr(upper_z))
        elif upper_z is None:
            term = float(log_ndtr(-lower_z))
        else:
            term = _log_difference(
                float(log_ndtr(upper_z)), float(log_ndtr(lower_z))
            )
        if not isfinite(term):
            return float("inf")
        total_log_likelihood += term

    if count == 0:
        raise ValueError("at least one observation is required")
    return -total_log_likelihood


def fit_censored_lognormal(
    observations: Iterable[CollapseObservation],
    *,
    initial_median: float | None = None,
    initial_beta_rtr: float = 0.4,
) -> CensoredFit:
    """Fit a lognormal fragility without replacing IDA brackets by midpoints."""

    try:
        from scipy.optimize import minimize
    except ImportError as exc:  # pragma: no cover - exercised by environment
        raise RuntimeError("SciPy is required for censored fragility fitting") from exc

    observations = tuple(observations)
    if not observations:
        raise ValueError("at least one observation is required")
    finite_bounds = [
        value
        for observation in observations
        for value in (observation.lower_im, observation.upper_im)
        if value is not None
    ]
    if initial_median is None:
        sorted_bounds = sorted(finite_bounds)
        middle = len(sorted_bounds) // 2
        if len(sorted_bounds) % 2:
            initial_median = sorted_bounds[middle]
        else:
            initial_median = sqrt(
                sorted_bounds[middle - 1] * sorted_bounds[middle]
            )
    initial_median = _require_positive("initial_median", initial_median)
    initial_beta_rtr = _require_positive("initial_beta_rtr", initial_beta_rtr)

    result = minimize(
        lambda values: censored_negative_log_likelihood(
            values[0], values[1], observations
        ),
        x0=(log(initial_median), log(initial_beta_rtr)),
        method="L-BFGS-B",
        bounds=((log(min(finite_bounds) / 100.0), log(max(finite_bounds) * 100.0)),
                (log(0.02), log(2.5))),
    )
    return CensoredFit(
        median=exp(float(result.x[0])),
        beta_rtr=exp(float(result.x[1])),
        negative_log_likelihood=float(result.fun),
        converged=bool(result.success),
        message=str(result.message),
        observation_count=len(observations),
    )


if __name__ == "__main__":
    print("beta_eff  FRR_1%")
    for beta_eff in (0.20, 0.30, 0.40, 0.50, 0.60):
        print(f"{beta_eff:0.2f}      {reserve_ratio(beta_eff):0.2f}")


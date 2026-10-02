import math
import unittest

from research.c2_deficiency_fragility.fragility import (
    CollapseObservation,
    censored_negative_log_likelihood,
    effective_beta_after_amplification,
    effective_beta_from_reserve,
    fit_censored_lognormal,
    fragility_probability,
    observed_amplification,
    reserve_ratio,
)


class FrameworkCalculationTests(unittest.TestCase):
    def test_paper_reserve_table(self) -> None:
        expected = {
            0.20: 1.59,
            0.30: 2.01,
            0.40: 2.54,
            0.50: 3.20,
            0.60: 4.04,
        }
        for beta_eff, rounded_reserve in expected.items():
            self.assertAlmostEqual(
                reserve_ratio(beta_eff), rounded_reserve, delta=0.01
            )

    def test_reserve_round_trip(self) -> None:
        for beta_eff in (0.2, 0.4, 0.6):
            self.assertAlmostEqual(
                effective_beta_from_reserve(reserve_ratio(beta_eff)),
                beta_eff,
                places=12,
            )

    def test_framework_amplification_examples(self) -> None:
        examples = {
            1.00: 0.60,
            1.20: 0.52,
            1.30: 0.49,
            1.50: 0.43,
            1.56: 0.41,
            1.99: 0.30,
            2.25: 0.25,
        }
        for amplification, rounded_beta in examples.items():
            self.assertAlmostEqual(
                effective_beta_after_amplification(0.60, amplification),
                rounded_beta,
                delta=0.006,
            )

    def test_fragility_is_one_half_at_median(self) -> None:
        self.assertAlmostEqual(fragility_probability(1.3, 1.3, 0.4), 0.5)

    def test_observed_amplification_with_common_anchor(self) -> None:
        self.assertAlmostEqual(observed_amplification(2.0, 1.0), 2.0)

    def test_observed_amplification_with_different_anchors(self) -> None:
        self.assertAlmostEqual(
            observed_amplification(
                2.0, 1.2, reference_anchor=0.5, case_anchor=0.4
            ),
            4.0 / 3.0,
        )


class CensoringTests(unittest.TestCase):
    def setUp(self) -> None:
        self.observations = (
            CollapseObservation("r1", 0.8, 1.0),
            CollapseObservation("r2", 1.0, 1.2),
            CollapseObservation("r3", 1.2, None),
            CollapseObservation("r4", None, 0.9),
        )

    def test_likelihood_is_finite(self) -> None:
        value = censored_negative_log_likelihood(
            math.log(1.0), math.log(0.4), self.observations
        )
        self.assertTrue(math.isfinite(value))
        self.assertGreater(value, 0.0)

    def test_fit_preserves_observation_count(self) -> None:
        result = fit_censored_lognormal(self.observations)
        self.assertTrue(result.converged, result.message)
        self.assertEqual(result.observation_count, 4)
        self.assertGreater(result.median, 0.0)
        self.assertGreater(result.beta_rtr, 0.0)

    def test_invalid_interval_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            CollapseObservation("bad", 1.0, 1.0)

    def test_empty_fit_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            fit_censored_lognormal(())


if __name__ == "__main__":
    unittest.main()


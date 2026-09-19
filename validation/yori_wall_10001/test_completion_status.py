import unittest

from completion_status import completion_status


class CompletionStatusTests(unittest.TestCase):
    def test_all_solver_outcomes_are_exclusive(self):
        cases = [
            ("completed", "completed", 100, 100, "both_completed"),
            ("completed", "numerical_noncompletion", 100, 80, "opensees_reference_limited"),
            ("numerical_noncompletion", "completed", 80, 100, "quakecore_limited"),
            ("numerical_noncompletion", "numerical_noncompletion", 80, 80, "shared_numerical_noncompletion_same_step"),
            ("numerical_noncompletion", "numerical_noncompletion", 2374, 1981, "dual_numerical_noncompletion_different_steps"),
            ("numerical_noncompletion", "numerical_noncompletion", 4567, 7283, "dual_numerical_noncompletion_different_steps"),
        ]
        for q, o, nq, no, label in cases:
            with self.subTest(label=label, nq=nq):
                result = completion_status(q, o, nq, no)
                self.assertTrue(result[label])
                self.assertEqual(sum(result.values()), 1)

    def test_reject_unknown_outcome(self):
        with self.assertRaises(ValueError):
            completion_status("collapse", "completed", 10, 20)

    def test_reject_invalid_counts(self):
        for count in [-1, 2.5, True]:
            with self.assertRaises(ValueError):
                completion_status("completed", "completed", count, 20)


if __name__ == "__main__":
    unittest.main()

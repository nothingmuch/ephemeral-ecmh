"""Independent small-field checks for the Edwards parameter search.

Run with ``sage -python sage/test_edwards_search.py``.
"""

import unittest
from unittest.mock import patch

from edwards_search import candidate, cheap_parameters, prime_fields, random_d, small_d


def affine_count(p, d):
    return sum(
        (-x * x + y * y - 1 - d * x * x * y * y) % p == 0
        for x in range(p)
        for y in range(p)
    )


class ParameterSearchTests(unittest.TestCase):
    def test_dense_sampling_accepts_success_on_its_last_attempt(self):
        with patch("edwards_search.hashlib.shake_256") as shake:
            shake.return_value.digest.side_effect = [b"\x00"] * 143 + [b"\x02"]
            self.assertEqual(list(random_d(5, b"boundary", 1)), [2])

    def test_dense_sampling_is_reproducible_and_uses_seed_and_field(self):
        values = list(random_d(109, b"namespace one", 8))
        self.assertEqual(values, list(random_d(109, b"namespace one", 8)))
        self.assertNotEqual(values, list(random_d(109, b"namespace two", 8)))
        self.assertNotEqual(values, list(random_d(101, b"namespace one", 8)))
        self.assertEqual(len(set(values)), 8)
        for d in values:
            self.assertTrue(0 < d < 109)
            self.assertEqual(pow(d, 54, 109), 108)
        self.assertEqual(list(random_d(109, b"namespace", 0)), [])
        with self.assertRaises(ValueError):
            list(random_d(13, b"namespace", 7))

    def test_prime_enumeration_is_bounded_and_deterministic(self):
        self.assertEqual(list(prime_fields(7, 100, 3)), [(109, 19), (101, 27)])
        self.assertEqual(list(prime_fields(7, 18, 3)), [])
        self.assertEqual(list(prime_fields(7, 100, 1)), [(109, 19)])

    def test_small_parameters_are_distinct_complete_and_nonsingular(self):
        for p in (5, 13, 29, 37):
            values = list(small_d(p, p))
            self.assertEqual(len(values), (p - 1) // 2)
            self.assertEqual(len({d % p for d in values}), len(values))
            for d in values:
                self.assertNotIn(d % p, (0, p - 1))
                self.assertEqual(pow(d, (p - 1) // 2, p), p - 1)

    def test_cheap_twice_d_includes_half_integer_parameters(self):
        self.assertEqual(list(cheap_parameters(13, 1)), [(-6, "2d", 1), (6, "2d", -1)])
        values = list(cheap_parameters(13, 5))
        self.assertEqual(len({d % 13 for d, _, _ in values}), len(values))
        for d, coefficient, value in values:
            self.assertEqual((d if coefficient == "d" else 2 * d) % 13, value % 13)

    def test_exact_orders_match_exhaustive_edwards_counts(self):
        for p in (5, 13, 29, 37):
            for d in small_d(p, p):
                row = candidate(p, d, min_prime_bits=1, embedding_bound=0)
                n = affine_count(p, d)
                self.assertEqual(int(row["order"]), n)
                product = 1
                for prime, exponent in row["factors"]:
                    product *= int(prime) ** exponent
                self.assertEqual(product, n)
                self.assertEqual(int(row["quotient_order"]), n // 4)
                self.assertTrue(row["primality_proven"])

    def test_encoding_widths_do_not_claim_an_implemented_codec(self):
        row = candidate(109, 2, min_prime_bits=1, embedding_bound=0)
        self.assertEqual(row["coordinate_sign_bits"], 8)
        self.assertEqual(row["quotient_field_bits"], 7)
        self.assertEqual(row["quotient_field_bytes"], 1)
        self.assertEqual(
            row["quotient_codec"], "Decaf capacity; codec not implemented here"
        )

    def test_large_factor_and_embedding_filters_are_reported(self):
        row = candidate(13, 2, min_prime_bits=112, embedding_bound=100)
        self.assertFalse(row["selected"])
        self.assertIn("prime_factor_below_target", row["rejections"])
        if int(row["largest_odd_prime"]) > 1:
            self.assertIn("embedding_degree_within_bound", row["rejections"])

    def test_invalid_domains_and_curves_are_rejected(self):
        for args in ((2, 8, 1), (128, 0, 1), (128, 275, 0)):
            with self.assertRaises(ValueError):
                list(prime_fields(*args))
        for p, d in ((15, 2), (7, 3), (13, 0), (13, -1), (13, 1)):
            with self.assertRaises(ValueError):
                candidate(p, d, min_prime_bits=1, embedding_bound=0)


if __name__ == "__main__":
    unittest.main()

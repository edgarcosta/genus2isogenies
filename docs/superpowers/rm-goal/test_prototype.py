#!/usr/bin/env python3
"""Public-boundary tests for the maximal-RM arithmetic prototype."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import unittest


PROTOTYPE = Path(__file__).with_name("prototype.sage")
SAGE = shutil.which("sage")
BRUTE_FORCE_MAX_DISCRIMINANT = 80


def run_cli(*arguments):
    """Run only the frozen public process seam."""
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = "0"
    return subprocess.run(
        [SAGE, "-python", str(PROTOTYPE), *map(str, arguments)],
        check=False,
        capture_output=True,
        env=env,
        text=True,
    )


def run_enumerate(discriminant, conductor=1):
    """Run the frozen enumerate command and decode its JSON document."""
    completed = run_cli(
        "enumerate",
        "--field-discriminant",
        discriminant,
        "--order-conductor",
        conductor,
    )
    return completed, json.loads(completed.stdout)


def prime_divisors(integer):
    """Independent elementary factorization for literal count expectations."""
    divisors = []
    candidate = 2
    while candidate * candidate <= integer:
        if integer % candidate == 0:
            divisors.append(candidate)
            while integer % candidate == 0:
                integer //= candidate
        candidate += 1
    if integer > 1:
        divisors.append(integer)
    return divisors


def is_squarefree(integer):
    return all(integer % (prime * prime) for prime in prime_divisors(integer))


def is_positive_fundamental_discriminant(integer):
    if integer <= 1:
        return False
    if integer % 4 == 1:
        return is_squarefree(integer)
    if integer % 4 == 0:
        radicand = integer // 4
        return radicand % 4 in (2, 3) and is_squarefree(radicand)
    return False


class PrototypeCLITests(unittest.TestCase):
    def test_D5_emits_the_literal_identity_arithmetic_request(self):
        completed, document = run_enumerate(5)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["schema"], "rm-hm-arithmetic-v1")
        self.assertEqual(document["status"], "OK_ARITHMETIC_ONLY")
        self.assertEqual(document["surface_scope"], "NOT_EVALUATED")
        self.assertEqual(
            document["theorem_status"],
            "GEOMETRIC_INTERPRETATION_CONDITIONAL",
        )
        self.assertEqual(document["field"]["defining_polynomial"], [-1, -1, 1])
        self.assertEqual(document["field"]["discriminant"], 5)
        self.assertEqual(document["order"]["basis"], [[1, 0], [0, 1]])
        self.assertEqual(document["ordinary_class_group"]["invariants"], [])
        self.assertEqual(document["narrow_class_group"]["invariants"], [])
        self.assertEqual(document["unit_signature_image"], [[-1, -1], [-1, 1], [1, -1], [1, 1]])
        self.assertEqual(
            document["polarization_unit_classes"],
            [{"element": [1, 0], "id": "unit-0", "signatures": [1, 1]}],
        )
        self.assertEqual(len(document["hm_ideal_classes"]), 1)
        self.assertEqual(len(document["hm_records"]), 1)
        self.assertEqual(document["hm_records"][0]["alpha"], [1, 0])

    def test_D12_emits_both_positive_unit_square_classes(self):
        completed, document = run_enumerate(12)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["field"]["defining_polynomial"], [-3, 0, 1])
        self.assertEqual(document["ordinary_class_group"]["invariants"], [])
        self.assertEqual(document["narrow_class_group"]["invariants"], [2])
        self.assertEqual(document["unit_signature_image"], [[-1, -1], [1, 1]])
        self.assertEqual(
            document["polarization_unit_classes"],
            [
                {"element": [1, 0], "id": "unit-0", "signatures": [1, 1]},
                {"element": [2, 1], "id": "unit-1", "signatures": [1, 1]},
            ],
        )
        self.assertEqual(
            [record["alpha"] for record in document["hm_records"]],
            [[1, 0], [2, 1]],
        )

    def test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq(self):
        completed, document = run_enumerate(40)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["ordinary_class_group"]["invariants"], [2])
        self.assertEqual(document["narrow_class_group"]["invariants"], [2])
        self.assertEqual(
            document["hm_ideal_classes"],
            [
                {
                    "id": "ideal-0",
                    "norm": 1,
                    "ordinary_class_coordinates": [0],
                    "row_hnf": [[1, 0], [0, 1]],
                },
                {
                    "id": "ideal-1",
                    "norm": 2,
                    "ordinary_class_coordinates": [1],
                    "row_hnf": [[2, 0], [0, 1]],
                },
            ],
        )
        self.assertEqual(
            [record["alpha"] for record in document["hm_records"]],
            [[1, 0], [2, 0]],
        )
        self.assertEqual(
            [record["expected_quotient_degree"] for record in document["hm_records"]],
            [1, 4],
        )
        self.assertTrue(
            all(all(record["checks"].values()) for record in document["hm_records"])
        )

    def test_D60_emits_the_literal_cartesian_HM_representatives(self):
        completed, document = run_enumerate(60)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["ordinary_class_group"]["invariants"], [2])
        self.assertEqual(document["narrow_class_group"]["invariants"], [2, 2])
        self.assertEqual(
            [record["element"] for record in document["polarization_unit_classes"]],
            [[1, 0], [4, 1]],
        )
        self.assertEqual(
            document["hm_ideal_classes"],
            [
                {
                    "id": "ideal-0",
                    "norm": 1,
                    "ordinary_class_coordinates": [0],
                    "row_hnf": [[1, 0], [0, 1]],
                },
                {
                    "id": "ideal-1",
                    "norm": 2,
                    "ordinary_class_coordinates": [1],
                    "row_hnf": [[1, 1], [0, 2]],
                },
            ],
        )
        self.assertEqual(
            [record["alpha"] for record in document["hm_records"]],
            [[1, 0], [4, 1], [2, 0], [8, 2]],
        )
        self.assertEqual(
            document["hm_structure"],
            {
                "eligible_ordinary_class_count": 2,
                "ideal_class_condition": "ker(Sq:Cl(F)->Cl+(F))",
                "positive_unit_class_count": 2,
                "record_count": 4,
            },
        )

    def test_D145_selects_only_the_order_two_class_from_ordinary_C4(self):
        completed, document = run_enumerate(145)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["ordinary_class_group"]["invariants"], [4])
        self.assertEqual(document["narrow_class_group"]["invariants"], [4])
        self.assertEqual(document["minkowski_enumeration"]["bound"], 6)
        self.assertEqual(document["minkowski_enumeration"]["ordinary_classes_covered"], 4)
        self.assertEqual(
            document["hm_ideal_classes"],
            [
                {
                    "id": "ideal-0",
                    "norm": 1,
                    "ordinary_class_coordinates": [0],
                    "row_hnf": [[1, 0], [0, 1]],
                },
                {
                    "id": "ideal-1",
                    "norm": 4,
                    "ordinary_class_coordinates": [2],
                    "row_hnf": [[1, 3], [0, 4]],
                },
            ],
        )
        self.assertEqual(
            [record["alpha"] for record in document["hm_records"]],
            [[1, 0], [17, 3]],
        )

    def test_D205_rejects_the_ineligible_ordinary_two_class(self):
        completed, document = run_enumerate(205)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["ordinary_class_group"]["invariants"], [2])
        self.assertEqual(document["narrow_class_group"]["invariants"], [4])
        self.assertEqual(document["minkowski_enumeration"]["ordinary_classes_covered"], 2)
        self.assertEqual(document["hm_structure"]["eligible_ordinary_class_count"], 1)
        self.assertEqual(
            [record["element"] for record in document["polarization_unit_classes"]],
            [[1, 0], [20, 3]],
        )
        self.assertEqual(
            document["hm_ideal_classes"],
            [
                {
                    "id": "ideal-0",
                    "norm": 1,
                    "ordinary_class_coordinates": [0],
                    "row_hnf": [[1, 0], [0, 1]],
                }
            ],
        )
        self.assertEqual(
            [record["alpha"] for record in document["hm_records"]],
            [[1, 0], [20, 3]],
        )

    def test_D229_uses_the_kernel_not_the_surjective_image_of_squaring(self):
        completed, document = run_enumerate(229)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(document["ordinary_class_group"]["invariants"], [3])
        self.assertEqual(document["narrow_class_group"]["invariants"], [3])
        self.assertEqual(document["minkowski_enumeration"]["ordinary_classes_covered"], 3)
        self.assertEqual(document["hm_structure"]["ideal_class_condition"], "ker(Sq:Cl(F)->Cl+(F))")
        self.assertEqual(document["hm_structure"]["eligible_ordinary_class_count"], 1)
        self.assertEqual(document["hm_structure"]["record_count"], 1)
        self.assertEqual(
            document["hm_ideal_classes"],
            [
                {
                    "id": "ideal-0",
                    "norm": 1,
                    "ordinary_class_coordinates": [0],
                    "row_hnf": [[1, 0], [0, 1]],
                }
            ],
        )

    def test_nonmaximal_order_is_a_validated_rejection_without_HM_data(self):
        completed, document = run_enumerate(5, conductor=2)

        self.assertEqual(completed.returncode, 2, completed.stderr)
        self.assertEqual(document["status"], "UNSUPPORTED_NONMAXIMAL_ORDER")
        self.assertEqual(document["field"]["discriminant"], 5)
        self.assertEqual(
            document["order"],
            {
                "basis": [[1, 0], [0, 2]],
                "conductor": 2,
                "discriminant": 20,
                "index_in_maximal_order": 2,
                "is_maximal": False,
            },
        )
        self.assertNotIn("hm_ideal_classes", document)
        self.assertNotIn("hm_records", document)

    def test_order_discriminant_is_not_accepted_as_a_field_discriminant(self):
        completed, document = run_enumerate(20)

        self.assertEqual(completed.returncode, 2, completed.stderr)
        self.assertEqual(document["status"], "INVALID_FIELD_DISCRIMINANT")
        self.assertEqual(
            document["input"],
            {"field_discriminant": 20, "order_conductor": 1},
        )
        self.assertNotIn("field", document)
        self.assertNotIn("hm_records", document)

    def test_success_stdout_is_byte_stable(self):
        first, first_document = run_enumerate(60)
        second, second_document = run_enumerate(60)

        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(first.stderr, "")
        self.assertEqual(second.stderr, "")
        self.assertEqual(first.stdout, second.stdout)
        self.assertEqual(first_document, second_document)
        self.assertTrue(first.stdout.endswith("\n"))
        self.assertEqual(first.stdout.count("\n"), 1)

    def test_success_document_exposes_every_mandatory_run_check(self):
        completed, document = run_enumerate(60)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            sorted(document["checks"]),
            [
                "all_record_checks",
                "canonical_field_basis",
                "canonical_polynomial_discriminant",
                "eligible_ideals_pass_direct_kernel_of_Sq",
                "field_discriminant",
                "hm_pairs_distinct_under_equivalence",
                "identity_included",
                "maximal_order_discriminant",
                "minkowski_class_coverage",
                "order_equals_maximal_order",
                "proof_mode",
                "record_count_is_cartesian_product",
                "selected_ideals_integral_invertible",
                "unit_representatives_distinct_modulo_squares",
                "unit_representatives_totally_positive",
            ],
        )
        self.assertTrue(all(document["checks"].values()))

    def test_invalid_CLI_arguments_have_an_explicit_JSON_status(self):
        completed = run_cli(
            "enumerate", "--field-discriminant", "not-an-integer"
        )
        document = json.loads(completed.stdout)

        self.assertEqual(completed.returncode, 2)
        self.assertEqual(document["status"], "INVALID_ARGUMENTS")
        self.assertNotIn("hm_records", document)
        self.assertEqual(completed.stdout.count("\n"), 1)

    def test_bounded_bruteforce_counts_match_genus_theory_and_narrow_two_torsion(self):
        discriminants = [
            D
            for D in range(2, BRUTE_FORCE_MAX_DISCRIMINANT + 1)
            if is_positive_fundamental_discriminant(D)
        ]
        self.assertEqual(
            discriminants,
            [5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 53, 56, 57, 60, 61, 65, 69, 73, 76, 77],
        )

        for D in discriminants:
            with self.subTest(field_discriminant=D):
                completed, document = run_enumerate(D)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                genus_theory_count = 2 ** (len(prime_divisors(D)) - 1)
                narrow_two_torsion_count = 1
                for invariant in document["narrow_class_group"]["invariants"]:
                    narrow_two_torsion_count *= 2 if invariant % 2 == 0 else 1
                self.assertEqual(len(document["hm_records"]), genus_theory_count)
                self.assertEqual(len(document["hm_records"]), narrow_two_torsion_count)
                self.assertTrue(
                    all(
                        all(record["checks"].values())
                        for record in document["hm_records"]
                    )
                )


if __name__ == "__main__":
    unittest.main()

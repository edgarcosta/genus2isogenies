#!/usr/bin/env sage -python
"""Black-box validation of the maximal real-quadratic HM prototype.

This script deliberately invokes ``prototype.sage`` only through its frozen
CLI.  It never imports implementation functions.  Every JSON arithmetic pair
is reconstructed in fresh proof-enabled Sage objects before it is accepted.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from sage.all import AA, Matrix, NumberField, ZZ, polygen, proof


PROTOTYPE = Path(__file__).with_name("prototype.sage")
DIAGNOSTIC_DISCRIMINANTS = (5, 12, 40, 60, 145, 205, 229)


class ValidationFailure(AssertionError):
    """A black-box result failed an independent check."""


def require(condition, message):
    if not condition:
        raise ValidationFailure(message)


def prime_divisors(integer):
    """Elementary trial division, independent of Sage's factorization API."""
    integer = int(integer)
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
    """Elementary squarefreeness test for the bounded input range."""
    integer = int(integer)
    candidate = 2
    while candidate * candidate <= integer:
        if integer % (candidate * candidate) == 0:
            return False
        candidate += 1
    return True


def is_positive_fundamental_discriminant(integer):
    """Classify positive quadratic field discriminants elementarily."""
    integer = int(integer)
    if integer <= 1:
        return False
    if integer % 4 == 1:
        return is_squarefree(integer)
    if integer % 4 == 0:
        radicand = integer // 4
        return radicand % 4 in (2, 3) and is_squarefree(radicand)
    return False


def canonical_field(discriminant):
    """Recreate the public [1,w] presentation without prototype internals."""
    x = polygen(ZZ, "x")
    polynomial = (
        x**2 - x + (1 - discriminant) // 4
        if discriminant % 2
        else x**2 - discriminant // 4
    )
    return polynomial, NumberField(polynomial, "v")


def embeddings(field):
    result = list(field.embeddings(AA))
    require(len(result) == 2, "field is not real quadratic")
    return sorted(result, key=lambda embedding: embedding(field.gen()))


def signatures(element, real_embeddings):
    values = [embedding(element) for embedding in real_embeddings]
    require(all(value != 0 for value in values), "unexpected zero real embedding")
    return tuple(1 if value > 0 else -1 for value in values)


def unit_parity_representatives(unit_group):
    invariants = tuple(ZZ(value) for value in unit_group.invariants())
    require(invariants == (2, 0), "unexpected real-quadratic unit group")
    torsion, free = [generator.value() for generator in unit_group.gens()]
    return [
        torsion**torsion_exponent * free**free_exponent
        for torsion_exponent in (0, 1)
        for free_exponent in (0, 1)
    ]


def positive_unit_quotient_size(unit_group, real_embeddings):
    signature_image = {
        signatures(unit, real_embeddings)
        for unit in unit_parity_representatives(unit_group)
    }
    require(4 % len(signature_image) == 0, "invalid unit signature image")
    return 4 // len(signature_image), signature_image


def narrow_principal(ideal, parity_units, real_embeddings):
    """Test narrow principality directly by generator signs."""
    if not ideal.is_principal(proof=True):
        return False
    generators = ideal.gens_reduced(proof=True)
    require(len(generators) == 1, "principal ideal has no unique reduced generator")
    generator = generators[0]
    return any(
        signatures(unit * generator, real_embeddings) == (1, 1)
        for unit in parity_units
    )


def direct_sq_data(class_group, parity_units, real_embeddings):
    """Compute kernel and image size of Sq: Cl -> Cl+ from ideal signs."""
    squared_ideals = [element.ideal() ** 2 for element in class_group]
    kernel_size = sum(
        narrow_principal(ideal, parity_units, real_embeddings)
        for ideal in squared_ideals
    )

    image_representatives = []
    for squared_ideal in squared_ideals:
        is_new = True
        for representative in image_representatives:
            quotient = squared_ideal * (~representative)
            if narrow_principal(quotient, parity_units, real_embeddings):
                is_new = False
                break
        if is_new:
            image_representatives.append(squared_ideal)
    require(
        kernel_size * len(image_representatives) == int(class_group.order()),
        "direct Sq kernel/image sizes violate the first isomorphism theorem",
    )
    return kernel_size, len(image_representatives)


def run_public_cli(sage_command, discriminant):
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = "0"
    command = [
        sage_command,
        "-python",
        str(PROTOTYPE),
        "enumerate",
        "--field-discriminant",
        str(discriminant),
        "--order-conductor",
        "1",
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        env=env,
        text=True,
    )
    require(completed.returncode == 0, "D={} CLI exit {}: {}".format(
        discriminant, completed.returncode, completed.stderr.strip()
    ))
    require(completed.stderr == "", "D={} emitted stderr".format(discriminant))
    require(
        completed.stdout.count("\n") == 1 and completed.stdout.endswith("\n"),
        "D={} did not emit exactly one JSON line".format(discriminant),
    )
    try:
        document = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValidationFailure("D={} emitted invalid JSON".format(discriminant)) from error
    require(document.get("status") == "OK_ARITHMETIC_ONLY", "D={} bad status".format(discriminant))
    return document


def element_from_coefficients(field, coefficients):
    require(len(coefficients) == 2, "element does not have two [1,w] coefficients")
    return ZZ(coefficients[0]) + ZZ(coefficients[1]) * field.gen()


def ideal_from_hnf(field, rows):
    require(
        len(rows) == 2 and all(len(row) == 2 for row in rows),
        "ideal row-HNF is not 2 by 2",
    )
    matrix = Matrix(ZZ, rows)
    require(matrix.det() > 0, "ideal row-HNF has nonpositive determinant")
    require(matrix.hermite_form() == matrix, "reported ideal matrix is not row-HNF")
    generators = [
        ZZ(row[0]) + ZZ(row[1]) * field.gen()
        for row in rows
    ]
    return field.ideal(generators), matrix


def unit_is_square(unit_group, unit):
    coordinates = [ZZ(value) for value in unit_group(unit).list()]
    return all(value % 2 == 0 for value in coordinates)


def pairs_equivalent(left, right, unit_group):
    """Exact test for (a,alpha)~(x*a,x^2*alpha)."""
    left_ideal, left_alpha = left
    right_ideal, right_alpha = right
    quotient_ideal = right_ideal * (~left_ideal)
    if not quotient_ideal.is_principal(proof=True):
        return False
    generators = quotient_ideal.gens_reduced(proof=True)
    require(len(generators) == 1, "principal ideal quotient has no generator")
    beta = generators[0]
    residual_unit = right_alpha / (beta**2 * left_alpha)
    require(
        residual_unit.is_integral() and abs(ZZ(residual_unit.norm())) == 1,
        "equivalent ideal classes produced a nonunit residual",
    )
    return unit_is_square(unit_group, residual_unit)


def reconstruct_document(discriminant, document):
    polynomial, field = canonical_field(discriminant)
    real_embeddings = embeddings(field)
    class_group = field.class_group(proof=True)
    narrow_class_group = field.narrow_class_group(proof=True)
    unit_group = field.unit_group(proof=True)
    parity_units = unit_parity_representatives(unit_group)

    require(int(field.discriminant()) == discriminant, "wrong reconstructed field discriminant")
    require(
        document["field"]["defining_polynomial"] == [int(value) for value in polynomial.list()],
        "D={} defining polynomial mismatch".format(discriminant),
    )
    require(document["order"]["is_maximal"] is True, "D={} order is not maximal".format(discriminant))
    require(
        document["ordinary_class_group"]["invariants"]
        == [int(value) for value in class_group.invariants()],
        "D={} ordinary class invariants mismatch".format(discriminant),
    )
    require(
        document["narrow_class_group"]["invariants"]
        == [int(value) for value in narrow_class_group.invariants()],
        "D={} narrow class invariants mismatch".format(discriminant),
    )

    positive_unit_count, signature_image = positive_unit_quotient_size(
        unit_group, real_embeddings
    )
    reported_signature_image = {
        tuple(value) for value in document["unit_signature_image"]
    }
    require(reported_signature_image == signature_image, "D={} unit signatures mismatch".format(discriminant))
    unit_records = document["polarization_unit_classes"]
    require(len(unit_records) == positive_unit_count, "D={} positive-unit count mismatch".format(discriminant))
    require(len({record["id"] for record in unit_records}) == len(unit_records), "duplicate unit ID")
    exact_units = []
    for record in unit_records:
        unit = element_from_coefficients(field, record["element"])
        require(field.ideal(unit) == field.ideal(1), "reported unit is not a unit")
        require(signatures(unit, real_embeddings) == (1, 1), "reported unit is not totally positive")
        require(record["signatures"] == [1, 1], "reported unit signatures are false")
        exact_units.append(unit)
    for left_index, left in enumerate(exact_units):
        for right in exact_units[left_index + 1 :]:
            require(not unit_is_square(unit_group, left / right), "duplicate positive-unit square classes")

    ideal_records = document["hm_ideal_classes"]
    require(len({record["id"] for record in ideal_records}) == len(ideal_records), "duplicate ideal ID")
    ideal_records_by_id = {record["id"]: record for record in ideal_records}
    ideals = {}
    ordinary_coordinates = set()
    for record in ideal_records:
        ideal, hnf = ideal_from_hnf(field, record["row_hnf"])
        require(ideal.is_integral(), "reported ideal is not integral")
        require(ideal * (~ideal) == field.ideal(1), "reported ideal is not invertible")
        require(int(ideal.norm()) == record["norm"], "reported ideal norm mismatch")
        require(abs(int(hnf.det())) == record["norm"], "row-HNF determinant/norm mismatch")
        coordinates = tuple(int(value) for value in class_group(ideal).list())
        require(
            list(coordinates) == record["ordinary_class_coordinates"],
            "reported ordinary class coordinates mismatch",
        )
        require(coordinates not in ordinary_coordinates, "duplicate reported ordinary ideal class")
        ordinary_coordinates.add(coordinates)
        ideals[record["id"]] = ideal

    direct_kernel_size, direct_image_size = direct_sq_data(
        class_group, parity_units, real_embeddings
    )
    require(len(ideals) == direct_kernel_size, "D={} incomplete Sq kernel".format(discriminant))

    records = document["hm_records"]
    require(len({record["id"] for record in records}) == len(records), "duplicate HM record ID")
    require(
        len({(record["ideal_id"], record["unit_class_id"]) for record in records})
        == len(records),
        "duplicate ideal/unit record coordinates",
    )
    exact_pairs = []
    for record in records:
        require(record["ideal_id"] in ideals, "HM record references unknown ideal")
        require(record["unit_class_id"] in {item["id"] for item in unit_records}, "HM record references unknown unit")
        require(
            record["ideal"] == ideal_records_by_id[record["ideal_id"]],
            "HM record-local ideal data do not match its ideal ID",
        )
        ideal = ideals[record["ideal_id"]]
        alpha = element_from_coefficients(field, record["alpha"])
        norm = ZZ(ideal.norm())
        actual_signatures = signatures(alpha, real_embeddings)
        require(ideal**2 == field.ideal(alpha), "HM record fails a^2=(alpha)")
        require(actual_signatures == (1, 1), "HM alpha is not totally positive")
        require(ZZ(alpha.norm()) == norm**2, "HM alpha norm does not equal ideal norm squared")
        require(record["alpha_signatures"] == [1, 1], "HM alpha signature claim mismatch")
        require(record["expected_graph_weight"] == int(norm), "HM graph-weight norm mismatch")
        require(record["expected_quotient_degree"] == int(norm**2), "HM degree norm mismatch")
        exact_pairs.append((ideal, alpha))

    for left_index, left in enumerate(exact_pairs):
        for right in exact_pairs[left_index + 1 :]:
            require(not pairs_equivalent(left, right, unit_group), "duplicate HM arithmetic pairs")

    genus_count = 2 ** (len(prime_divisors(discriminant)) - 1)
    narrow_two_torsion_count = 1
    for invariant in narrow_class_group.invariants():
        narrow_two_torsion_count *= 2 if ZZ(invariant) % 2 == 0 else 1
    require(len(records) == genus_count, "D={} genus-theory count mismatch".format(discriminant))
    require(len(records) == narrow_two_torsion_count, "D={} narrow 2-torsion mismatch".format(discriminant))
    require(
        len(records) == direct_kernel_size * positive_unit_count,
        "D={} HM exact-sequence cardinality mismatch".format(discriminant),
    )
    return {
        "class_invariants": [int(value) for value in class_group.invariants()],
        "direct_sq_image_size": direct_image_size,
        "direct_sq_kernel_size": direct_kernel_size,
        "hm_record_count": len(records),
        "narrow_class_invariants": [int(value) for value in narrow_class_group.invariants()],
        "narrow_two_torsion_count": narrow_two_torsion_count,
        "positive_unit_class_count": positive_unit_count,
        "signature_image_size": len(signature_image),
    }


def assert_diagnostic_coverage(results):
    require(results[5]["signature_image_size"] == 4, "D=5 must have the full unit-signature image")
    require(results[12]["signature_image_size"] == 2, "D=12 must have the same-sign unit regime")
    require(results[12]["class_invariants"] == [], "D=12 ordinary class group changed")
    require(results[12]["narrow_class_invariants"] == [2], "D=12 narrow class group changed")
    require(results[40]["class_invariants"] == [2], "D=40 must have ordinary 2-torsion")
    require(results[40]["direct_sq_kernel_size"] == 2, "D=40 nontrivial 2-class must be admissible")
    require(results[60]["hm_record_count"] == 4, "D=60 must combine two ideals and two units")
    require(results[145]["class_invariants"] == [4], "D=145 must exhibit ordinary C4")
    require(results[145]["direct_sq_kernel_size"] == 2, "D=145 Sq kernel must have order two")
    require(results[205]["class_invariants"] == [2], "D=205 must exhibit ordinary C2")
    require(results[205]["narrow_class_invariants"] == [4], "D=205 narrow group must be C4")
    require(results[205]["direct_sq_kernel_size"] == 1, "D=205 ordinary 2-class must be ineligible")
    require(results[229]["class_invariants"] == [3], "D=229 must exhibit ordinary C3")
    require(results[229]["direct_sq_image_size"] == 3, "D=229 Sq image must have order three")
    require(results[229]["direct_sq_kernel_size"] == 1, "D=229 Sq kernel must be trivial")


def magma_class_invariants(magma_command, discriminants):
    ds = ",".join(str(value) for value in discriminants)
    source = "\n".join(
        [
            "for D in [{}] do".format(ds),
            "  K<a> := QuadraticField(D);",
            "  O := Integers(K);",
            "  C, cmap := ClassGroup(O);",
            "  N, nmap := NarrowClassGroup(O);",
            "  cs := Join([IntegerToString(Integers()!x) : x in Invariants(C)], \",\");",
            "  ns := Join([IntegerToString(Integers()!x) : x in Invariants(N)], \",\");",
            "  printf \"MAGMA|%o|%o|%o\\n\", D, cs, ns;",
            "end for;",
            "quit;",
        ]
    )
    completed = subprocess.run(
        [magma_command, "-b"],
        input=source + "\n",
        check=False,
        capture_output=True,
        text=True,
    )
    require(completed.returncode == 0, "Magma class-group run failed: {}".format(completed.stderr.strip()))
    parsed = {}
    for line in completed.stdout.splitlines():
        if not line.startswith("MAGMA|"):
            continue
        unused, d_text, ordinary_text, narrow_text = line.split("|", 3)
        parsed[int(d_text)] = {
            "class_invariants": [] if ordinary_text == "" else [int(value) for value in ordinary_text.split(",")],
            "narrow_class_invariants": [] if narrow_text == "" else [int(value) for value in narrow_text.split(",")],
        }
    require(set(parsed) == set(discriminants), "Magma did not report every diagnostic field")
    return parsed


def parse_arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-discriminant", type=int, default=260)
    parser.add_argument(
        "--require-magma",
        action="store_true",
        help="fail instead of reporting a skipped independent Magma cross-check",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_arguments(argv)
    proof.all(True)
    require(args.max_discriminant >= max(DIAGNOSTIC_DISCRIMINANTS), "range omits a required diagnostic field")
    sage_command = shutil.which("sage")
    require(sage_command is not None, "sage command not found")
    require(PROTOTYPE.is_file(), "prototype CLI not found")

    discriminants = [
        value
        for value in range(2, args.max_discriminant + 1)
        if is_positive_fundamental_discriminant(value)
    ]
    require(
        all(ZZ(value).is_fundamental_discriminant() for value in discriminants),
        "elementary and Sage fundamental-discriminant classifiers disagree",
    )
    results = {}
    for discriminant in discriminants:
        document = run_public_cli(sage_command, discriminant)
        results[discriminant] = reconstruct_document(discriminant, document)
    assert_diagnostic_coverage(results)

    magma_command = shutil.which("magma")
    magma_status = "SKIPPED_NOT_INSTALLED"
    if magma_command is not None:
        magma_results = magma_class_invariants(magma_command, DIAGNOSTIC_DISCRIMINANTS)
        for discriminant, magma_result in magma_results.items():
            require(
                magma_result["class_invariants"] == results[discriminant]["class_invariants"],
                "D={} Magma ordinary class invariants mismatch".format(discriminant),
            )
            require(
                magma_result["narrow_class_invariants"] == results[discriminant]["narrow_class_invariants"],
                "D={} Magma narrow class invariants mismatch".format(discriminant),
            )
        magma_status = "PASS_7_DIAGNOSTIC_FIELDS"
    elif args.require_magma:
        raise ValidationFailure("--require-magma was set but magma was not found")

    summary = {
        "diagnostic_fields": {
            str(discriminant): results[discriminant]
            for discriminant in DIAGNOSTIC_DISCRIMINANTS
        },
        "fundamental_discriminants_checked": len(discriminants),
        "magma_cross_check": magma_status,
        "max_discriminant": args.max_discriminant,
        "proof_mode": True,
        "prototype_access": "PUBLIC_CLI_ONLY",
        "reconstructed_hm_records": sum(
            result["hm_record_count"] for result in results.values()
        ),
        "status": "PASS",
    }
    sys.stdout.write(json.dumps(summary, separators=(",", ":"), sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        sys.stderr.write("{}: {}\n".format(type(error).__name__, error))
        sys.exit(1)

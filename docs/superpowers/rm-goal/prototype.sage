#!/usr/bin/env sage -python
"""Exact maximal real-quadratic Hurwitz--Maass arithmetic enumerator.

Public interface::

    PYTHONHASHSEED=0 sage -python prototype.sage enumerate \
        --field-discriminant D --order-conductor f

The output is arithmetic data only.  No abelian surface or quotient is
constructed or deduplicated here.
"""

import argparse
import json
import sys

from sage.all import AA, Matrix, NumberField, ZZ, polygen, proof


SCHEMA = "rm-hm-arithmetic-v1"
SUCCESS = "OK_ARITHMETIC_ONLY"
SURFACE_SCOPE = "NOT_EVALUATED"
THEOREM_STATUS = "GEOMETRIC_INTERPRETATION_CONDITIONAL"


class ArithmeticInconsistency(RuntimeError):
    """A mandatory exact-arithmetic check failed."""


class CLIValidationError(ValueError):
    """The command line did not match the public grammar."""


class ProofComputationFailure(RuntimeError):
    """A required proof-enabled class/unit computation did not certify."""


class PositiveGeneratorFailure(RuntimeError):
    """A principal-square/sign computation did not certify."""


class JSONArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise CLIValidationError(message)


def _emit(document):
    """Write exactly one canonical JSON document to stdout."""
    sys.stdout.write(
        json.dumps(document, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    )


def _integer_coefficients(element):
    coefficients = list(element)
    if len(coefficients) < 2:
        coefficients += [ZZ(0)] * (2 - len(coefficients))
    if any(coefficient.denominator() != 1 for coefficient in coefficients):
        raise ArithmeticInconsistency("an asserted integral element has nonintegral coordinates")
    return [int(coefficient) for coefficient in coefficients]


def _embeddings(K):
    embeddings = list(K.embeddings(AA))
    if len(embeddings) != 2:
        raise ArithmeticInconsistency("the field does not have two real embeddings")
    return sorted(embeddings, key=lambda embedding: embedding(K.gen()))


def _signatures(element, embeddings):
    values = [embedding(element) for embedding in embeddings]
    if any(value == 0 for value in values):
        raise ArithmeticInconsistency("a unit or positive generator has a zero embedding")
    return [1 if value > 0 else -1 for value in values]


def _row_hnf(ideal):
    rows = []
    for basis_element in ideal.basis():
        coordinates = list(basis_element)
        if any(coordinate.denominator() != 1 for coordinate in coordinates):
            raise ArithmeticInconsistency("an integral ideal has a nonintegral lattice basis")
        rows.append([ZZ(coordinate) for coordinate in coordinates])
    hnf = Matrix(ZZ, rows).hermite_form()
    return [[int(entry) for entry in row] for row in hnf.rows()]


def _canonical_positive_modulo_square_units(element, free_unit, embeddings):
    """Normalize a positive element in the fixed pair of real embeddings."""
    if _signatures(element, embeddings) != [1, 1]:
        raise ArithmeticInconsistency("attempted to normalize a nonpositive element")
    eta = free_unit**2
    if embeddings[1](eta) < embeddings[0](eta):
        eta = 1 / eta
    if embeddings[1](eta) == embeddings[0](eta):
        raise ArithmeticInconsistency("the oriented square unit has equal embeddings")

    steps = 0
    while embeddings[1](element) < embeddings[0](element):
        element *= eta
        steps += 1
        if steps > 10000:
            raise ArithmeticInconsistency("positive normalization did not converge")
    while embeddings[1](element / eta) >= embeddings[0](element / eta):
        element /= eta
        steps += 1
        if steps > 10000:
            raise ArithmeticInconsistency("positive normalization did not converge")
    if not element.is_integral():
        raise ArithmeticInconsistency("unit-square normalization lost integrality")
    return element


def _positive_unit_classes(unit_group, embeddings):
    """Return canonical exact representatives of U^+/U^2."""
    invariants = tuple(ZZ(value) for value in unit_group.invariants())
    if invariants != (2, 0):
        raise ArithmeticInconsistency("unexpected real-quadratic unit-group structure")
    torsion_unit, free_unit = [generator.value() for generator in unit_group.gens()]
    representatives = []
    parity_coordinates = []
    for torsion_exponent in (0, 1):
        for free_exponent in (0, 1):
            unit = torsion_unit**torsion_exponent * free_unit**free_exponent
            if _signatures(unit, embeddings) != [1, 1]:
                continue
            normalized = _canonical_positive_modulo_square_units(
                unit, free_unit, embeddings
            )
            representatives.append(normalized)
            parity_coordinates.append((torsion_exponent, free_exponent))
    coefficient_map = {
        tuple(_integer_coefficients(unit)): unit for unit in representatives
    }
    if len(coefficient_map) != len(representatives):
        raise ArithmeticInconsistency("positive unit classes collided after normalization")
    if len(set(parity_coordinates)) != len(parity_coordinates):
        raise ArithmeticInconsistency("positive unit parity classes were not distinct")
    return [coefficient_map[key] for key in sorted(coefficient_map)]


def _unit_parity_elements(unit_group):
    invariants = tuple(ZZ(value) for value in unit_group.invariants())
    if invariants != (2, 0):
        raise ArithmeticInconsistency("unexpected real-quadratic unit-group structure")
    torsion_unit, free_unit = [generator.value() for generator in unit_group.gens()]
    return (
        [
            torsion_unit**torsion_exponent * free_unit**free_exponent
            for torsion_exponent in (0, 1)
            for free_exponent in (0, 1)
        ],
        free_unit,
    )


def _ideal_sort_key(ideal):
    hnf = _row_hnf(ideal)
    return (int(ideal.norm()), tuple(entry for row in hnf for entry in row))


def _bounded_class_representatives(K, class_group, D):
    """Choose the least (norm,row-HNF) ideal in every ordinary class."""
    bound = ZZ(D).isqrt() // 2
    candidates = []
    ideals_by_norm = K.ideals_of_bdd_norm(bound)
    for norm in sorted(ideals_by_norm):
        candidates.extend(ideals_by_norm[norm])
    candidates.sort(key=_ideal_sort_key)

    representatives_by_coordinates = {}
    for ideal in candidates:
        coordinates = tuple(int(value) for value in class_group(ideal).list())
        representatives_by_coordinates.setdefault(coordinates, ideal)
    if len(representatives_by_coordinates) != int(class_group.order()):
        raise ArithmeticInconsistency(
            "the Minkowski enumeration did not cover every ordinary ideal class"
        )
    representatives = sorted(representatives_by_coordinates.values(), key=_ideal_sort_key)
    unit_ideal = K.ideal(1)
    for ideal in representatives:
        if ideal.norm() <= 0 or not ideal.is_integral() or ideal * (~ideal) != unit_ideal:
            raise ArithmeticInconsistency(
                "a selected ideal is zero, nonintegral, or noninvertible"
            )
    return int(bound), representatives


def _positive_generators_modulo_square_units(
    squared_ideal, parity_units, free_unit, embeddings
):
    """Directly test whether an ideal square is narrow-principal."""
    try:
        is_principal = squared_ideal.is_principal(proof=True)
    except Exception as error:
        raise PositiveGeneratorFailure("principal-ideal proof failed") from error
    if not is_principal:
        return []
    try:
        generators = squared_ideal.gens_reduced(proof=True)
    except Exception as error:
        raise PositiveGeneratorFailure("principal generator proof failed") from error
    if len(generators) != 1:
        raise ArithmeticInconsistency("a principal ideal lacked one reduced generator")
    generator = generators[0]
    normalized = {}
    for unit in parity_units:
        candidate = unit * generator
        if _signatures(candidate, embeddings) != [1, 1]:
            continue
        candidate = _canonical_positive_modulo_square_units(
            candidate, free_unit, embeddings
        )
        normalized[tuple(_integer_coefficients(candidate))] = candidate
    return [normalized[key] for key in sorted(normalized)]


def _distinct_modulo_unit_squares(unit_group, elements):
    """Check exact pairwise distinction in unit-group exponent coordinates."""
    for left_index, left in enumerate(elements):
        for right in elements[left_index + 1 :]:
            quotient_coordinates = [
                ZZ(value) for value in unit_group(left / right).list()
            ]
            if all(value % 2 == 0 for value in quotient_coordinates):
                return False
    return True


def _make_hm_record(
    K, ideal, ideal_record, unit_id, alpha, record_index, embeddings
):
    alpha_signatures = _signatures(alpha, embeddings)
    norm = ZZ(ideal.norm())
    record_checks = {
        "alpha_integral": bool(alpha.is_integral()),
        "alpha_norm_equals_ideal_norm_squared": bool(alpha.norm() == norm**2),
        "alpha_totally_positive": alpha_signatures == [1, 1],
        "ideal_square_equals_alpha_ideal": bool(ideal**2 == K.ideal(alpha)),
    }
    if not all(record_checks.values()):
        raise ArithmeticInconsistency("an HM record failed an exact check")
    return {
        "alpha": _integer_coefficients(alpha),
        "alpha_signatures": alpha_signatures,
        "checks": record_checks,
        "expected_graph_weight": int(norm),
        "expected_quotient_degree": int(norm**2),
        "geometric_label_status": "THEOREM_CONDITIONAL",
        "id": "hm-{}".format(record_index),
        "ideal": {
            "id": ideal_record["id"],
            "norm": ideal_record["norm"],
            "ordinary_class_coordinates": list(
                ideal_record["ordinary_class_coordinates"]
            ),
            "row_hnf": [list(row) for row in ideal_record["row_hnf"]],
        },
        "ideal_id": ideal_record["id"],
        "unit_class_id": unit_id,
    }


def _status_vocabulary():
    return {
        "success": [SUCCESS],
        "validated_rejection": [
            "INVALID_ARGUMENTS",
            "INVALID_FIELD_DISCRIMINANT",
            "INVALID_ORDER_CONDUCTOR",
            "UNSUPPORTED_NONMAXIMAL_ORDER",
        ],
        "unknown_or_internal": [
            "PROOF_FAILURE_CLASS_GROUP",
            "UNKNOWN_POSITIVE_GENERATOR",
            "INTERNAL_ARITHMETIC_INCONSISTENCY",
        ],
        "future_surface_or_analytic": [
            "UNSUPPORTED_SURFACE_SCOPE",
            "UNKNOWN_GEOMETRIC_ENDOMORPHISMS",
            "UNKNOWN_ANALYTIC_QUOTIENT",
            "UNKNOWN_PPAS_COLLISION",
        ],
    }


def _base_document(status):
    return {
        "schema": SCHEMA,
        "status": status,
        "status_vocabulary": _status_vocabulary(),
        "surface_scope": SURFACE_SCOPE,
        "theorem_status": THEOREM_STATUS,
    }


def _theorem_boundary():
    return {
        "arithmetic_result": "EXACT_PROOF_ENABLED",
        "geometric_interpretation": "CONDITIONAL_ON_HURWITZ_MAASS_THEOREM",
        "quotient_ppas": "NOT_COMPUTED",
        "target_deduplication": "NOT_COMPUTED",
        "future_analytic_success_requires": [
            "CERTIFIED_TARGET_PPAS_AND_POLARIZATION",
            "TRANSPORTED_O_F_EMBEDDING",
            "EXACT_KERNEL_DEGREE_AND_PULLBACK_IDENTITY",
            "RATIONAL_MODEL_WITH_TWIST_AND_DESCENT_DATA",
            "POLARIZED_K_ISOMORPHISM_CERTIFICATE",
        ],
    }


def _order_rejection_theorem_boundary():
    return {
        "arithmetic_result": "EXACT_ORDER_REJECTION_ONLY",
        "class_unit_enumeration": "NOT_RUN",
        "geometric_interpretation": "NOT_EVALUATED",
        "quotient_ppas": "NOT_COMPUTED",
        "target_deduplication": "NOT_COMPUTED",
    }


def _fixed_surface_statement():
    return {
        "status": "CONDITIONAL_ON_POLARIZATION_CORRESPONDENCE",
        "input_scope": [
            "A_GEOMETRICALLY_SIMPLE",
            "CHOSEN_IOTA_O_F_EQUALS_END_K_A",
            "REFERENCE_PRINCIPAL_POLARIZATION_LAMBDA",
            "O_F_MAXIMAL_REAL_QUADRATIC",
        ],
        "polarizations": "lambda_epsilon=lambda o iota(epsilon)",
        "classification": "U_POSITIVE_MODULO_U_SQUARED",
        "changed_rosati": "x |-> eta^-1 x^dagger eta",
        "rosati_on_commutative_End_k0_equals_F": "IDENTITY",
        "larger_geometric_endomorphisms": "OUT_OF_SCOPE",
    }


def _canonical_field(D):
    x = polygen(ZZ, "x")
    polynomial = x**2 - x + (1 - D) // 4 if D % 2 else x**2 - D // 4
    K = NumberField(polynomial, "w")
    return polynomial, K


def enumerate_arithmetic(D, conductor):
    """Return the frozen arithmetic-only result for a validated maximal order."""
    polynomial, K = _canonical_field(D)
    w = K.gen()
    supplied_order = K.order([K(1), conductor * w])
    maximal_order = K.maximal_order()
    if supplied_order != maximal_order:
        document = _base_document("UNSUPPORTED_NONMAXIMAL_ORDER")
        document.update(
            {
                "field": {
                    "defining_polynomial": [int(c) for c in polynomial.list()],
                    "discriminant": int(K.discriminant()),
                },
                "order": {
                    "basis": [[1, 0], [0, int(conductor)]],
                    "conductor": int(conductor),
                    "discriminant": int(supplied_order.discriminant()),
                    "index_in_maximal_order": int(
                        supplied_order.index_in(maximal_order)
                    ),
                    "is_maximal": False,
                },
                "theorem_boundary": _order_rejection_theorem_boundary(),
            }
        )
        return document, 2

    try:
        class_group = K.class_group(proof=True)
        narrow_class_group = K.narrow_class_group(proof=True)
        unit_group = K.unit_group(proof=True)
    except Exception as error:
        raise ProofComputationFailure(
            "class, narrow-class, or unit computation did not certify"
        ) from error
    embeddings = _embeddings(K)

    parity_units, free_unit = _unit_parity_elements(unit_group)
    signature_image = [
        list(signatures)
        for signatures in sorted(
            {tuple(_signatures(unit, embeddings)) for unit in parity_units}
        )
    ]

    positive_units = _positive_unit_classes(unit_group, embeddings)
    unit_classes = [
        {
            "element": _integer_coefficients(unit),
            "id": "unit-{}".format(index),
            "signatures": _signatures(unit, embeddings),
        }
        for index, unit in enumerate(positive_units)
    ]
    minkowski_bound, ordinary_representatives = _bounded_class_representatives(
        K, class_group, D
    )
    eligible = []
    for ideal in ordinary_representatives:
        positive_generators = _positive_generators_modulo_square_units(
            ideal**2, parity_units, free_unit, embeddings
        )
        if positive_generators:
            eligible.append((ideal, positive_generators[0]))

    ideal_classes = []
    for index, (ideal, unused_generator) in enumerate(eligible):
        ideal_classes.append(
            {
                "id": "ideal-{}".format(index),
                "norm": int(ideal.norm()),
                "ordinary_class_coordinates": [
                    int(c) for c in class_group(ideal).list()
                ],
                "row_hnf": _row_hnf(ideal),
            }
        )
    records = []
    record_alphas_by_ideal = []
    for ideal_index, ((ideal, base_generator), ideal_record) in enumerate(
        zip(eligible, ideal_classes)
    ):
        seen_alphas = set()
        exact_alphas = []
        for unit_record, unit in zip(unit_classes, positive_units):
            alpha = _canonical_positive_modulo_square_units(
                base_generator * unit, free_unit, embeddings
            )
            alpha_coefficients = tuple(_integer_coefficients(alpha))
            if alpha_coefficients in seen_alphas:
                raise ArithmeticInconsistency(
                    "two positive-unit classes produced equivalent HM records"
                )
            seen_alphas.add(alpha_coefficients)
            exact_alphas.append(alpha)
            records.append(
                _make_hm_record(
                    K,
                    ideal,
                    ideal_record,
                    unit_record["id"],
                    alpha,
                    len(records),
                    embeddings,
                )
            )
        record_alphas_by_ideal.append(exact_alphas)

    unit_ideal = K.ideal(1)
    canonical_basis = [
        _integer_coefficients(basis_element)
        for basis_element in maximal_order.basis()
    ]
    selected_ideals_ok = all(
        ideal.norm() > 0
        and ideal.is_integral()
        and ideal * (~ideal) == unit_ideal
        for ideal in ordinary_representatives
    )
    eligible_direct_checks = all(
        ideal**2 == K.ideal(generator)
        and _signatures(generator, embeddings) == [1, 1]
        for ideal, generator in eligible
    )
    eligible_class_coordinates = [
        tuple(int(value) for value in class_group(ideal).list())
        for ideal, unused_generator in eligible
    ]
    hm_pairs_distinct = (
        len(set(eligible_class_coordinates)) == len(eligible_class_coordinates)
        and all(
            _distinct_modulo_unit_squares(unit_group, exact_alphas)
            for exact_alphas in record_alphas_by_ideal
        )
    )
    identity_included = (
        any(
            ideal.norm() == 1 and _row_hnf(ideal) == [[1, 0], [0, 1]]
            for ideal, unused_generator in eligible
        )
        and any(_integer_coefficients(unit) == [1, 0] for unit in positive_units)
        and any(record["alpha"] == [1, 0] for record in records)
    )
    run_checks = {
        "all_record_checks": all(
            all(record["checks"].values()) for record in records
        ),
        "canonical_field_basis": canonical_basis == [[1, 0], [0, 1]],
        "canonical_polynomial_discriminant": int(polynomial.discriminant()) == D,
        "eligible_ideals_pass_direct_kernel_of_Sq": eligible_direct_checks,
        "field_discriminant": int(K.discriminant()) == D,
        "hm_pairs_distinct_under_equivalence": hm_pairs_distinct,
        "identity_included": identity_included,
        "maximal_order_discriminant": int(maximal_order.discriminant()) == D,
        "minkowski_class_coverage": len(ordinary_representatives)
        == int(class_group.order()),
        "order_equals_maximal_order": supplied_order == maximal_order,
        "proof_mode": True,
        "record_count_is_cartesian_product": len(records)
        == len(eligible) * len(positive_units),
        "selected_ideals_integral_invertible": selected_ideals_ok,
        "unit_representatives_distinct_modulo_squares": _distinct_modulo_unit_squares(
            unit_group, positive_units
        ),
        "unit_representatives_totally_positive": all(
            _signatures(unit, embeddings) == [1, 1] for unit in positive_units
        ),
    }

    document = _base_document(SUCCESS)
    document.update(
        {
            "arithmetic_pair_equivalence": "(a,alpha)~(beta*a,beta^2*alpha)",
            "checks": run_checks,
            "field": {
                "defining_polynomial": [int(c) for c in polynomial.list()],
                "discriminant": int(K.discriminant()),
            },
            "fixed_surface_polarizations": _fixed_surface_statement(),
            "hm_ideal_classes": ideal_classes,
            "hm_records": records,
            "hm_structure": {
                "eligible_ordinary_class_count": len(eligible),
                "ideal_class_condition": "ker(Sq:Cl(F)->Cl+(F))",
                "positive_unit_class_count": len(positive_units),
                "record_count": len(records),
            },
            "identity_data": {"ideal_id": "ideal-0", "unit_class_id": "unit-0"},
            "minkowski_enumeration": {
                "bound": minkowski_bound,
                "ordinary_class_number": int(class_group.order()),
                "ordinary_classes_covered": len(ordinary_representatives),
            },
            "narrow_class_group": {
                "invariants": [int(n) for n in narrow_class_group.invariants()],
                "order": int(narrow_class_group.order()),
            },
            "order": {
                "basis": [[1, 0], [0, 1]],
                "conductor": 1,
                "discriminant": int(maximal_order.discriminant()),
                "index_in_maximal_order": 1,
                "is_maximal": True,
            },
            "ordinary_class_group": {
                "invariants": [int(n) for n in class_group.invariants()],
                "order": int(class_group.order()),
            },
            "polarization_unit_classes": unit_classes,
            "proof_policy": {
                "class_group": True,
                "narrow_class_group": True,
                "principal_ideal": True,
                "unit_group": True,
            },
            "theorem_boundary": _theorem_boundary(),
            "unit_signature_image": signature_image,
        }
    )
    if not all(document["checks"].values()):
        raise ArithmeticInconsistency("field/order setup failed an exact check")
    return document, 0


def _parse(argv):
    parser = JSONArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    enumerate_parser = subparsers.add_parser("enumerate")
    enumerate_parser.add_argument("--field-discriminant", required=True, type=int)
    enumerate_parser.add_argument("--order-conductor", default=1, type=int)
    return parser.parse_args(argv)


def main(argv=None):
    proof.all(True)
    try:
        args = _parse(sys.argv[1:] if argv is None else argv)
    except CLIValidationError as error:
        print("invalid arguments: {}".format(error), file=sys.stderr)
        _emit(_base_document("INVALID_ARGUMENTS"))
        return 2
    D = ZZ(args.field_discriminant)
    conductor = ZZ(args.order_conductor)
    if D <= 1 or not D.is_fundamental_discriminant():
        document = _base_document("INVALID_FIELD_DISCRIMINANT")
        document["input"] = {
            "field_discriminant": int(D),
            "order_conductor": int(conductor),
        }
        _emit(document)
        return 2
    if conductor <= 0:
        document = _base_document("INVALID_ORDER_CONDUCTOR")
        document["input"] = {
            "field_discriminant": int(D),
            "order_conductor": int(conductor),
        }
        _emit(document)
        return 2
    try:
        document, exit_code = enumerate_arithmetic(D, conductor)
    except ProofComputationFailure as error:
        print("proof failure: {}".format(error), file=sys.stderr)
        document = _base_document("PROOF_FAILURE_CLASS_GROUP")
        exit_code = 3
    except PositiveGeneratorFailure as error:
        print("positive-generator failure: {}".format(error), file=sys.stderr)
        document = _base_document("UNKNOWN_POSITIVE_GENERATOR")
        exit_code = 3
    except ArithmeticInconsistency as error:
        print("internal arithmetic inconsistency: {}".format(error), file=sys.stderr)
        document = _base_document("INTERNAL_ARITHMETIC_INCONSISTENCY")
        exit_code = 3
    except Exception as error:
        print("internal failure: {}: {}".format(type(error).__name__, error), file=sys.stderr)
        document = _base_document("INTERNAL_ARITHMETIC_INCONSISTENCY")
        exit_code = 3
    _emit(document)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

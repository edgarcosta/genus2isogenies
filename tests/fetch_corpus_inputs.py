"""Fetch the LMFDB corpus inputs for the isogeny-primes engine.

Writes tests/corpus_inputs.txt: the elliptic factors, over their quadratic
splitting fields, of the LMFDB genus 2 classes whose Jacobian is geometrically
Q x Q but splits only over a quadratic field. One factor per line, in the
schema documented in that file's header and shared with the reader: the rows
are handed to generate_corpus.write_inputs, so the fetcher and
generate_corpus.read_inputs stay in lockstep on the format.

Data comes straight from the LMFDB over lmfdb-lite (roed314/lmfdb-lite, module
`lmf`), reproducing this query as two searches joined on the label in Python
(g2c_curves supplies the geom_end_alg filter, g2c_endomorphisms the fields):

    SELECT c.label, e.spl_fod_label, e.spl_fod_coeffs, e.spl_facs_coeffs,
           e.spl_facs_labels
    FROM g2c_curves c JOIN g2c_endomorphisms e ON c.label = e.label
    WHERE c.geom_end_alg = 'Q x Q' AND e.spl_fod_label != '1.1.1.1'
    ORDER BY c.label

Run from the repo root (needs network access to the LMFDB):

    sage -pip install git+https://github.com/roed314/lmfdb-lite   # once
    sage -python tests/fetch_corpus_inputs.py [out_path]

out_path defaults to tests/corpus_inputs.txt.
"""

import sys

from lmf import db

from generate_corpus import INPUTS_PATH, _blist, _vec, write_inputs

# Provenance recorded in the committed corpus_inputs.txt header. FETCHED is the
# canonical snapshot date of the committed dataset; bump it only when you
# deliberately re-fetch and re-commit changed LMFDB content.
FETCHED = "2026-07-20"
SOURCE = "LMFDB g2c_curves + g2c_endomorphisms"
SQL = (
    "SELECT c.label, e.spl_fod_label, e.spl_fod_coeffs, e.spl_facs_coeffs, "
    "e.spl_facs_labels FROM g2c_curves c JOIN g2c_endomorphisms e "
    "ON c.label = e.label WHERE c.geom_end_alg = 'Q x Q' "
    "AND e.spl_fod_label != '1.1.1.1' ORDER BY c.label"
)

PROJECTION = ["label", "spl_fod_label", "spl_fod_coeffs", "spl_facs_coeffs",
              "spl_facs_labels"]


def fetch_rows():
    """The factor rows for write_inputs, in SQL (ORDER BY label) order: one
    tuple (label, factor_index, spl_fod_label, factor_label, field, curve) per
    elliptic factor. field is spl_fod_coeffs and curve is the factor's
    [A, B] coordinate-vector pair, both serialized as bracketed lists.
    factor_label is '' where the LMFDB has none: spl_facs_labels is then a NULL
    column, which lmfdb-lite drops from the row dict."""
    qxq = list(db.g2c_curves.search({"geom_end_alg": "Q x Q"}, "label"))
    endo = db.g2c_endomorphisms.search(
        {"label": {"$in": qxq}, "spl_fod_label": {"$ne": "1.1.1.1"}}, PROJECTION)
    rows = []
    for r in sorted(endo, key=lambda row: row["label"]):
        field = r["spl_fod_coeffs"]
        assert len(field) == 3, "%s: splitting field not quadratic: %r" % (
            r["label"], field)
        labels = r.get("spl_facs_labels") or []
        for i, fac in enumerate(r["spl_facs_coeffs"]):
            assert len(fac) == 2, "%s: factor %d is not an [A, B] pair" % (
                r["label"], i)
            A, B = fac
            faclabel = (labels[i] if i < len(labels) else None) or ""
            rows.append((r["label"], str(i), r["spl_fod_label"], faclabel,
                         _vec(field), _blist([_vec(A), _vec(B)])))
    return rows


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else INPUTS_PATH
    write_inputs(fetch_rows(), FETCHED, SOURCE, SQL, path)


if __name__ == "__main__":
    main()

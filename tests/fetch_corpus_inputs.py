"""Fetch and transform the LMFDB corpus inputs for the isogeny-primes engine.

The corpus's `lmfdb-qcurve` stratum is the set of elliptic factors, over their
quadratic splitting fields, of the 125 LMFDB genus 2 classes whose Jacobian is
geometrically Q x Q but splits only over a quadratic field (Weil restrictions).

The data comes from this LMFDB SQL query (verified 2026-07-20, 125 rows):

    SELECT c.label, e.spl_fod_label, e.spl_fod_coeffs, e.spl_facs_coeffs,
           e.spl_facs_labels
    FROM g2c_curves c JOIN g2c_endomorphisms e ON c.label = e.label
    WHERE c.geom_end_alg = 'Q x Q' AND e.spl_fod_label != '1.1.1.1'
    ORDER BY c.label

The query is executed out-of-band (orchestrator MCP access or the LMFDB SQL
mirror) and its JSONL result is fed to transform() below, which writes
tests/corpus_curves.json in the schema consumed by generate_corpus.py:
each factor becomes one entry {id, stratum, field, curve} where `field` is the
ascending defining polynomial of K and `curve` is the short Weierstrass model
y^2 = x^3 + A x + B with A, B as coordinate vectors in the power basis of K.
"""

import datetime
import json
import sys

SQL = (
    "SELECT c.label, e.spl_fod_label, e.spl_fod_coeffs, e.spl_facs_coeffs, "
    "e.spl_facs_labels FROM g2c_curves c JOIN g2c_endomorphisms e "
    "ON c.label = e.label WHERE c.geom_end_alg = 'Q x Q' "
    "AND e.spl_fod_label != '1.1.1.1' ORDER BY c.label"
)


def transform(jsonl_path, out_path):
    entries = []
    with open(jsonl_path) as f:
        rows = [json.loads(line) for line in f if line.strip()]
    if len(rows) != 125:
        raise SystemExit(f"expected 125 rows, got {len(rows)}")
    for row in rows:
        field = row["spl_fod_coeffs"]
        if len(field) != 3:
            raise SystemExit(f"{row['label']}: splitting field not quadratic: {field}")
        labels = row.get("spl_facs_labels") or []
        for i, fac in enumerate(row["spl_facs_coeffs"]):
            if len(fac) != 2:
                raise SystemExit(f"{row['label']}: factor {i} is not an [A, B] pair")
            A, B = fac
            entries.append(
                {
                    "id": f"lmfdb-{row['label']}-{i}",
                    "stratum": "lmfdb-qcurve",
                    "lmfdb_factor_label": labels[i] if i < len(labels) else None,
                    "field": field,
                    "curve": {"model": "AB", "A": A, "B": B},
                }
            )
    data = {
        "_provenance": {
            "date": datetime.date.today().isoformat(),
            "source": "LMFDB g2c_curves + g2c_endomorphisms",
            "sql": SQL,
            "magma_version": "2.29-8",
            "sage_version": "10.8",
        },
        "entries": entries,
    }
    with open(out_path, "w") as f:
        json.dump(data, f, indent=1)
    print(f"wrote {out_path}: {len(entries)} entries")


def main():
    if len(sys.argv) == 3:
        transform(sys.argv[1], sys.argv[2])
    else:
        raise SystemExit(
            "usage: fetch_corpus_inputs.py <rows.jsonl> <corpus_curves.json>\n"
            "obtain rows.jsonl by running the SQL in the module docstring "
            "against the LMFDB (orchestrator MCP or SQL mirror)"
        )


if __name__ == "__main__":
    main()

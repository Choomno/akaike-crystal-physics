"""Audit every archived XLS run, then write the canonical figure data.

Run from any working directory. A failed numerical check exits nonzero.
"""

import csv
import hashlib
import json
import platform

import matplotlib
import numpy as np
import scipy
import xlrd
from models import ROOT, aicc, datasets, debye, einstein, nonlinear, polynomial, weights

OUT = ROOT / "docs" / "generated"


def sheet(path):
    return xlrd.open_workbook(str(path)).sheet_by_index(0)


def close(actual, expected, tolerance, description):
    if not np.isclose(actual, expected, rtol=0.0, atol=tolerance):
        raise AssertionError(
            f"{description}: {actual} != {expected} (atol={tolerance})"
        )


def ranked(records):
    scores = [record["aicc"] for record in records]
    ranks = np.argsort(np.argsort(scores)) + 1
    for record, weight, rank in zip(records, weights(scores), ranks):
        record.update(weight=float(weight), rank=int(rank))
    return records


def main():
    result = {
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "xlrd": xlrd.__version__,
        },
        "datasets": {},
        "checks": [],
    }
    for name, x, y, maximum in datasets():
        path = ROOT / "results" / "2020-04_polynomial" / name
        table = sheet(path / "AIC.xls")
        archived, corrected = [], []
        for degree in range(1, maximum + 1):
            parameters = sheet(path / f"degree_{degree}.xls")
            beta = [parameters.cell_value(i, 0) for i in range(1, degree + 2)]
            score = float(table.cell_value(degree, 0))
            residual = y - np.polynomial.polynomial.polyval(x, beta)
            close(
                aicc(residual @ residual, len(y), degree + 2),
                score,
                1e-7,
                f"{name} degree {degree}: stored parameters reproduce stored AICc",
            )
            fit = polynomial(x, y, degree)
            # Normal equations in the historical program become inaccurate at high degree.
            tolerance = (
                1e-7
                if degree <= 6
                else 1e-5
                if degree == 7
                else 1e-2
                if degree == 8
                else None
            )
            if tolerance is not None:
                close(
                    fit["aicc"],
                    score,
                    tolerance,
                    f"{name} degree {degree}: independent fit",
                )
            else:
                # The scaled solution should minimize RSS at least as well as the historical fit.
                if fit["aicc"] > score + 1e-7:
                    raise AssertionError(
                        f"{name} degree {degree}: scaled fit worsened RSS"
                    )
            archived.append(
                {
                    "model": fit["model"],
                    "degree": degree,
                    "parameters": beta,
                    "aicc": score,
                    "adjusted_r2": float(parameters.cell_value(1, 2)),
                }
            )
            corrected.append(fit)
            result["checks"].append(
                {
                    "dataset": name,
                    "degree": degree,
                    "archived_aicc": score,
                    "refit_aicc": fit["aicc"],
                    "difference": fit["aicc"] - score,
                    "reproduction_tolerance": tolerance,
                    "status": "verified"
                    if tolerance
                    else "historical numerical discrepancy",
                }
            )
        ranked(archived)
        for degree, record in enumerate(archived, 1):
            close(
                record["weight"],
                table.cell_value(degree, 1),
                1e-10,
                f"{name}: archived weights",
            )
        result["datasets"][name] = {
            "x": x.tolist(),
            "y": y.tolist(),
            "archived": archived,
            "corrected": ranked(corrected),
        }

    path = ROOT / "results" / "2020-11_modular_swaminathan"
    table, parameters = sheet(path / "AIC.xls"), sheet(path / "nonlinear.xls")
    data = result["datasets"]["swaminathan_a"]
    x, y = np.array(data["x"]), np.array(data["y"])
    archived, corrected = [], []
    for degree in range(1, 4):
        parameter_sheet = sheet(path / f"degree_{degree}.xls")
        beta = [parameter_sheet.cell_value(i, 0) for i in range(1, degree + 2)]
        residual = y - np.polynomial.polynomial.polyval(x, beta)
        score = table.cell_value(degree, 1)
        close(
            aicc(residual @ residual, len(y), degree + 2),
            score,
            1e-7,
            "modular polynomial parameters",
        )
        fit = polynomial(x, y, degree)
        close(fit["aicc"], score, 1e-7, "modular polynomial fit")
        archived.append({"model": fit["model"], "parameters": beta, "aicc": score})
        corrected.append(fit)
    for column, (name, function) in enumerate(
        (("Einstein", einstein), ("Debye", debye))
    ):
        beta = [parameters.cell_value(i, column) for i in range(1, 4)]
        score = table.cell_value(4 + column, 1)
        residual = y - function(beta, x)
        close(
            aicc(residual @ residual, len(y), 4),
            score,
            1e-7,
            f"archived {name} parameters",
        )
        archived.append({"model": name, "parameters": beta, "aicc": score})
        fit = nonlinear(name, x, y)
        if name == "Einstein":
            close(fit["aicc"], score, 1e-4, "Einstein independent refit")
        elif fit["aicc"] >= score - 1:
            raise AssertionError(
                "Debye refit did not resolve the documented early stop"
            )
        corrected.append(fit)
    ranked(archived)
    for i, row in enumerate(archived, 1):
        close(row["weight"], table.cell_value(i, 2), 1e-10, "modular weights")
    result["urea_models"] = {"archived": archived, "corrected": ranked(corrected)}
    # Fingerprint all immutable inputs; figures reject a stale canonical data file.
    inputs = sorted((ROOT / "data").glob("*.txt")) + sorted(
        (ROOT / "results").rglob("*.xls")
    )
    sources = [ROOT / "docs/scripts" / name for name in ("verify.py", "models.py")]
    result["inputs"] = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in inputs + sources
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "audit.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    with (OUT / "polynomial_audit.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(result["checks"][0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(result["checks"])
    with (OUT / "urea_rankings.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "model",
                "archived_aicc",
                "archived_weight",
                "archived_rank",
                "refit_aicc",
                "refit_weight",
                "refit_rank",
            ]
        )
        for old, new in zip(
            result["urea_models"]["archived"], result["urea_models"]["corrected"]
        ):
            writer.writerow(
                [
                    old["model"],
                    old["aicc"],
                    old["weight"],
                    old["rank"],
                    new["aicc"],
                    new["weight"],
                    new["rank"],
                ]
            )
    print(
        f"PASS: {len(result['checks'])} polynomial fits across 8 datasets; all {len(inputs) - 3} XLS files read; modular fits checked."
    )
    print(
        "Canonical results: docs/generated/audit.json; comparison tables: docs/generated/*.csv"
    )


if __name__ == "__main__":
    main()

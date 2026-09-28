"""Sequential exact-certificate runner.

Run: python -B -O verify_all.py --group all
The runner uses only the standard library; lattice and histogram checks
additionally require NumPy. Children run sequentially with -B -O from
this directory. Reports go under results/ by default. The analytic and
operator arguments remain separate parts of the accompanying manuscript.
"""

import argparse
import ast
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import monotonic


ROOT = Path(__file__).resolve().parent
GROUPS = {
    "arithmetic": [
        ("circle_constants", "verify_circle_tail_constants.py", []),
        ("rotation_gram", "verify_rotation_gram.py", []),
        ("angular_ninth_eleventh", "verify_low_mode_rotation_repairs.py", []),
        ("angular_fifth", "verify_fifth_mode_upper_repair.py", []),
        ("small_histogram_audit", "audit_square_corridor_histogram_small.py", []),
        ("lattice_048_070", "replay_all_angles_band.py", []),
        ("lattice_070_080", "replay_midband_independent.py", []),
        ("lattice_080_083", "replay_upper_midband_independent.py", []),
        ("lattice_083_090", "replay_compact_band_independent.py", []),
        ("lattice_090_098", "replay_final_compact_independent.py", []),
    ],
    "geometry": [
        ("radial_scalar", "verify_radial_scalar.py", []),
        ("cap_098", "verify_cap_closure_constants.py", []),
        ("cap_092_x8", "verify_cap_q92_x8_constants.py", []),
    ],
    "inherited": [
        ("small_hole_deficits", "verify_deficits.py", []),
        ("radius_201_corridor", "verify_corridor.py", []),
        ("rectangular_holes", "verify_rectangular.py", []),
    ],
    "histograms": [
        ("histogram_50000", "certify_square_corridor_histogram.py",
         ["--radius", "50000", "--bins-per-unit", "100", "--paired",
          "--eta", "1/14", "--eta", "1/19"]),
        ("histogram_fine_5000", "certify_square_corridor_histogram.py",
         ["--radius", "5000", "--bins-per-unit", "1000", "--paired",
          "--eta", "1/19", "--eta", "1/24"]),
        ("histogram_90000", "certify_square_corridor_histogram.py",
         ["--radius", "90000", "--bins-per-unit", "100", "--paired",
          "--eta", "1/24"]),
    ],
}
HISTOGRAM_REPORTS = {
    "histogram_50000": "square_corridor_paired_50000.json",
    "histogram_fine_5000": "square_corridor_paired_fine_5000.json",
    "histogram_90000": "square_corridor_paired_90000.json",
}
DERIVED_REPORTS = {
    "angular_ninth_eleventh": "low_mode_rotation_repairs.json",
    "angular_fifth": "fifth_mode_upper_repair.json",
}
WITNESS_FILES = [
    "all_angles_certificate.json", "certificate_midband.json",
    "certificate_upper_midband.json", "certificate_compact.json",
    "certificate_final.json",
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True) + "\n",
                    encoding="utf-8")


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_audit():
    """Reject assert statements and undeclared nonstandard imports.

    Explicit raise AssertionError remains active under -O and is allowed.
    This is a source/import inventory, not an algorithmic proof.
    """
    require(hasattr(sys, "stdlib_module_names"), "Python 3.10 or newer is required")
    sources = sorted(ROOT.glob("*.py"))
    allowed = set(sys.stdlib_module_names) | {p.stem for p in sources} | {"numpy"}
    result = []
    for path in sources:
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=path.name)
        assertions = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)]
        require(not assertions, path.name + ": assert disabled by -O: " + str(assertions))
        imports = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(item.name.split(".")[0] for item in node.names)
            elif isinstance(node, ast.ImportFrom):
                require(node.level == 0, path.name + ": unexpected relative import")
                if node.module:
                    imports.add(node.module.split(".")[0])
        require(imports <= allowed,
                path.name + ": undeclared imports " + str(sorted(imports - allowed)))
        result.append({"file": path.name, "sha256": file_hash(path),
                       "assert_statements": assertions, "imports": sorted(imports)})
    return result


def parsed_output(text):
    """Collect complete JSON objects printed among progress/plain-text lines."""
    decoder = json.JSONDecoder()
    objects, offset = [], 0
    while offset < len(text):
        start = text.find("{", offset)
        if start < 0:
            break
        try:
            value, length = decoder.raw_decode(text[start:])
        except json.JSONDecodeError:
            offset = start + 1
        else:
            objects.append(value)
            offset = start + length
    return objects


def mathematical_report(report):
    return {k: v for k, v in report.items() if k != "seconds"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", action="append", choices=["all", *GROUPS],
                        help="Repeat for selected groups; default all includes large histograms")
    parser.add_argument("--output-root", type=Path, default=Path("results"),
                        help="Relative paths are resolved beside this runner")
    parser.add_argument("--list", action="store_true", help="List checks without executing")
    args = parser.parse_args()
    chosen = args.group or ["all"]
    require("all" not in chosen or chosen == ["all"], "Use --group all alone")
    chosen = list(GROUPS) if chosen == ["all"] else list(dict.fromkeys(chosen))
    jobs = [job for group in chosen for job in GROUPS[group]]
    if args.list:
        print(json.dumps({"groups": chosen, "checks": [j[0] for j in jobs]}, indent=2))
        return 0

    output_root = args.output_root
    if not output_root.is_absolute():
        output_root = ROOT / output_root
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "verification_results.json"
    started = monotonic()
    try:
        numpy_version = metadata.version("numpy")
    except metadata.PackageNotFoundError:
        numpy_version = None
    summary = {
        "status": "RUNNING", "started_utc": stamp(), "groups": chosen,
        "python": sys.version, "python_executable": sys.executable,
        "platform": platform.platform(), "numpy_version": numpy_version,
        "child_flags": ["-B", "-O"],
        "execution_directory": "directory containing verify_all.py",
        "analytic_scope": "Computational certificates only; analytic and operator proofs are in the manuscript.",
        "checks": [],
    }
    write_json(summary_path, summary)
    try:
        summary["source_audit"] = source_audit()
        required = WITNESS_FILES + list(HISTOGRAM_REPORTS.values())
        for filename in required:
            require((ROOT / filename).is_file(), "Missing release input: " + filename)
        summary["input_sha256"] = {name: file_hash(ROOT / name) for name in required}
        write_json(summary_path, summary)
        for check_id, script, options in jobs:
            script_path = ROOT / script
            require(script_path.is_file(), "Missing verifier: " + script)
            options, report_path = list(options), None
            if check_id in HISTOGRAM_REPORTS:
                report_path = output_root / HISTOGRAM_REPORTS[check_id]
                require(report_path != ROOT / HISTOGRAM_REPORTS[check_id],
                        "Histogram output must not overwrite the supplied reference report")
                options += ["--output", str(report_path)]
            if check_id in DERIVED_REPORTS:
                require(not (ROOT / DERIVED_REPORTS[check_id]).exists(),
                        "Move the pre-existing derived report before running: "
                        + DERIVED_REPORTS[check_id])
            command = [sys.executable, "-B", "-O", str(script_path), *options]
            print("START " + check_id, flush=True)
            record = {"id": check_id, "script": script, "status": "RUNNING",
                      "command": command, "started_utc": stamp()}
            summary["checks"].append(record)
            write_json(summary_path, summary)
            job_start = monotonic()
            process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True,
                                       encoding="utf-8", errors="replace", bufsize=1)
            captured = []
            for line in process.stdout:
                print(line, end="", flush=True)
                captured.append(line)
            code = process.wait()
            record.update({"exit_code": code,
                           "elapsed_seconds": round(monotonic()-job_start, 6),
                           "finished_utc": stamp(), "stdout": "".join(captured)})
            record["reported_json"] = parsed_output(record["stdout"])
            require(code == 0, check_id + " returned nonzero exit status " + str(code))
            if check_id in DERIVED_REPORTS:
                generated = ROOT / DERIVED_REPORTS[check_id]
                require(generated.is_file() and generated.parent == ROOT,
                        "Missing generated angular report")
                report_path = output_root / generated.name
                if generated != report_path:
                    os.replace(generated, report_path)
                record["derived_report"] = report_path.name
                record["derived_sha256"] = file_hash(report_path)
            if check_id in HISTOGRAM_REPORTS:
                fresh = json.loads(report_path.read_text(encoding="utf-8-sig"))
                reference = json.loads((ROOT / HISTOGRAM_REPORTS[check_id]).read_text(encoding="utf-8-sig"))
                require(mathematical_report(fresh) == mathematical_report(reference),
                        check_id + ": fresh mathematical report differs from supplied reference")
                record["reference_report_match_excluding_seconds"] = True
                record["derived_report"] = report_path.name
                record["derived_sha256"] = file_hash(report_path)
            record["status"] = "PASS"
            write_json(output_root / (check_id + ".json"), record)
            write_json(summary_path, summary)
            print("PASS " + check_id + " (" + str(record["elapsed_seconds"]) + " s)", flush=True)
        summary["status"] = "PASS"
        summary["completed_checks"] = len(jobs)
    except Exception as error:
        summary["status"] = "FAIL"
        summary["error"] = type(error).__name__ + ": " + str(error)
        if summary["checks"] and summary["checks"][-1]["status"] == "RUNNING":
            summary["checks"][-1]["status"] = "FAIL"
        print(summary["error"], file=sys.stderr, flush=True)
    finally:
        summary["finished_utc"] = stamp()
        summary["elapsed_seconds"] = round(monotonic()-started, 6)
        write_json(summary_path, summary)
    print(json.dumps({"status": summary["status"],
                      "completed_checks": sum(r["status"] == "PASS" for r in summary["checks"]),
                      "requested_checks": len(jobs),
                      "elapsed_seconds": summary["elapsed_seconds"],
                      "results": str(summary_path)}, indent=2), flush=True)
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())

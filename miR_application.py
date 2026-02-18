#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path
import uuid

def runner(inp_name: str, out_name: str) -> None:
    # Cartella dove si trova miR_application.py (main/)
    main_dir = Path(__file__).resolve().parent

    # Cartella da cui lanci il comando (output qui)
    run_dir = Path.cwd()

    utilities_dir = main_dir / "utilities" / "prediction"
    py_script = utilities_dir / "PY_miR_features_extraction.py"
    r_script = utilities_dir / "make_miR_pred.R"

    # INPUT: supporta path assoluto (~/...) oppure relativo a main/input/ ---
    inp_path = Path(inp_name).expanduser()
    if not inp_path.is_absolute():
        inp_path = (run_dir / inp_path).resolve()
    else:
        inp_path = inp_path.resolve()

    # OUTPUT nella cartella da cui lanci il codice ---
    out_path = (run_dir / out_name).resolve()

    # Run id unico + cartella dedicata per intermedi ---
    run_id = uuid.uuid4().hex[:10]
    run_workdir = utilities_dir / "runs" / run_id
    run_workdir.mkdir(parents=True, exist_ok=True)

    # File intermedio unico (usato da R) ---
    intermediate_file = run_workdir / f"{inp_path.name}_featureTable.tsv"

    # Check minimi
    if not utilities_dir.is_dir():
        raise FileNotFoundError(f"utilities/prediction not found: {utilities_dir}")
    if not py_script.is_file():
        raise FileNotFoundError(f"Python script not found: {py_script}")
    if not r_script.is_file():
        raise FileNotFoundError(f"R script not found: {r_script}")
    if not inp_path.is_file():
        raise FileNotFoundError(f"Input file not found: {inp_path}")

    # 1) Step python
    # Adesso passo i 3 argomenti:
    #   - input file
    #   - run_workdir (dove scrivere intermedi)
    #   - intermediate_file (feature table finale)
    subprocess.run(
        ["python3", str(py_script), str(inp_path), str(run_workdir), str(intermediate_file)],
        check=True,
        cwd=str(utilities_dir)
    )

    # sicurezza: se non esiste, meglio fermarsi qui con errore chiaro
    if not intermediate_file.is_file():
        raise FileNotFoundError(f"Feature table not created: {intermediate_file}")

    # 2) Step R
    subprocess.run(
        ["Rscript", str(r_script), str(intermediate_file), str(out_path)],
        check=True,
        cwd=str(utilities_dir)
    )

    #print(f"[] Output prediction: {out_path}")

        # Sotto per cancellare anche la feature table dopo R:
    try:
        intermediate_file.unlink()
    except FileNotFoundError:
        pass

    # Cancella anche la cartella run/ se è vuota
    try:
        run_workdir.rmdir()
    except OSError:
        pass

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 miR_application.py <INPUT> <OUTPUT>")
        print("Example (relative): python3 miR_application.py RNAfold_prova.txt outpred.tsv")
        print("Example (absolute): python3 miR_application.py ~/proveConsistent_path/input/RNAfold_prova.txt outpred.tsv")
        sys.exit(1)

    try:
        runner(sys.argv[1], sys.argv[2])
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    sys.exit(0)

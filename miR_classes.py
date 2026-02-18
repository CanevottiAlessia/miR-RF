import sys
import argparse
import subprocess
from pathlib import Path
import uuid


def resolve_path(run_dir: Path, user_out: str) -> Path:
    p = Path(user_out).expanduser()
    if not p.is_absolute():
        p = (run_dir / p).resolve()
    else:
        p = p.resolve()
    return p


def runner(rnafold: str, fasta: str, base_pred: str, out: str) -> None:
    main_dir = Path(__file__).resolve().parent
    run_dir = Path.cwd()

    classes_dir = main_dir / "utilities" / "classes"
    if not classes_dir.is_dir():
        raise FileNotFoundError(f"utilities/classes not found: {classes_dir}")

    # Script paths (assoluti)
    get_fasta_py = classes_dir / "get_fasta_for_insertion.py"
    snp_py = classes_dir / "SNP_insertion.py"
    feat_py = classes_dir / "PY_miR_features_extraction.py"
    merge_py = classes_dir / "merge_table.py"
    get_lens_py = classes_dir / "get_lens.py"
    fisher_py = classes_dir / "make_fisher_test.py"
    first_py = classes_dir / "make_first_classes.py"
    final_py = classes_dir / "make_final_classes.py"
    pred_r = classes_dir / "make_miR_pred.R"

    for s in [get_fasta_py, snp_py, feat_py, merge_py, get_lens_py, fisher_py, first_py, final_py, pred_r]:
        if not s.is_file():
            raise FileNotFoundError(f"Required file not found: {s}")

    # Risolvo input reali
    rnafold_path = resolve_path(run_dir, rnafold)
    fasta_path = resolve_path(run_dir, fasta)
    base_pred_path = resolve_path(run_dir, base_pred)

    if not rnafold_path.is_file():
        raise FileNotFoundError(f"RNAfold input not found: {rnafold_path}")
    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA input not found: {fasta_path}")
    if not base_pred_path.is_file():
        raise FileNotFoundError(f"Base prediction file not found: {base_pred_path}")

    # Risolvo output
    out_path = resolve_path(run_dir, out)

    # Run workdir dedicata
    run_id = uuid.uuid4().hex[:10]
    runs_root = classes_dir / "runs"
    run_workdir = runs_root / run_id
    run_workdir.mkdir(parents=True, exist_ok=True)

    # 1) get_fasta_for_insertion.py
    subprocess.run(
        ["python3", str(get_fasta_py), str(rnafold_path), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    intermediate_hairpins = run_workdir / f"{rnafold_path.name}_final_hairpins.txt"
    if not intermediate_hairpins.is_file():
        raise FileNotFoundError(f"Expected intermediate not created: {intermediate_hairpins}")

    # 2) SNP_insertion.py <fasta> <hairpins_file>
    subprocess.run(
        ["python3", str(snp_py), str(fasta_path), str(intermediate_hairpins), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    potential_rnafold = run_workdir / f"{intermediate_hairpins.name}_potential_RNAfold.txt"
    if not potential_rnafold.is_file():
        raise FileNotFoundError(f"Expected intermediate not created: {potential_rnafold}")

    # 3) Feature extraction su potential_rnafold
    subprocess.run(
        ["python3", str(feat_py), str(potential_rnafold), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    feature_table = run_workdir / f"{potential_rnafold.name}_features_table.tsv"
    if not feature_table.is_file():
        raise FileNotFoundError(f"Feature table not created: {feature_table}")

    # 4) Predizione R su feature_table (scrive *_pred.tsv in run_workdir)
    subprocess.run(
        ["Rscript", str(pred_r), str(feature_table), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    potential_pred = run_workdir / f"{potential_rnafold.name}_pred.tsv"
    if not potential_pred.is_file():
        raise FileNotFoundError(f"Potential prediction not created: {potential_pred}")

    # 5) merge_table.py <base_pred> <potential_pred>
    subprocess.run(
        ["python3", str(merge_py), str(base_pred_path), str(potential_pred), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    classes_table = run_workdir / f"{base_pred_path.name}_temp_classes_table.tsv"
    if not classes_table.is_file():
        raise FileNotFoundError(f"Classes table not created: {classes_table}")

    # 6) get_lens.py <rnafold>
    subprocess.run(
        ["python3", str(get_lens_py), str(rnafold_path), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    lens_file = run_workdir / f"{rnafold_path.name}_temp_len_hairpins.txt"
    if not lens_file.is_file():
        raise FileNotFoundError(f"Lens file not created: {lens_file}")

    # 7) fisher test
    subprocess.run(
        ["python3", str(fisher_py), str(classes_table), str(lens_file), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    dr_class = run_workdir / f"{classes_table.name}_temp_DR_class.tsv"
    is_class = run_workdir / f"{classes_table.name}_temp_IS_class.tsv"
    if not dr_class.is_file() or not is_class.is_file():
        raise FileNotFoundError("DR/IS class tables not created as expected.")

    # 8) first classes
    subprocess.run(
        ["python3", str(first_py), str(classes_table), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    classes_mirnas = run_workdir / f"{classes_table.name}_temp_classes_mirna.txt"
    if not classes_mirnas.is_file():
        raise FileNotFoundError(f"Classes mirnas file not created: {classes_mirnas}")

    # 9) final classes -> out_path (assoluto, quindi va dove lanci il comando)
    subprocess.run(
        ["python3", str(final_py), str(dr_class), str(is_class), str(classes_mirnas), str(out_path)],
        check=True,
        cwd=str(run_workdir)
    )

    # --- CLEANUP: cancella SOLO i file dentro run_workdir (ricorsivo), lasciando le cartelle ---
    for p in run_workdir.rglob("*"):
        if p.is_file() or p.is_symlink():
            try:
                p.unlink()
            except FileNotFoundError:
                pass

    # Non rimuovo run_workdir: resta vuota (come richiesto)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="miRF classes wrapper (Strada A)")
    parser.add_argument("rnafold", type=str, help="RNAfold file (relative to cwd or absolute)")
    parser.add_argument("fasta", type=str, help="FASTA file (relative to cwd or absolute)")
    parser.add_argument("base_pred", type=str, help="Prediction file from previous step (relative to cwd or absolute)")
    parser.add_argument("out", type=str, help="Output file for final classes (relative to cwd or absolute)")
    args = parser.parse_args()

    out_preview = resolve_path(Path.cwd(), args.out)
    if out_preview.exists():
        resp = input(f"Warning: '{out_preview}' exists. Overwrite? (y/n): ").strip().lower()
        if resp != "y":
            print("Operation cancelled by the user.")
            sys.exit(1)

    try:
        runner(args.rnafold, args.fasta, args.base_pred, args.out)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    sys.exit(0)

import sys
import argparse
import shutil
import subprocess
from pathlib import Path
import uuid
import os


"""def resolve_input(main_dir: Path, user_path: str) -> Path:
    p = Path(user_path).expanduser()
    if not p.is_absolute():
        p = (main_dir / "input" / p).resolve()
    else:
        p = p.resolve()
    return p
"""

def resolve_path(run_dir: Path, user_out: str) -> Path:
    p = Path(user_out).expanduser()
    if not p.is_absolute():
        p = (run_dir / p).resolve()
    else:
        p = p.resolve()
    return p


'''def safe_link_or_copy(src: Path, dst: Path) -> None:
    """
    Prova symlink (più veloce), se non possibile copia.
    """
    if dst.exists():
        return
    try:
        os.symlink(src, dst)
    except Exception:
        shutil.copy2(src, dst)
'''

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
    #model_rds = classes_dir / "trained_model_new.RDS"

    for s in [get_fasta_py, snp_py, feat_py, merge_py, get_lens_py, fisher_py, first_py, final_py, pred_r]: #, model_rds]:
        if not s.is_file():
            raise FileNotFoundError(f"Required file not found: {s}")

    # Risolvo input reali
    rnafold_path = resolve_path(run_dir, rnafold)
    fasta_path = resolve_path(run_dir, fasta)
    base_pred_path = resolve_path(run_dir,base_pred)
    
    '''if not base_pred_path.is_absolute():
        # 1) prova relativo a dove lanci
        cand = (run_dir / base_pred_path).resolve()
        if cand.exists():
            base_pred_path = cand
        else:
            # 2) fallback su main/input/
            base_pred_path = (main_dir / "input" / base_pred_path).resolve()
    else:
        base_pred_path = base_pred_path.resolve()'''


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

    '''# Copio input nella run dir con nomi “puliti”
    rnafold_local = run_workdir / rnafold_path.name
    fasta_local = run_workdir / fasta_path.name
    base_pred_local = run_workdir / base_pred_path.name

    shutil.copy2(rnafold_path, rnafold_local)
    shutil.copy2(fasta_path, fasta_local)
    shutil.copy2(base_pred_path, base_pred_local)

    # R deve trovare trained_model_new.RDS “in cwd”
    #safe_link_or_copy(model_rds, run_workdir / model_rds.name)

    cwd = str(run_workdir)'''

    # 1) get_fasta_for_insertion.py
    #print("python3",str(get_fasta_py),str(rnafold_path))
    subprocess.run(
        ["python3", str(get_fasta_py), str(rnafold_path), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    intermediate_hairpins = run_workdir / f"{rnafold_path.name}_final_hairpins.txt" #check
    if not intermediate_hairpins.is_file():
        raise FileNotFoundError(f"Expected intermediate not created: {intermediate_hairpins}")

    # 2) SNP_insertion.py <fasta_local.name> <hairpins_file.name>
    subprocess.run(
        ["python3", str(snp_py), str(fasta_path), str(intermediate_hairpins), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    potential_rnafold = run_workdir / f"{intermediate_hairpins.name}_potential_RNAfold.txt" #check
    if not potential_rnafold.is_file():
        raise FileNotFoundError(f"Expected intermediate not created: {potential_rnafold}")

    # 3) Feature extraction su potential_rnafold
    subprocess.run(
        ["python3", str(feat_py), str(potential_rnafold), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    feature_table = run_workdir / f"{potential_rnafold.name}_features_table.tsv" #check
    if not feature_table.is_file():
        raise FileNotFoundError(f"Feature table not created: {feature_table}")

    # 4) Predizione R su feature_table (scrive *_temp_pred_potential.tsv)
    subprocess.run(
        ["Rscript", str(pred_r), str(feature_table), str(run_workdir)], 
        check=True,
        cwd=str(run_workdir)
    )
    potential_pred = run_workdir / f"{potential_rnafold.name}_pred.tsv" #check il suffisso _features_table.tsv è tolto
    if not potential_pred.is_file():
        raise FileNotFoundError(f"Potential prediction not created: {potential_pred}")

    # 5) merge_table.py <base_pred_local> <potential_pred>
    subprocess.run(
        ["python3", str(merge_py), str(base_pred_path), str(potential_pred),str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    classes_table = run_workdir / f"{base_pred_path.name}_temp_classes_table.tsv" #check
    if not classes_table.is_file():
        raise FileNotFoundError(f"Classes table not created: {classes_table}")

    # 6) get_lens.py <rnafold_local>
    subprocess.run(
        ["python3", str(get_lens_py), str(rnafold_path), str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    lens_file = run_workdir / f"{rnafold_path.name}_temp_len_hairpins.txt" #check
    if not lens_file.is_file():
        raise FileNotFoundError(f"Lens file not created: {lens_file}")

    # 7) fisher test
    subprocess.run(
        ["python3", str(fisher_py), str(classes_table), str(lens_file),str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    dr_class = run_workdir / f"{classes_table.name}_temp_DR_class.tsv" #check
    is_class = run_workdir / f"{classes_table.name}_temp_IS_class.tsv" #check
    if not dr_class.is_file() or not is_class.is_file():
        raise FileNotFoundError("DR/IS class tables not created as expected.")

    # 8) first classes
    subprocess.run(
        ["python3", str(first_py), str(classes_table),str(run_workdir)],
        check=True,
        cwd=str(run_workdir)
    )
    classes_mirnas = run_workdir / f"{classes_table.name}_temp_classes_mirna.txt" #check
    if not classes_mirnas.is_file():
        raise FileNotFoundError(f"Classes mirnas file not created: {classes_mirnas}")

    # 9) final classes -> out_path (assoluto, quindi va dove lanci il comando)
    subprocess.run(
        ["python3", str(final_py), str(dr_class), str(is_class), str(classes_mirnas), str(out_path)],
        check=True,
        cwd=str(run_workdir)
    )

    print(f"[OK] Run id: {run_id}")
    print(f"[OK] Workdir: {run_workdir}")
    print(f"[OK] Output classes: {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="miRF classes wrapper (Strada A)")
    parser.add_argument("rnafold", type=str, help="RNAfold file (relative to main/input/ or absolute)")
    parser.add_argument("fasta", type=str, help="FASTA file (relative to main/input/ or absolute)")
    parser.add_argument("base_pred", type=str, help="Prediction file from previous step (relative to cwd OR main/input/ OR absolute)")
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

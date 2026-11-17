#!/usr/bin/env python3
"""Runsheet ERCC_mix gate + run combined_ercc_analysis.ipynb (python → R → python)."""

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

MIX_COL = "ERCC_mix"
EXPECTED_MIXES = ("Mix 1", "Mix 2")
USES_KERNEL = re.compile(r"Uses \*\*(python|R)\*\* kernel", re.I)
PY_KERNEL = "python3"
R_KERNEL = "ir"


def check_runsheet_mixes(runsheet):
    with open(runsheet, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if MIX_COL not in (reader.fieldnames or []):
            return (
                f"Runsheet {Path(runsheet).name} has no '{MIX_COL}' column; it is required when has_ERCC is true. "
                f"Columns: {list(reader.fieldnames or [])}"
            )
        rows = list(reader)
    blank = [r.get("Sample Name") for r in rows if not (r.get(MIX_COL) or "").strip()]
    if blank:
        return f"Runsheet column '{MIX_COL}' is blank for samples: {blank}"
    observed = {(r.get(MIX_COL) or "").strip() for r in rows}
    missing = [m for m in EXPECTED_MIXES if m not in observed]
    extra = sorted(observed - set(EXPECTED_MIXES))
    if not missing and not extra:
        return None
    bits = []
    if missing:
        bits.append(f"missing {missing}")
    if extra:
        bits.append(f"unexpected values {extra}")
    return (
        f"Runsheet column '{MIX_COL}' must be 'Mix 1' or 'Mix 2': "
        f"{'; '.join(bits)}. Observed: {sorted(observed)}"
    )


def write_error(msg, assay_suffix):
    Path(f"ERCC_analysis_error{assay_suffix}.txt").write_text(msg.rstrip() + "\n", encoding="utf-8")
    print(msg, flush=True)


def _cell_text(source):
    return "".join(source) if isinstance(source, list) else (source or "")


def kernel_from_markdown(text):
    m = USES_KERNEL.search(text)
    if not m:
        return None
    return R_KERNEL if m.group(1).lower() == "r" else PY_KERNEL


def _py_assign(src, name, value):
    return re.sub(
        rf"^(\s*){re.escape(name)}\s*=\s*.*$",
        lambda m: f"{m.group(1)}{name} = {repr(value)}",
        src,
        flags=re.M,
    )


def _r_assign(src, name, value):
    return re.sub(
        rf"^(\s*){re.escape(name)}\s*<-\s*.*$",
        lambda m: f"{m.group(1)}{name} <- {json.dumps(value)}",
        src,
        flags=re.M,
    )


def patch_cell_source(src, ctx):
    src = _py_assign(src, "accession", ctx["accession"])
    src = _py_assign(src, "runsheetPath", ctx["runsheet"])
    src = _py_assign(src, "UnnormalizedCountsPath", ctx["counts"])
    src = _py_assign(src, "GENE_ID_PREFIX", ctx["gene_id_prefix"])
    src = _py_assign(src, "assay_suffix", ctx["assay_suffix"])
    src = _r_assign(src, "assay_suffix", ctx["assay_suffix"])
    if ctx.get("ercc_table"):
        src = _py_assign(src, "ercc_url", ctx["ercc_table"])
    return src


def _iopub_to_nb(msg):
    kind = msg["header"]["msg_type"]
    content = msg["content"]
    if kind == "stream":
        return {"output_type": "stream", "name": content.get("name", "stdout"), "text": content.get("text", "")}
    if kind in ("display_data", "execute_result"):
        out = {"output_type": kind, "data": content.get("data", {}), "metadata": content.get("metadata") or {}}
        if kind == "execute_result":
            out["execution_count"] = content.get("execution_count")
        return out
    if kind == "error":
        return {
            "output_type": "error",
            "ename": content.get("ename", "Error"),
            "evalue": content.get("evalue", ""),
            "traceback": content.get("traceback", []),
        }
    return None


def _start_kernel(name, timeout):
    from jupyter_client import KernelManager

    print(f"starting kernel {name}", flush=True)
    km = KernelManager(kernel_name=name)
    km.start_kernel()
    kc = km.client()
    kc.start_channels()
    try:
        kc.wait_for_ready(timeout=timeout)
        if name == PY_KERNEL:
            kc.execute_interactive("%matplotlib inline", timeout=timeout)
    except Exception:
        _stop_kernel(km, kc)
        raise
    return km, kc


def _stop_kernel(km, kc):
    if kc is not None:
        try:
            kc.stop_channels()
        except Exception:
            pass
    if km is not None:
        try:
            km.shutdown_kernel(now=True)
        except Exception:
            pass


def _run_cell(kc, code, timeout):
    outputs = []

    def hook(msg):
        out = _iopub_to_nb(msg)
        if out is not None:
            outputs.append(out)

    reply = kc.execute_interactive(code, timeout=timeout, output_hook=hook)
    if reply["content"]["status"] != "ok":
        err = reply["content"]
        raise RuntimeError(f"{err.get('ename', 'error')}: {err.get('evalue', '')}")
    return outputs


def execute_notebook(nb_path, out_path, ctx, timeout=600):
    nb = json.loads(Path(nb_path).read_text(encoding="utf-8"))
    kernel = PY_KERNEL
    km = kc = None
    current = None
    n = 0
    err = None
    try:
        for cell in nb.get("cells", []):
            src = _cell_text(cell.get("source"))
            if cell.get("cell_type") == "markdown":
                marked = kernel_from_markdown(src)
                if marked:
                    kernel = marked
                continue
            if cell.get("cell_type") != "code":
                continue
            patched = patch_cell_source(src, ctx)
            cell["source"] = patched
            if current != kernel:
                _stop_kernel(km, kc)
                km = kc = None
                km, kc = _start_kernel(kernel, timeout)
                current = kernel
            n += 1
            print(f"executing cell {n} on {kernel}", flush=True)
            try:
                cell["outputs"] = _run_cell(kc, patched, timeout)
                cell["execution_count"] = n
            except Exception as e:
                err = f"Notebook cell {n} ({kernel}) failed: {e}"
                break
    finally:
        Path(out_path).write_text(json.dumps(nb) + "\n", encoding="utf-8")
        _stop_kernel(km, kc)
    return err


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--accession", default="")
    p.add_argument("--runsheet", default="")
    p.add_argument("--assay-suffix", default="_GLbulkRNAseq")
    p.add_argument("--counts", default="")
    p.add_argument("--gene-id-prefix", default="")
    p.add_argument("--notebook", default="combined_ercc_analysis.ipynb")
    p.add_argument("--ercc-table", default="")
    p.add_argument("--cell-timeout", type=int, default=600)
    args = p.parse_args()

    def _empty(v):
        return v if v not in ("", "[]", "null") else ""

    runsheet = _empty(args.runsheet)
    counts = _empty(args.counts)
    ercc_table = _empty(args.ercc_table)
    if not runsheet:
        mix_error = "No runsheet; cannot read ERCC_mix."
    else:
        try:
            mix_error = check_runsheet_mixes(runsheet)
        except Exception as e:
            mix_error = str(e) or type(e).__name__

    lines = [
        f"accession={args.accession}",
        f"runsheetPath={runsheet}",
        f"assay_suffix={args.assay_suffix}",
        f"UnnormalizedCountsPath={counts}",
        f"GENE_ID_PREFIX={args.gene_id_prefix}",
        f"mix_ok={mix_error is None}",
    ]
    Path("ercc_inputs.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)

    if mix_error:
        write_error(mix_error, args.assay_suffix)
        return

    nb_path = Path(args.notebook)
    if not nb_path.is_file():
        write_error(f"notebook not found: {nb_path}", args.assay_suffix)
        return
    if counts and not Path(counts).is_file():
        write_error(f"counts file not found: {counts}", args.assay_suffix)
        return

    ctx = {
        "accession": args.accession,
        "runsheet": runsheet,
        "assay_suffix": args.assay_suffix,
        "counts": counts,
        "gene_id_prefix": args.gene_id_prefix,
        "ercc_table": ercc_table,
    }
    out_nb = Path(f"combined_ercc_analysis{args.assay_suffix}.ipynb")
    html_name = f"ERCC_analysis{args.assay_suffix}.html"
    err = execute_notebook(nb_path, out_nb, ctx, timeout=args.cell_timeout)
    if err:
        write_error(err, args.assay_suffix)
        return
    try:
        subprocess.run(
            ["jupyter", "nbconvert", "--to", "html", f"--output={html_name}", str(out_nb)],
            check=True,
        )
    except Exception as e:
        write_error(f"nbconvert failed: {e}", args.assay_suffix)
        return
    print(f"wrote {out_nb} and {html_name}", flush=True)


if __name__ == "__main__":
    main()

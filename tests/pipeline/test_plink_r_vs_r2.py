"""
Proves that PLINK is invoked with the signed `--r` statistic (not `--r2`),
and that signed, negative LD values survive the parsing/merge logic intact
(in particular, they are not confused with the "missing LD" sentinel).

subprocess.run is faked so this does not depend on a working PLINK binary;
it still exercises the real argument-building and output-parsing code in
app/colocalization/plink.py.
"""
import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from app.colocalization.plink import plink_ld_pairwise, plink_ldmat


@pytest.fixture
def fake_plink_run(monkeypatch, tmp_path):
    """Capture the argv passed to PLINK, and fake its output files."""
    captured = {}

    def fake_run(args, stdout=None, stderr=None, **kwargs):
        captured["args"] = args
        out_prefix = args[args.index("--out") + 1]

        if "square" in args:
            # plink_ldmat: write a fake .bim (15 SNPs) and a signed, square .ld matrix
            bim_df = pd.DataFrame(
                {
                    0: [1] * 3,
                    1: ["chr1:100", "chr1:200", "chr1:300"],
                    2: [0, 0, 0],
                    3: [100, 200, 300],
                    4: ["A", "A", "A"],
                    5: ["G", "G", "G"],
                }
            )
            bim_df.to_csv(out_prefix + ".bim", sep="\t", header=False, index=False)
            # A real negative correlation (r=-1) on the off-diagonal: only possible
            # with --r (signed), impossible with --r2 (always >= 0).
            ld_matrix = np.array([[1.0, -1.0, 0.3], [-1.0, 1.0, -0.3], [0.3, -0.3, 1.0]])
            np.savetxt(out_prefix + ".ld", ld_matrix, delimiter="\t")
        else:
            # plink_ld_pairwise: write a fake pairwise .ld report
            pairwise_df = pd.DataFrame(
                {
                    "CHR_A": [1, 1],
                    "BP_A": [200, 200],
                    "SNP_A": ["chr1:200", "chr1:200"],
                    "CHR_B": [1, 1],
                    "BP_B": [100, 300],
                    "SNP_B": ["chr1:100", "chr1:300"],
                    "R": [-1.0, 0.3],
                }
            )
            pairwise_df.to_csv(out_prefix + ".ld", sep="\t", index=False)

        class FakeCompletedProcess:
            returncode = 0
            stdout = b""

        return FakeCompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)
    return captured


def test_plink_ldmat_uses_signed_r(flask_app, fake_plink_run, tmp_path, monkeypatch):
    """plink_ldmat must invoke PLINK with --r (signed), not --r2, and the
    resulting matrix must preserve real negative (signed) LD values."""
    with flask_app.app_context():
        monkeypatch.setattr(
            "app.colocalization.plink.resolve_plink_filepath",
            lambda build, pop, chrom: "fake_plink_path",
        )
        _, ldmat = plink_ldmat(
            build="hg38",
            pop="EUR",
            chrom=1,
            snp_positions=[100, 200, 300],
            outfilename=str(tmp_path / "mat"),
        )

    args = fake_plink_run["args"]
    assert "--r" in args, "PLINK must be called with --r (signed), the whole point of this change"
    assert "--r2" not in args, "PLINK must NOT be called with --r2 anymore"

    assert (ldmat < 0).any(), "a real negative (signed) r value must survive parsing"
    assert ldmat[0, 1] == -1.0


def test_plink_ld_pairwise_returns_signed_r_column(
    flask_app, fake_plink_run, tmp_path, monkeypatch
):
    """plink_ld_pairwise must produce an 'R' column (not 'R2'), and a real
    r=-1 must not be conflated with the "missing LD" sentinel."""
    with flask_app.app_context():
        monkeypatch.setattr(
            "app.colocalization.plink.resolve_plink_filepath",
            lambda build, pop, chrom: "fake_plink_path",
        )
        monkeypatch.setattr(
            "app.colocalization.plink.find_plink_1kg_overlap",
            lambda plink_filepath, snp_positions, snp_pvalues: pd.DataFrame(
                {"pos": snp_positions, "p": snp_pvalues, 3: snp_positions}
            ),
        )
        df, lead_snp_position = plink_ld_pairwise(
            build="hg38",
            pop="EUR",
            chrom=1,
            snp_positions=[100, 200, 300, 999],  # 999 has no LD info at all
            snp_pvalues=[1e-8, 0.3, 0.5, 0.5],
            outfilename=str(tmp_path / "pw"),
        )

    args = fake_plink_run["args"]
    assert "--r" in args and "--r2" not in args

    assert "R" in df.columns
    by_pos = df.set_index("pos")["R"]

    # Real, signed r=-1 (perfect negative LD) must be reported as -1, not
    # treated as "missing" (which used to be sentinel-encoded as -1 too).
    assert by_pos[100] == -1.0
    assert by_pos[300] == 0.3
    # A SNP with no LD data at all must be genuinely missing (NaN), which is
    # distinguishable from a real r=-1 -- the bug this change fixes.
    assert pd.isna(by_pos[999])

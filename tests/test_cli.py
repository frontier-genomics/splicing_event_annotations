import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run_cli(input_path: Path, output_path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "splicing_event_annotator.main",
            str(input_path),
            str(output_path),
            "--dataset",
            "refseq",
            "--genome",
            "hg38",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        check=False,
    )


def _read_tsv_rows(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_cli_generates_annotations(tmp_path: Path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    input_path.write_text(
        "\t".join(
            [
                "chromosome",
                "start",
                "end",
                "strand",
                "transcript",
                "event_type",
            ]
        )
        + "\n"
        + "\t".join(["chr9", "34646787", "34647088", "+", "NM_000155.4", "sj"])
        + "\n"
    )

    result = _run_cli(input_path, output_path)

    assert result.returncode == 0, result.stderr
    rows = _read_tsv_rows(output_path)
    assert rows, "Expected CLI to produce output rows"
    assert rows[0]["event"] == "canonical exon 1-2 splicing"
    assert rows[0]["event_type"] == "canonical"
    assert rows[0]["introns"] == "1"


def test_cli_rejects_missing_columns(tmp_path: Path):
    input_path = tmp_path / "bad_input.tsv"
    output_path = tmp_path / "output.tsv"
    input_path.write_text(
        "chromosome\tstart\tend\tstrand\ttranscript\n"
        "chr9\t34646787\t34647088\t+\tNM_000155.4\n"
    )

    result = _run_cli(input_path, output_path)

    assert result.returncode != 0
    assert "missing required columns" in result.stderr.lower()


def test_cli_warns_on_overwrite(tmp_path: Path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    input_path.write_text(
        "chromosome\tstart\tend\tstrand\ttranscript\tevent_type\n"
        "chr9\t34646787\t34647088\t+\tNM_000155.4\tsj\n"
    )
    output_path.write_text("existing data")

    result = _run_cli(input_path, output_path)

    assert result.returncode == 0, result.stderr
    assert "already exists" in result.stderr.lower()
    rows = _read_tsv_rows(output_path)
    assert rows[0]["event_type"] == "canonical"

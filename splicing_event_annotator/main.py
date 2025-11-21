import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate
from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList
from splicing_event_annotator.tsvAnnotateProcessor import TsvAnnotate

logger = logging.getLogger(__name__)


REQUIRED_COLUMNS: Sequence[str] = (
    "chromosome",
    "start",
    "end",
    "strand",
    "transcript",
    "event_type",
)


def run_workflow(
    input_data,
    dataset: str = "refseq",
    genome: str = "hg38",
    tsv: bool = False,
    columns: List[int] | None = None,
    output_file: str | Path | None = None,
):
    """Run the annotation workflow."""

    if columns is None:
        columns = [0, 1, 2, 3, 4, 5]

    if tsv:
        tsv_processor = TsvAnnotate(input_data, dataset, columns)
        input_data = tsv_processor.tsv

    if isinstance(input_data, list):
        processor = EventAnnotateList(input_data, dataset, genome)
        annotations = processor.process()
    else:
        processor = EventAnnotate(
            input_data["chrom"],
            input_data["start"],
            input_data["end"],
            input_data["strand"],
            input_data["transcript"],
            input_data["type"],
        )
        annotations = processor.process(dataset, genome)

    if output_file:
        write_output(annotations, output_file)

    return annotations


def write_output(output: Iterable[dict], output_file: str | Path):
    output_path = Path(output_file)
    rows = list(output)
    if not rows:
        raise ValueError("No annotations were generated for the provided input")

    keys = list(rows[0].keys())

    with output_path.open("w", newline="") as file:
        dict_writer = csv.DictWriter(file, keys, delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(rows)


def _validate_required_columns(fieldnames: Sequence[str]) -> None:
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing_columns:
        raise ValueError(
            f"Input file missing required columns: {', '.join(missing_columns)}"
        )


def _load_input_rows(input_file: Path) -> list[dict]:
    with input_file.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input file is empty or missing a header row")
        _validate_required_columns(reader.fieldnames)

        events: list[dict] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid integer value for start/end on line {line_number}"
                ) from exc

            event_type = row["event_type"].strip()
            if not event_type:
                raise ValueError(f"Missing event_type on line {line_number}")

            strand = row["strand"].strip()
            if not strand:
                raise ValueError(f"Missing strand on line {line_number}")

            events.append(
                {
                    "chrom": row["chromosome"],
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "transcript": row["transcript"],
                    "type": event_type.lower(),
                }
            )

    if not events:
        raise ValueError("Input file does not contain any rows to process")

    return events


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate splicing events from a TSV file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", help="Path to a TSV file with splicing events")
    parser.add_argument("output", help="Where to write the annotated TSV")
    parser.add_argument(
        "--dataset",
        default="refseq",
        choices=["refseq", "ensembl"],
        help="Annotation dataset to use",
    )
    parser.add_argument(
        "--genome",
        default="hg38",
        help="Genome assembly matching the annotation dataset",
    )

    if argv is None:
        argv = sys.argv[1:]

    if len(argv) == 0:
        parser.print_help()
        raise SystemExit(0)

    args = parser.parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists():
        parser.error(f"Input file not found: {input_path}")

    try:
        args.events = _load_input_rows(input_path)
    except ValueError as exc:
        parser.error(str(exc))

    args.output_path = Path(args.output)
    return args


def main(argv: Sequence[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO)
    args = parse_args(argv)

    if args.output_path.exists():
        logger.warning("Output file %s already exists and will be overwritten", args.output_path)

    annotations = run_workflow(args.events, args.dataset, args.genome)
    write_output(annotations, args.output_path)


if __name__ == "__main__":
    main()

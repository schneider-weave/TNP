#!/usr/bin/env python3
"""Filter nanobody TXT files with TNP only and write passed sequences.

Input formats supported:
  - One raw sequence per line
  - Optional `label sequence` lines, where the last whitespace-delimited token
    is treated as the sequence
  - FASTA-like headers starting with `>` are ignored

Behavior:
  - Accepts a single TXT file or an input directory of TXT files
  - Runs TNP only, without diversity or nativeness checks
  - Writes passed sequences to the output directory using the same file name
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class InputRecord:
    """A parsed sequence row from the input TXT file."""

    line_no: int
    label: Optional[str]
    sequence: str


def normalize_sequence(raw: str) -> str:
    """Normalize a raw sequence string to uppercase amino-acid text."""
    sequence = re.sub(r"\s+", "", raw).upper().strip()
    invalid = set(sequence) - AMINO_ACIDS
    if invalid:
        raise ValueError(
            f"Invalid amino acid characters: {sorted(invalid)} in sequence {raw!r}"
        )
    return sequence


def parse_input_file(path: Path) -> List[InputRecord]:
    """Parse plain text or simple FASTA-like input into sequence records."""
    records: List[InputRecord] = []
    lines = path.read_text(encoding="utf-8").splitlines()

    for line_no, raw_line in enumerate(lines, start=1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith(";"):
            continue
        if stripped.startswith(">"):
            continue

        parts = stripped.split()
        label = None

        if len(parts) > 1:
            try:
                # Prefer explicit label + sequence rows when whitespace is present.
                label = " ".join(parts[:-1])
                sequence = normalize_sequence(parts[-1])
            except ValueError:
                sequence = normalize_sequence(stripped)
                label = None
        else:
            sequence = normalize_sequence(stripped)

        records.append(InputRecord(line_no=line_no, label=label, sequence=sequence))

    if not records:
        raise ValueError("No sequences found in input file.")

    return records


def write_sequences(path: Path, sequences: Sequence[str]) -> None:
    """Write passed sequences to a text file, one per line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sequences) + ("\n" if sequences else ""), encoding="utf-8")


def _analyze_sequence(sequence: str) -> bool:
    """Return True when the sequence passes TNP thresholds."""
    from metanano.config import DevelopabilityConfig
    from metanano.filters.developability import DevelopabilityFilter

    dev_filter = DevelopabilityFilter(DevelopabilityConfig())
    result = dev_filter.analyze(sequence)
    return result.passed


def filter_passed_sequences(input_path: Path, output_path: Path) -> Tuple[int, int]:
    """Run TNP-only filtering for a single TXT file."""
    records = parse_input_file(input_path)
    passed_sequences: List[str] = []
    for record in records:
        if _analyze_sequence(record.sequence):
            passed_sequences.append(record.sequence)

    write_sequences(output_path, passed_sequences)
    return len(passed_sequences), len(records) - len(passed_sequences)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser."""
    parser = argparse.ArgumentParser(
        description="Filter nanobody TXT files with TNP and write passed sequences."
    )
    parser.add_argument("input", type=Path, help="Input TXT file or directory")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output TXT file or directory",
    )
    return parser


def iter_input_txt_files(input_path: Path) -> Iterable[Path]:
    """Yield TXT files from a file or directory input."""
    if input_path.is_file():
        yield input_path
        return

    if input_path.is_dir():
        for path in sorted(input_path.iterdir()):
            if path.is_file() and path.suffix.lower() == ".txt":
                yield path
        return

    raise ValueError(f"Input path does not exist: {input_path}")


def resolve_output_path(input_file: Path, output_root: Path) -> Path:
    """Preserve the input file name under the output location."""
    if input_file.is_dir():
        if output_root.suffix.lower() == ".txt" and not output_root.is_dir():
            raise ValueError(
                "When input is a directory, output must be a directory path."
            )
        output_root.mkdir(parents=True, exist_ok=True)
        return output_root / input_file.name

    if output_root.exists() and output_root.is_dir():
        output_root.mkdir(parents=True, exist_ok=True)
        return output_root / input_file.name

    if output_root.suffix.lower() == ".txt":
        output_root.parent.mkdir(parents=True, exist_ok=True)
        return output_root

    output_root.mkdir(parents=True, exist_ok=True)
    return output_root / input_file.name


def main() -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    input_path: Path = args.input

    try:
        output_path: Path = args.output
        total_passed = 0
        total_failed = 0

        for input_file in iter_input_txt_files(input_path):
            file_output = resolve_output_path(input_file, output_path)
            passed_count, failed_count = filter_passed_sequences(
                input_path=input_file,
                output_path=file_output,
            )
            total_passed += passed_count
            total_failed += failed_count
    except Exception as exc:  # noqa: BLE001
        parser.error(str(exc))
        return 2

    print(f"Passed: {total_passed}, Failed: {total_failed}, Output: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

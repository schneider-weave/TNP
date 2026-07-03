from pathlib import Path

import pytest

from scripts.filter_nanobody_txt import (
    filter_passed_sequences,
    iter_input_txt_files,
    normalize_sequence,
    parse_input_file,
    resolve_output_path,
    write_sequences,
)


def test_normalize_sequence_removes_whitespace_and_uppercases() -> None:
    assert normalize_sequence(" evqlvesggg lvqpg ") == "EVQLVESGGGLVQPG"


def test_parse_input_file_supports_comments_and_labels(tmp_path: Path) -> None:
    input_path = tmp_path / "input.txt"
    input_path.write_text(
        "\n".join(
            [
                "# comment line",
                ">fasta_header",
                "seq_a EVQLVESGGGLVQPG",
                "QVQLVESGGGLVQPG",
                "",
            ]
        ),
        encoding="utf-8",
    )

    records = parse_input_file(input_path)

    assert len(records) == 2
    assert records[0].label == "seq_a"
    assert records[0].sequence == "EVQLVESGGGLVQPG"
    assert records[1].label is None
    assert records[1].sequence == "QVQLVESGGGLVQPG"


def test_write_sequences_writes_trailing_newline(tmp_path: Path) -> None:
    output_path = tmp_path / "passed.txt"
    write_sequences(output_path, ["SEQ1", "SEQ2"])

    assert output_path.read_text(encoding="utf-8") == "SEQ1\nSEQ2\n"


def test_normalize_sequence_rejects_invalid_characters() -> None:
    with pytest.raises(ValueError):
        normalize_sequence("EVQLV3S")


def test_iter_input_txt_files_yields_txt_files_only(tmp_path: Path) -> None:
    (tmp_path / "a.txt").write_text("AAA", encoding="utf-8")
    (tmp_path / "b.fasta").write_text("BBB", encoding="utf-8")
    (tmp_path / "c.txt").write_text("CCC", encoding="utf-8")

    files = list(iter_input_txt_files(tmp_path))

    assert [path.name for path in files] == ["a.txt", "c.txt"]


def test_resolve_output_path_preserves_filename_for_directory_mode(
    tmp_path: Path,
) -> None:
    input_file = tmp_path / "input.txt"
    input_file.write_text("AAA", encoding="utf-8")
    output_dir = tmp_path / "output"

    resolved = resolve_output_path(input_file, output_dir)

    assert resolved == output_dir / "input.txt"


def test_filter_passed_sequences_writes_only_passed_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_file = tmp_path / "sample.txt"
    input_file.write_text("PASSSEQ\nFAILSEQ\n", encoding="utf-8")
    output_file = tmp_path / "output.txt"

    monkeypatch.setattr(
        "scripts.filter_nanobody_txt._analyze_sequence",
        lambda sequence: sequence == "PASSSEQ",
    )

    passed_count, failed_count = filter_passed_sequences(input_file, output_file)

    assert passed_count == 1
    assert failed_count == 1
    assert output_file.read_text(encoding="utf-8") == "PASSSEQ\n"

"""Build a VSEARCH SINTAX reference from official Greengenes2 artifacts."""

from __future__ import annotations

import argparse
import csv
import io
import zipfile
from pathlib import Path


RANKS = set("dpcofgs")


def _artifact_member(artifact: zipfile.ZipFile, suffix: str) -> str:
    matches = [name for name in artifact.namelist() if name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {suffix!r} file in the QIIME 2 artifact, "
            f"found {len(matches)}."
        )
    return matches[0]


def _sintax_taxonomy(lineage: str) -> str:
    converted = []
    for rank in lineage.split(";"):
        rank = rank.strip()
        if "__" not in rank:
            raise ValueError(f"Malformed Greengenes2 rank: {rank!r}")
        prefix, name = rank.split("__", 1)
        if prefix not in RANKS:
            raise ValueError(f"Unsupported Greengenes2 rank: {rank!r}")
        if not name:
            continue
        if not name.isascii() or "," in name or ";" in name:
            raise ValueError(f"Taxon name is not SINTAX-compatible: {name!r}")
        converted.append(f"{prefix}:{name}")
    if not converted:
        raise ValueError(f"Lineage has no named ranks: {lineage!r}")
    return ",".join(converted)


def _load_taxonomy(artifact: zipfile.ZipFile) -> dict[str, str]:
    member = _artifact_member(artifact, "/data/taxonomy.tsv")
    with artifact.open(member) as raw:
        rows = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"), delimiter="\t")
        if rows.fieldnames != ["Feature ID", "Taxon"]:
            raise ValueError(
                "Expected Greengenes2 taxonomy columns: Feature ID and Taxon."
            )
        taxonomy = {}
        for row in rows:
            feature_id = row["Feature ID"]
            if feature_id in taxonomy:
                raise ValueError(f"Duplicate taxonomy feature ID: {feature_id}")
            taxonomy[feature_id] = _sintax_taxonomy(row["Taxon"])
    return taxonomy


def prepare_reference(sequence_qza, taxonomy_qza, output_fasta) -> int:
    """Convert matching Greengenes2 sequence/taxonomy QZAs to SINTAX FASTA."""
    sequence_qza = Path(sequence_qza)
    taxonomy_qza = Path(taxonomy_qza)
    output_fasta = Path(output_fasta)

    if output_fasta.exists():
        raise FileExistsError(f"Output already exists: {output_fasta}")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    try:
        with zipfile.ZipFile(taxonomy_qza) as tax_artifact:
            taxonomy = _load_taxonomy(tax_artifact)

        with zipfile.ZipFile(sequence_qza) as seq_artifact:
            member = _artifact_member(seq_artifact, "/data/dna-sequences.fasta")
            with seq_artifact.open(member) as raw, output_fasta.open(
                    "w", encoding="ascii", newline="\n") as output:
                sequences = io.TextIOWrapper(raw, encoding="ascii")
                current_id = None
                for line in sequences:
                    if line.startswith(">"):
                        current_id = line[1:].strip()
                        taxon = taxonomy.pop(current_id, None)
                        if taxon is None:
                            raise ValueError(
                                "Sequence has no matching taxonomy or is duplicated: "
                                f"{current_id}"
                            )
                        output.write(f">{current_id};tax={taxon};\n")
                        count += 1
                    elif current_id is None and line.strip():
                        raise ValueError("Sequence artifact does not contain valid FASTA.")
                    else:
                        output.write(line)
        if count == 0:
            raise ValueError("Sequence artifact contains no FASTA records.")
    except Exception:
        if output_fasta.exists():
            output_fasta.unlink()
        raise

    return count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert official Greengenes2 QZA artifacts to VSEARCH SINTAX FASTA."
    )
    parser.add_argument("sequence_qza", help="GG2 backbone V4 or full-length sequence QZA")
    parser.add_argument("taxonomy_qza", help="Matching GG2 backbone taxonomy QZA")
    parser.add_argument("output_fasta", help="New SINTAX-compatible FASTA path")
    args = parser.parse_args()

    try:
        count = prepare_reference(
            args.sequence_qza, args.taxonomy_qza, args.output_fasta
        )
    except (FileExistsError, OSError, ValueError, zipfile.BadZipFile) as error:
        parser.error(str(error))
    print(f"Wrote {count:,} Greengenes2 reference sequences to {args.output_fasta}")


if __name__ == "__main__":
    main()

import gzip
from pathlib import Path


ENA_SINGLE_CONTIG_HELP = (
    "ENA accepts a one-sequence MAG as a chromosome assembly only when the "
    "sequence represents a fully assembled chromosome. Incomplete one-sequence "
    "MAGs require an exception from the ENA helpdesk."
)


def inspect_gzipped_fasta(path):
    """Read the basic information needed from a gzipped FASTA."""
    fasta_path = Path(path)
    sequence_count = 0
    first_sequence_name = ""
    first_raw_header = ""

    try:
        with gzip.open(fasta_path, "rt", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.startswith(">"):
                    continue

                raw_header = line[1:].strip()
                if not raw_header:
                    raise ValueError(
                        f"{fasta_path.name}: empty FASTA header at line {line_number}."
                    )

                sequence_name = raw_header.split()[0]
                sequence_count += 1
                if sequence_count == 1:
                    first_sequence_name = sequence_name
                    first_raw_header = raw_header
    except (OSError, UnicodeError) as exc:
        raise ValueError(
            f"{fasta_path.name}: could not read gzip-compressed FASTA ({exc})."
        ) from exc

    if sequence_count == 0:
        raise ValueError(f"{fasta_path.name}: no FASTA records were found.")

    return {
        "sequence_count": sequence_count,
        "first_sequence_name": first_sequence_name,
        "first_raw_header": first_raw_header,
    }


def inspect_fasta_collection(fasta_map):
    """Inspect all uploaded FASTA files."""
    return {
        str(sample_name): inspect_gzipped_fasta(path)
        for sample_name, path in fasta_map.items()
    }


def single_contig_samples(inspections):
    """Return samples with only one sequence."""
    return sorted(
        sample_name
        for sample_name, inspection in inspections.items()
        if inspection["sequence_count"] == 1
    )


def require_single_contig_confirmation(inspections, confirmed):
    """Require confirmation before treating one sequence as a chromosome."""
    samples = single_contig_samples(inspections)
    if samples and not confirmed:
        joined = ", ".join(samples)
        raise ValueError(
            f"Single-contig MAG(s) detected: {joined}. {ENA_SINGLE_CONTIG_HELP} "
            "Confirm that these sequences are complete chromosomes before submission."
        )
    return samples


def chromosome_list_text(inspection, chromosome_name="1"):
    """Build the ENA chromosome-list row for a one-sequence assembly."""
    if inspection["sequence_count"] != 1:
        raise ValueError(
            "A chromosome list can be generated only for a FASTA containing "
            "exactly one sequence."
        )

    # ENA expects the identifier before the first space, not the full header.
    try:
        inspection["first_sequence_name"].encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(
            "The FASTA sequence name must contain only ASCII characters for an "
            "ENA chromosome-list submission."
        ) from exc

    if chromosome_name != "1":
        raise ValueError("MAG2ENA uses the ENA-safe chromosome name '1'.")

    return f"{inspection['first_sequence_name']}\t{chromosome_name}\tchromosome\n"

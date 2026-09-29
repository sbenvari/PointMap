#!/usr/bin/env python3

import argparse
from Bio import AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Call amino-acid substitutions relative to a reference sequence."
    )
    parser.add_argument("--alignment", required=True, help="Aligned protein FASTA")
    parser.add_argument("--output", required=True, help="Output TSV file")
    return parser.parse_args()


def main():
    args = parse_args()
    alignment = AlignIO.read(args.alignment, "fasta")

    ref_index = None
    for i, record in enumerate(alignment):
        if record.id.lower() in {"reference", "ref"}:
            ref_index = i
            break

    if ref_index is None:
        raise ValueError("Reference sequence not found in alignment.")

    ref_seq = alignment[ref_index].seq

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("Sample\tMutations\n")

        for i, record in enumerate(alignment):
            if i == ref_index:
                continue

            mutations = []
            for pos, (ref_aa, sample_aa) in enumerate(
                zip(ref_seq, record.seq), start=1
            ):
                if ref_aa != sample_aa and ref_aa != "-" and sample_aa != "-":
                    mutations.append(f"{ref_aa}{pos}{sample_aa}")

            out.write(
                f"{record.id}\t{','.join(mutations) if mutations else 'No mutations'}\n"
            )


if __name__ == "__main__":
    main()

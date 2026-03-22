#!/usr/bin/env bash
# Usage: bash mauve_pseudochromosome.sh alignment.fa [gap_size] [output.fa]
# This script was used on the MCM output file to format it as a typical pseudochromosome

FASTA="${1:?Usage: $0 alignment.fa [gap_size] [output.fa]}"
GAP_SIZE="${2:-100}"
OUTPUT="${3:-$(basename "${FASTA%.*}").pseudochromosome.fa}"

awk -v gap_size="$GAP_SIZE" '
BEGIN {
    gap = ""
    for (i = 0; i < gap_size; i++) gap = gap "N"
    pseudo = ""
    in_contig = 0
    seq = ""
}

/^>/ {
    # Save accumulated contig sequence when we hit a new header
    if (in_contig && seq != "") {
        gsub(/-/, "", seq)
        pseudo = (pseudo == "") ? seq : pseudo gap seq
        seq = ""
    }
    # Determine if the upcoming sequence is a contig or reference
    # by peeking at the header — we flag it after reading the first seq line
    in_contig = 0
    next
}

# First sequence line after a header reveals case
in_contig == 0 && /^[a-z]/ { in_contig = 1; seq = seq $0; next }
in_contig == 0 && /^[A-Z]/ { in_contig = 0; next }

# Accumulate contig lines
in_contig == 1 { seq = seq $0 }

END {
    # Save the last contig if file does not end with a new header
    if (in_contig && seq != "") {
        gsub(/-/, "", seq)
        pseudo = (pseudo == "") ? seq : pseudo gap seq
    }

    print ">pseudochromosome"
    while (length(pseudo) > 60) {
        print substr(pseudo, 1, 60)
        pseudo = substr(pseudo, 61)
    }
    if (length(pseudo) > 0) print pseudo
}
' "$FASTA" > "$OUTPUT"

echo "Output written to: $OUTPUT"

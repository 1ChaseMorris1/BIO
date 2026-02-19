MATCH_SCORE = 1
MISMATCH_PENALTY = -1
GAP_PENALTY = -2


def clean_sequence_text(text):
    #return only uppercase letters from text.
    cleaned_letters = []
    for ch in text.upper():
        if ch.isalpha():
            cleaned_letters.append(ch)
    return "".join(cleaned_letters)


def sequence_from_fasta_record(record_text):
    # extract a sequence from one FASTA record

    sequence_parts = []
    for line in record_text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            continue
        sequence_parts.append(stripped)
    return clean_sequence_text("".join(sequence_parts))


def read_sequence_input(prompt):
    # sequence from user input
    first_line = input(prompt).strip()
    if not first_line:
        return ""

    # If the input starts with > read FASTA lines until there is a blank line
    if first_line.startswith(">"):
        lines = [first_line]
        while True:
            line = input().strip()
            if line == "":
                break
            lines.append(line)
        return sequence_from_fasta_record("\n".join(lines))

    return clean_sequence_text(first_line)


def needleman_wunsch(seq1, seq2):
    # Build and return score matrix
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score = []
    trace = []
    for _ in range(rows):
        score.append([0] * cols)
        trace.append(["" for _ in range(cols)])
    trace[0][0] = "start"

    # First column every step is a gap in sequence 2
    for i in range(1, rows):
        score[i][0] = i * GAP_PENALTY
        trace[i][0] = "gap_up"

    # First row every step is a gap in sequence 1
    for j in range(1, cols):
        score[0][j] = j * GAP_PENALTY
        trace[0][j] = "gap_left"

    for i in range(1, rows):
        for j in range(1, cols):
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            diag_score = MATCH_SCORE if char1 == char2 else MISMATCH_PENALTY

            diag = score[i - 1][j - 1] + diag_score
            up = score[i - 1][j] + GAP_PENALTY
            left = score[i][j - 1] + GAP_PENALTY

            best = max(diag, up, left)
            score[i][j] = best

            # Store exactly one move for one optimal path.
            # Tie-break order: diagonal, then up, then left.
            if diag == best:
                if char1 == char2:
                    trace[i][j] = "match_diag"
                else:
                    trace[i][j] = "mismatch_diag"
            elif up == best:
                trace[i][j] = "gap_up"
            else:
                trace[i][j] = "gap_left"

    return score, trace


def traceback_alignment(seq1, seq2, trace):
    # Follow the traceback matrix from bottom-right and build one optimal alignment
    i = len(seq1)
    j = len(seq2)
    aligned1 = []
    aligned2 = []

    while i > 0 or j > 0:
        move = trace[i][j]

        if i > 0 and j > 0 and move in ("match_diag", "mismatch_diag"):
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and move == "gap_up":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2)


def traceback_path_cells(seq1, seq2, trace):
    # Return matrix coordinates that belong to the single chosen traceback path
    i = len(seq1)
    j = len(seq2)
    path = {(i, j)}

    while i > 0 or j > 0:
        move = trace[i][j]

        if i > 0 and j > 0 and move in ("match_diag", "mismatch_diag"):
            i -= 1
            j -= 1
        elif i > 0 and move == "gap_up":
            i -= 1
        else:
            j -= 1
        path.add((i, j))

    return path


def traceback_labels(move, on_path=False):
    # Convert one internal move name to a short label for matrix
    if not move or move == "start":
        label = "start"
    elif move in ("match_diag", "mismatch_diag"):
        label = "D"
    elif move == "gap_up":
        label = "U"
    else:
        label = "L"

    return f"{label}+" if on_path else label


def print_matrix(seq1, seq2, score, trace):
    #Print matrix as a table with score and traceback labels in each cell
    path_cells = traceback_path_cells(seq1, seq2, trace)
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)

    table = [[" "] + col_labels]
    for i, row_label in enumerate(row_labels):
        row = [row_label]
        for j in range(len(col_labels)):
            row.append(f"{score[i][j]}:{traceback_labels(trace[i][j], (i, j) in path_cells)}")
        table.append(row)

    widths = []
    for col in range(len(table[0])):
        widths.append(max(len(row[col]) for row in table))

    def border_line():
        return "+" + "+".join("-" * (w + 2) for w in widths) + "+"

    def format_row(row):
        cells = [f" {cell.ljust(widths[i])} " for i, cell in enumerate(row)]
        return "|" + "|".join(cells) + "|"

    print(border_line())
    for row in table:
        print(format_row(row))
        print(border_line())


def alignment_stats(aligned1, aligned2):
    #count identities mismatches and gaps 
    identities = 0
    mismatches = 0
    gaps = 0

    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            identities += 1
        else:
            mismatches += 1

    return identities, mismatches, gaps


def main():
    print("Enter sequence 1:")
    seq1 = read_sequence_input("Sequence 1: ")

    print("\nEnter sequence 2:")
    seq2 = read_sequence_input("Sequence 2: ")

    if not seq1 or not seq2:
        print("Error: both inputs must contain letters.")
        raise SystemExit(1)

    score, trace = needleman_wunsch(seq1, seq2)
    aligned1, aligned2 = traceback_alignment(seq1, seq2, trace)
    identities, mismatches, gaps = alignment_stats(aligned1, aligned2)

    print(f"match = {MATCH_SCORE}")
    print(f"mismatch = {MISMATCH_PENALTY}")
    print(f"gap = {GAP_PENALTY}")

    print("\nSequence 1:", seq1)
    print("Sequence 2:", seq2)

    print("\nMatrix:")
    print_matrix(seq1, seq2, score, trace)

    print("\nOptimal Global Alignment:")
    print(aligned1)
    print(aligned2)

    print(f"\nAlignment Score: {score[len(seq1)][len(seq2)]}")
    print(f"identities = {identities}")
    print(f"mismatches = {mismatches}")
    print(f"gaps = {gaps}")


if __name__ == "__main__":
    main()

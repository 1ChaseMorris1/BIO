MATCH_SCORE = 1
MISMATCH_PENALTY = -1
GAP_PENALTY = -2


# take raw user input and return only uppercase letters.
def clean_sequence(text):    # keep the letters only and make them uppercase
    upper_text = text.upper()
    cleaned_letters = []
    for ch in upper_text:
        if ch.isalpha():
            cleaned_letters.append(ch)
    cleaned_text = "".join(cleaned_letters)
    return cleaned_text


# build the Needleman-Wunsch score matrix and traceback matrix.
def needleman_wunsch(seq1, seq2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # build empty score and traceback matrices
    score = []
    trace = []
    for _ in range(rows):
        score_row = []
        trace_row = []
        for _ in range(cols):
            score_row.append(0)
            trace_row.append([])
        score.append(score_row)
        trace.append(trace_row)

    # initialize first column and first row with gap penalties
    for i in range(1, rows):
        score[i][0] = i * GAP_PENALTY
        trace[i][0] = ["gap_up"]
    for j in range(1, cols):
        score[0][j] = j * GAP_PENALTY
        trace[0][j] = ["gap_left"]

    # Fill the rest of the matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]

            if char1 == char2:
                diag_add = MATCH_SCORE
            else:
                diag_add = MISMATCH_PENALTY

            diag = score[i - 1][j - 1] + diag_add
            up = score[i - 1][j] + GAP_PENALTY
            left = score[i][j - 1] + GAP_PENALTY

            best = max(diag, up, left)
            score[i][j] = best

            # save all best moves
            moves: list[str] = []
            if diag == best:
                if char1 == char2:
                    moves.append("match_diag")
                else:
                    moves.append("mismatch_diag")
            if up == best:
                moves.append("gap_up")
            if left == best:
                moves.append("gap_left")

            trace[i][j] = moves

    return score, trace

# follow traceback moves to build one optimal global alignment
def traceback_alignment(seq1, seq2, trace):
    i = len(seq1)
    j = len(seq2)
    aligned1: list[str] = []
    aligned2: list[str] = []

    # Walk backward from bottom right to the  top left.
    while i > 0 or j > 0:
        moves = trace[i][j]
        take_diag = "match_diag" in moves or "mismatch_diag" in moves
        take_up = "gap_up" in moves

        if i > 0 and j > 0 and take_diag:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and take_up:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()
    final1 = "".join(aligned1)
    final2 = "".join(aligned2)
    return final1, final2

# convert internal move names into labels.
def move_label(moves):
    if not moves:
        return "start"

    # Show one label only (simple view).
    if "match_diag" in moves:
        return "match"
    if "mismatch_diag" in moves:
        return "mismatch"
    return "gap"

# print the matrix in an ASCII SQL-style table (+ and |).
def print_matrix(seq1, seq2, score, trace):
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)

    table = [[" "] + col_labels]
    for i, row_label in enumerate(row_labels):
        row = [row_label]
        for j in range(len(col_labels)):
            row.append(f"{score[i][j]}:{move_label(trace[i][j])}")
        table.append(row)

    widths = []
    col_count = len(table[0])
    for col in range(col_count):
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


def main():
    seq1 = clean_sequence(input("Enter sequence 1: "))
    seq2 = clean_sequence(input("Enter sequence 2: "))
    if not seq1 or not seq2:
        print("Error: both inputs must contain letters.")
        raise SystemExit(1)

    score, trace = needleman_wunsch(seq1, seq2)
    aligned1, aligned2 = traceback_alignment(seq1, seq2, trace)

    print("Scoring:")
    print(f"  match = {MATCH_SCORE}")
    print(f"  mismatch = {MISMATCH_PENALTY}")
    print(f"  gap = {GAP_PENALTY}")
    print("\nSequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("\nN-W Matrix (score:label):")
    print_matrix(seq1, seq2, score, trace)
    print("\nOptimal Global Alignment:")
    print(aligned1)
    print(aligned2)
    print(f"\nAlignment Score: {score[len(seq1)][len(seq2)]}")


if __name__ == "__main__":
    main()

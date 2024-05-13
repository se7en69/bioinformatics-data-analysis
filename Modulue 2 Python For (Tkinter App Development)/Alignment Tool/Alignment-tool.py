import tkinter as tk

def needleman_wunsch(seq1, seq2):
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    gap_penalty = -1

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = 1 if seq1[i - 1] == seq2[j - 1] else -1
            dp[i][j] = max(dp[i - 1][j - 1] + match_score,
                           dp[i - 1][j] + gap_penalty,
                           dp[i][j - 1] + gap_penalty)

    aligned_seq1, aligned_seq2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + gap_penalty:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
        else:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1

    alignment_score = dp[n][m]
    return aligned_seq1, aligned_seq2, alignment_score

def align_sequences():
    seq1 = seq1_entry.get()
    seq2 = seq2_entry.get()
    aligned_seq1, aligned_seq2, score = needleman_wunsch(seq1, seq2)
    result_label.config(text=f"Sequence 1: {aligned_seq1}\nSequence 2: {aligned_seq2}\nAlignment Score: {score}")

# Create the main window
root = tk.Tk()
root.title("Needleman-Wunsch Alignment Tool")

# Create input fields for sequences
seq1_label = tk.Label(root, text="Sequence 1:")
seq1_label.pack()
seq1_entry = tk.Entry(root)
seq1_entry.pack()

seq2_label = tk.Label(root, text="Sequence 2:")
seq2_label.pack()
seq2_entry = tk.Entry(root)
seq2_entry.pack()

# Button to trigger the alignment
align_button = tk.Button(root, text="Align Sequences", command=align_sequences)
align_button.pack()

# Label to display the result
result_label = tk.Label(root, text="")
result_label.pack()

# Start the GUI event loop
root.mainloop()

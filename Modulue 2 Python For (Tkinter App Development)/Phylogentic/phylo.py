from tkinter import *
from tkinter import filedialog
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
import matplotlib.pyplot as plt

def generate_phylogenetic_tree():
    input_file_path = filedialog.askopenfilename(title="Select aligned sequences file", filetypes=[("FASTA files", "*.fasta")])
    if not input_file_path:
        return

    output_file_path = filedialog.asksaveasfilename(title="Save phylogenetic tree as", defaultextension=".png", filetypes=[("PNG files", "*.png")])
    if not output_file_path:
        return

    # Read the alignment file
    alignment = AlignIO.read(input_file_path, "fasta")

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct the tree using Neighbor-Joining method
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Draw the tree
    plt.figure(figsize=(12, 8))
    draw(tree, do_show=False)
    plt.axis('off')
    plt.savefig(output_file_path, bbox_inches='tight')
    plt.close()

    # Display success message
    success_label.config(text="Phylogenetic tree generated successfully!")

# Create a GUI window
root = Tk()
root.title("Phylogenetic Tree Generator")
root.geometry("400x200")

# Create a label
label = Label(root, text="Click below to select an aligned sequences file and generate a phylogenetic tree.", wraplength=380)
label.pack(pady=10)

# Create a button to select input file and generate tree
button = Button(root, text="Select File", command=generate_phylogenetic_tree)
button.pack(pady=10)

# Create a label for success message
success_label = Label(root, text="")
success_label.pack(pady=10)

# Run the Tkinter event loop
root.mainloop()

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox

class GenomeAnnotationTool:
    def __init__(self, master):
        self.master = master
        self.master.title("Genome Annotation Tool")
        self.master.geometry("800x600")

        self.label = tk.Label(self.master, text="Select a GFF file:")
        self.label.pack()

        self.select_button = tk.Button(self.master, text="Select File", command=self.select_file)
        self.select_button.pack()

        self.export_button = tk.Button(self.master, text="Export Annotations", command=self.export_annotations, state=tk.DISABLED)
        self.export_button.pack()

        self.gff_file = None
        self.genes = {}

    def select_file(self):
        self.gff_file = filedialog.askopenfilename(filetypes=[("GFF files", "*.gff")])
        if self.gff_file:
            self.export_button.config(state=tk.NORMAL)
            messagebox.showinfo("File Selected", f"Selected file: {self.gff_file}")

    def parse_gff(self):
        with open(self.gff_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                feature_type, start, end, strand, attributes = fields[2], int(fields[3]), int(fields[4]), fields[6], fields[8]
                if feature_type == "gene":
                    gene_id = attributes.split(";")[0].split("=")[1]
                    self.genes[gene_id] = {"start": start, "end": end, "strand": strand, "function": "Unknown", "exons": []}

    def annotate_genes(self):
        for gene_id, gene_info in self.genes.items():
            # Perform gene annotation logic here (e.g., gene function, exon identification)
            gene_info["function"] = "Hypothetical protein"
            gene_info["exons"] = [(gene_info["start"] + 10, gene_info["start"] + 20), (gene_info["start"] + 30, gene_info["start"] + 40)]

    def export_annotations(self):
        if self.gff_file:
            output_file = filedialog.asksaveasfilename(defaultextension=".gff", filetypes=[("GFF files", "*.gff")])
            if output_file:
                self.parse_gff()
                self.annotate_genes()
                with open(output_file, "w") as f:
                    f.write("## Genome Annotation Results\n")
                    for gene_id, gene_info in self.genes.items():
                        f.write(f"{gene_id}\tCustomAnnotation\tgene\t{gene_info['start']}\t{gene_info['end']}\t.\t{gene_info['strand']}\t.\tID={gene_id};Function={gene_info['function']}\n")
                        for exon_start, exon_end in gene_info["exons"]:
                            f.write(f"{gene_id}\tCustomAnnotation\texon\t{exon_start}\t{exon_end}\t.\t{gene_info['strand']}\t.\tParent={gene_id}\n")
                messagebox.showinfo("Export Complete", f"Annotations exported to {output_file}")
        else:
            messagebox.showerror("Error", "No file selected")

if __name__ == "__main__":
    root = tk.Tk()
    app = GenomeAnnotationTool(root)
    root.mainloop()

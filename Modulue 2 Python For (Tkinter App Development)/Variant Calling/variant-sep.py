import tkinter as tk
from tkinter import filedialog, messagebox

# Function to parse the VCF file and separate the variants
def separate_variants(vcf_file):
    # Initialize empty lists to store different variants
    snvs = []
    indels = []
    sv = []
    other = []
    
    # Iterate over each line in the VCF file
    with open(vcf_file, 'r') as file:
        for line in file:
            # Skip header lines
            if line.startswith('#'):
                continue
            
            # Split the line into fields
            fields = line.strip().split('\t')
            
            # Check the variant type
            if len(fields[3]) == len(fields[4]):  # SNV
                snvs.append(line)
            elif len(fields[3]) > len(fields[4]):  # Deletion
                indels.append(line)
            elif len(fields[3]) < len(fields[4]):  # Insertion
                indels.append(line)
            else:  # Structural variant
                sv.append(line)
    
    # Return the separated variants
    return snvs, indels, sv

# Function to handle file selection
def select_file():
    vcf_file = filedialog.askopenfilename(filetypes=[('VCF files', '*.vcf')])
    if vcf_file:
        snvs, indels, sv = separate_variants(vcf_file)
        messagebox.showinfo('Variants Separated', f'Number of SNVs: {len(snvs)}\nNumber of Indels: {len(indels)}\nNumber of Structural Variants: {len(sv)}')

# Create the main application window
app = tk.Tk()
app.title('VCF Variant Separator')

# Create a button to select the VCF file
select_button = tk.Button(app, text='Select VCF File', command=select_file)
select_button.pack(padx=20, pady=20)

# Start the Tkinter event loop
app.mainloop()

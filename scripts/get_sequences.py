import os
from Bio import SeqIO
import re
from collections import defaultdict

def rename_sequences(base_dir):
   """
   Rename sequences in .faa files under .out folders by appending folder and file names.
   """
   for folder in os.listdir(base_dir):
       if folder.endswith(".out"):
           folder_name = folder.replace(".fasta.out", "")
           faa_dir = os.path.join(base_dir, folder, "run_metazoa_odb12", "busco_sequences", "single_copy_busco_sequences")

           if not os.path.exists(faa_dir):
               print(f"Directory not found: {faa_dir}")
               continue

           for faa_file in os.listdir(faa_dir):
               if faa_file.endswith(".faa"):
                   faa_path = os.path.join(faa_dir, faa_file)
                   faa_name = faa_file.replace(".faa", "")

                   new_records = []
                   for record in SeqIO.parse(faa_path, "fasta"):
                       new_id = f"{folder_name}_{faa_name}"
                       record.id = new_id
                       record.description = ""
                       new_records.append(record)

                   SeqIO.write(new_records, faa_path, "fasta")
                   print(f"Processed file: {faa_path}")

   print("All faa sequences are renamed!")

def merge_sequences(base_dir, output_dir):
   """
   Merge renamed .faa files into a single file for each species (folder).
   """
   os.makedirs(output_dir, exist_ok=True)

   for folder in os.listdir(base_dir):
       if folder.endswith(".out"):
           source_dir = os.path.join(base_dir, folder, "run_metazoa_odb12/busco_sequences/single_copy_busco_sequences")
           if not os.path.exists(source_dir):
               print(f"Warning: {source_dir} does not exist, skipping.")
               continue

           output_filename = f"{folder.replace('.fasta.out', '')}.faa"
           output_file = os.path.join(output_dir, output_filename)

           with open(output_file, "w") as outfile:
               for faa_file in os.listdir(source_dir):
                   if faa_file.endswith(".faa"):
                       faa_path = os.path.join(source_dir, faa_file)
                       for record in SeqIO.parse(faa_path, "fasta"):
                           SeqIO.write(record, outfile, "fasta")

           print(f"Created {output_file}.")

   print("All faa files are merged!")

def group_by_gene_family(input_dir, output_dir):
   """
   Group sequences by gene family and save into separate .faa files.
   """
   os.makedirs(output_dir, exist_ok=True)
   gene_family_dict = defaultdict(list)

   for file_name in os.listdir(input_dir):
       if file_name.endswith(".faa"):
           file_path = os.path.join(input_dir, file_name)

           for record in SeqIO.parse(file_path, "fasta"):
               match = re.search(r'(\d+)at33208', record.id)
               if match:
                   family_id = match.group(1)
                   gene_family_dict[family_id].append(record)
               else:
                   print(f"Skipping sequence with unexpected format: {record.id}")

   for family_id, sequences in gene_family_dict.items():
       output_file = os.path.join(output_dir, f"{family_id}at33208.faa")
       with open(output_file, "w") as f:
           SeqIO.write(sequences, f, "fasta")

   print("FASTA files have been grouped by gene family and saved in the output directory.")

def main(base_dir, intermediate_dir, final_output_dir, keep_intermediates=False):
   """
   Main function to process .faa files: rename, merge, and group by gene family.
   """
   # Step 1: Rename sequences
   rename_sequences(base_dir)

   # Step 2: Merge sequences by species
   merge_sequences(base_dir, intermediate_dir)

   # Step 3: Group sequences by gene family
   group_by_gene_family(intermediate_dir, final_output_dir)

   # Clean up intermediate files if not needed
   if not keep_intermediates:
       for file in os.listdir(intermediate_dir):
           os.remove(os.path.join(intermediate_dir, file))
       print("Intermediate files have been removed.")

# Paths and settings
BASE_DIR = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025_copy"
INTERMEDIATE_DIR = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025_copy/intermediate"
FINAL_OUTPUT_DIR = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025_copy/final"
KEEP_INTERMEDIATES = True  # Change to True if intermediate files should be kept

# Run the pipeline
main(BASE_DIR, INTERMEDIATE_DIR, FINAL_OUTPUT_DIR, KEEP_INTERMEDIATES)

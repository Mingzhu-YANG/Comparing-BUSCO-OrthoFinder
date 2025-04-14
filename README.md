# Comparing common strategies for ortholog selection used in phylogenomics

These codes are for constructing core BUSCOs on metazoan phylogeny.

Basically, all the scripts would work on the output of the OrthoFinder directories, by doing MAFFT-Linsi sequences alignment and IQ-TREE2 tree building, we can get the topologies of different orthogroups, then we use a customized script to identify the out-paralog OGs and in-paralog sequences(here we only keep one copy to make them single copy, notice we only keep the one that busco have kept), and filtering both for further analysis.


## S1-Data collection

12 chromosome level genome assemblies (10 metazoans and 2 outgroups, and these species are Amphimedon_queenslandica, Mnemiopsis_leydi, Spizellomyces_punctatus, Capitella_teleta, Ciona_intestinalis, Aplysia_californica, Acropora_millepora, Hoilungia_hongkongensis, Priapulus_caudatus, Tribolium_castaneum, Salpingoeca_dolichothecata, and Strongylocentrotus_purpuratus) to perform orthologs inference using BUSCO and OrthoFinder with default parameters. 

![12 chromosome level genome assemblies](/final_combined_plot_reversed_fixed2.png)
this screenshot is showing the taxa we selected and their completeness

BUSCO: Metazoan 954 families

OrthoFinder: Orthogroup_Sequences 71249 OGs, Single_Copy_Orthologur_Sequences

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](/busco_of_dir_output.png)
this is the screenshot of original BUSCO and OrthoFinder output

------------------------------------------------------------------------------------------------------------------------------------------------------------------


## S2-Make comparisons
### S2.0- Get all the busco families using get_sequences.py (remember to create a copy before change)

Go to the output path of BUSCO, it should looks like:

```
qw23953@it103273:~/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025$ ls
Acropora_millepora.fasta.out        Mnemiopsis_leydi.fasta.out
Amphimedon_queenslandica.fasta.out  Priapulus_caudatus.fasta.out
Aplysia_californica.fasta.out       Salpingoeca_dolichothecata.fasta.out
busco_downloads                     Spizellomyces_punctatus.fasta.out
Capitella_teleta.fasta.out          Strongylocentrotus_purpuratus.fasta.out
Ciona_intestinalis.fasta.out        Tribolium_castaneum.fasta.out
Hoilungia_hongkongensis.fasta.out
```

***get_sequences.py***

```
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
```
By using the script above, you will get two new folders in your working directory, they are :

   -/path/to/intermediate
   
   -/path/to/final

You will need the data from /final folder, these are all the busco single-copy gene families, contains all the species sequences belong to that single-copy genes.
later on you would use this dataset to run comparisons.
   
------------------------------------------------------------------------------------------------------------------------------------------------------------------



### S2.1-Copy the Orthogroup_Sequences directory to the working directory (now we work on OrthoFinder output).

1. Using the bash script to group all OGs into different categories, for the Orthogroup_Sequences dir, we have 71249 OGs, we then group them into 2:
   
   - folder OF_1_3: 57459 OGs, using script **get_groups_1_3_seqs.sh**
     these are OGs with less than 3 sequences, cannot be used to build trees, so we have to take them out to folder OF_1_3_seqs.
   - folder OF_4_N: 13790 OGs, we separate them using script **find_single_copy.sh**
       - we further select the OGs with 4-12 sequences without duplication, there are 1659 OGs
       - so the rest 12131 OGs with duplications
       - we can also use **get_strict_single_copy.sh** to get single-copy orthogroups with all taxa present

**get_groups_1_3_seqs.sh**
```
#!/bin/bash

# Define directory
SOURCE_DIR="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451"   
TARGET_DIR="$SOURCE_DIR/OF_1_3_seqs"

# Create target dir (if not exist)
mkdir -p "$TARGET_DIR"

# Loop through all the fa file
for file in "$SOURCE_DIR"/*.fa; do
    # Caculate the number of seqs in each file
    count=$(grep -c "^>" "$file")

    # Check number of seqs, whether they are in between the defined number
    if [ "$count" -ge 1 ] && [ "$count" -le 3 ]; then
        # Move file to target folder
        mv "$file" "$TARGET_DIR/"
        # echo "移动文件: $file"
    fi
done

echo "Done！"
```
**find_single_copy.sh**
```
#!/bin/bash

input_directory="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/OF_4_N"

output_directory="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/OF_4_N/OF_4_12_single"

mkdir -p "$output_directory"

# Loop through all the fasta file
for fasta_file in "$input_directory"/*.fa; do
    # Get the file name
    filename=$(basename "$fasta_file")
    # Check if there are duplicated species ID
    duplicate=$(grep -o '^>[A-Za-z]\{4\}' "$fasta_file" | sort | uniq -d)
    # If no duplicated ID, then move file to target folder
    if [ -z "$duplicate" ]; then
        mv "$fasta_file" "$output_directory/$filename"
        echo "Moved $filename to $output_directory"
    fi
done
```
**get_strict_single_copy.sh**

```
#!/bin/bash

# Define directory
SOURCE_DIR="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/OF_4_N/OF_4_12_single"   
TARGET_DIR="$SOURCE_DIR/OF_12_seqs"

# Create target dir (if not exist)
mkdir -p "$TARGET_DIR"

# Loop through all the fa file
for file in "$SOURCE_DIR"/*.fa; do
    # Caculate the number of seqs in each file
    count=$(grep -c "^>" "$file")

    # Check number of seqs, whether they are in between the defined number
    if [ "$count" -ge 12 ]; then
        # Move file to target folder
        mv "$file" "$TARGET_DIR/"
        echo "Moved file: $file"
    fi
done
```




2. MAFFT alignment

```
nohup bash -c 'for i in *.fa; do mafft --localpair --maxiterate 1000 --preservecase "$i" > "$i.mafft"; done' > mafft.log 2>&1 &
```

       
3. IQtree tree building
   
    - 196 OGs are failed from MAFFT or IQtree, so copy all the successful .fa and .treefile to a new folder called 2_all_OG_multi_13594 (or 11935)

```
for i in *.mafft
do

iqtree2 -s $i -m MFP -mset LG+F+G -madd LG+F+G+C60,LG+F+G+C50,LG+C40+F+G,LG+C30+F+G,LG+C20+F+G,LG+C10+F+G,C10,C20,C30,C40,C50,C60 --score-diff all -mwopt  -bb 1000

done
```
 Or if there are too many files or too big files, you need to configure the parallel tasks to avoid disrupt the system
```
find . -maxdepth 1 -iname "*mafft" > job_array.txt
```
```
nohup bash -c '
    nt=4
    max_parallel_jobs=25
    cat job_array.txt | parallel -j $max_parallel_jobs --no-notice "
        echo \"Running {}\"
        iqtree2 -s {} \
            -m MFP \
            -mset LG+F+G \
            -madd LG+F+G+C60,LG+F+G+C50,LG+C40+F+G,LG+C30+F+G,LG+C20+F+G,LG+C10+F+G,C10,C20,C30,C40,C50,C60 \
            --score-diff all \
            -mwopt \
            -nt ${nt} \
            -bb 1000 > {}.log 2>&1
    "
' > parallel.log 2>&1 &
```

----------------------------------------------------------------------------------------------------------------------------------------------------------------

 
### S2.2-topology analysis (identify paralogs, i.e. in-paralogs and out-paralogs)
1. identify_paralogs_2.py
```
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jun  4 10:54:20 2024
@author: jesuslozanofernandez
Modified on July by Mingzhu Yang
"""

from ete3 import Tree
from collections import defaultdict
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python identify_paralogs.py <treefile>")
    sys.exit(1)

my_tree = sys.argv[1]
print(f"Input file: {my_tree}")

# Get file name and path
base_name = os.path.basename(my_tree)
file_name, file_ext = os.path.splitext(base_name)

# Load trees
tree = Tree(my_tree, format=1)
print("Tree successfully loaded")

my_dict = defaultdict(list)

# Build dictionaries of leafs
for leaf in tree.get_leaf_names():
    my_dict[leaf.split('_')[0]].append(leaf)

inparalogs = []
outparalogs = []

# Identify paralogs
for k, v in my_dict.items():
    if len(v) > 1:
        paralogs = tuple(v)
        is_monophyletic = tree.check_monophyly(paralogs, "name", unrooted=True)
        if is_monophyletic[0]:
            inparalogs.append(str(paralogs) + ': In-Paralogs')
        else:
            outparalogs.append(str(paralogs) + ': Out-Paralogs')


inparalogs_file = f"{file_name}_inparalogs.txt"
outparalogs_file = f"{file_name}_outparalogs.txt"

with open(inparalogs_file, "w") as f:
    for line in inparalogs:
        f.write(line + "\n")

with open(outparalogs_file, "w") as f:
    for line in outparalogs:
        f.write(line + "\n")

print(f"In-Paralogs written to {inparalogs_file}")
print(f"Out-Paralogs written to {outparalogs_file}")

# Check if files are empty and delete if so
for file in [inparalogs_file, outparalogs_file]:
    if os.path.getsize(file) == 0:
        os.remove(file)
        print(f"Deleted empty file: {file}")
```
 the output should like below:  
```
Input file: /user/work/qw23953/6_Compare_Software/8_Redo/OG0000000.fa.mafft.treefile
Tree successfully loaded
In-Paralogs written to OG0000000.fa.mafft_inparalogs.txt
Out-Paralogs written to OG0000000.fa.mafft_outparalogs.txt
Input file: /user/work/qw23953/6_Compare_Software/8_Redo/OG0000001.fa.mafft.treefile
Tree successfully loaded
In-Paralogs written to OG0000001.fa.mafft_inparalogs.txt
Out-Paralogs written to OG0000001.fa.mafft_outparalogs.txt
Input file: /user/work/qw23953/6_Compare_Software/8_Redo/OG0000002.fa.mafft.treefile
Tree successfully loaded
In-Paralogs written to OG0000002.fa.mafft_inparalogs.txt
Out-Paralogs written to OG0000002.fa.mafft_outparalogs.txt
```
2. Prepare the file of outparalog_family.txt
   
   It should looks like
```
OG0000000_tree_outparalogs.txt
OG0000001_tree_outparalogs.txt
OG0000002_tree_outparalogs.txt
OG0000003_tree_outparalogs.txt
OG0000004_tree_outparalogs.txt
...
```
3. Copy all the file into a new folder called 3_all_OG_no_outparalogs/
   
   Move out-paralog files to target_folder="${folder}outparalogs_file/"
   
   using the script of move_outparalog_file.sh
```
import os
import shutil

# Base directory containing all the files
base_dir = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/copy_OF_4_N_list/of_4_n_8139"
# Target directory to move the matched files into
output_dir = os.path.join(base_dir, "outparalog_family")
# Text file that contains the list of *_outparalogs.txt filenames
list_file = os.path.join(base_dir, "outparalog_family.txt")

# Create the output directory if it doesn't already exist
os.makedirs(output_dir, exist_ok=True)

# Read the list file line by line
with open(list_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.endswith("_outparalogs.txt"):
            # Remove the "_outparalogs.txt" suffix to get the prefix
            prefix = line.replace("_outparalogs.txt", "")
            # Construct the corresponding tree file name
            tree_file = f"{prefix}.txt"
            src_path = os.path.join(base_dir, tree_file)
            dst_path = os.path.join(output_dir, tree_file)
            
            # Move the file if it exists
            if os.path.exists(src_path):
                shutil.move(src_path, dst_path)
                print(f"✅ Moved: {tree_file}")
            else:
                print(f"❌ File not found: {tree_file}")        
```
and move all the *_outparalogs.txt file and *inparalogs.txt file by edit tree_file = f"{prefix}.txt" line of codes, the suffix part.

4. Further filtering the in-paralog sequences
   
   Here we going to filter the in-paralog sequences, we keep only once copy thus make it single copy, we keep the sequence that busco have kept.

   To do this, we first cat all the inparalog file (e.g. OG0000000.fa.mafft_inparalogs.txt) into one file called inparalogs_after_out.txt

   it looks like
   
```
('HHON_TR10104_c1_g1_i1_m_21166', 'HHON_TR22145_c0_g1_i1_m_43391', 'HHON_TR3035_c1_g3_i10_m_6505', 'HHON_TR25907_c2_g1_i1_m_51961', 'HHON_TR14192_c0_g2_i1_m_27580', 'HHON_TR28827_c0_g1_i1_m_61152', 'HHON_TR23333_c1_g1_i1_m_46165', 'HHON_TR2887_c0_g1_i5_m_5117', 'HHON_TR19679_c0_g1_i3_m_37965', 'HHON_TR2382_c1_g1_i1_m_3051', 'HHON_TR25907_c3_g1_i1_m_51963', 'HHON_TR3618_c2_g1_i1_m_8986', 'HHON_TR4576_c0_g1_i1_m_11202', 'HHON_TR10152_c0_g2_i3_m_21661', 'HHON_TR19743_c5_g1_i3_m_38508', 'HHON_TR19743_c5_g1_i4_m_38509', 'HHON_TR19743_c5_g1_i5_m_38510', 'HHON_TR1104_c0_g1_i2_m_1420', 'HHON_TR1104_c0_g1_i3_m_1422', 'HHON_TR27469_c1_g1_i1_m_56219', 'HHON_TR1104_c0_g1_i2_m_1421', 'HHON_TR19727_c9_g1_i1_m_38404', 'HHON_TR6491_c0_g1_i1_m_15298', 'HHON_TR28332_c1_g1_i3_m_59665', 'HHON_TR22310_c3_g2_i2_m_43400', 'HHON_TR22310_c3_g2_i3_m_43402', 'HHON_TR16572_c1_g1_i2_m_32193', 'HHON_TR27469_c3_g1_i1_m_56220', 'HHON_TR27533_c0_g2_i1_m_56300', 'HHON_TR6779_c0_g1_i2_m_15473', 'HHON_TR27533_c0_g2_i1_m_56301', 'HHON_TR10172_c0_g1_i1_m_21863', 'HHON_TR3244_c3_g1_i2_m_8220', 'HHON_TR3244_c3_g1_i1_m_8219', 'HHON_TR13760_c5_g4_i2_m_26183', 'HHON_TR13760_c3_g2_i1_m_26180', 'HHON_TR26683_c5_g9_i1_m_54935', 'HHON_TR18475_c0_g4_i1_m_34692', 'HHON_TR28332_c3_g12_i1_m_59672', 'HHON_TR22399_c0_g1_i2_m_44152', 'HHON_TR22399_c0_g1_i5_m_44156', 'HHON_TR22399_c0_g7_i1_m_44160', 'HHON_TR22399_c0_g3_i1_m_44154', 'HHON_TR28332_c3_g13_i1_m_59673', 'HHON_TR13760_c5_g9_i1_m_26186', 'HHON_TR13760_c5_g5_i1_m_26185', 'HHON_TR18475_c0_g8_i1_m_34694', 'HHON_TR18475_c0_g1_i2_m_34690', 'HHON_TR3244_c1_g1_i1_m_8217', 'HHON_TR20643_c0_g1_i1_m_40065', 'HHON_TR220_c1_g1_i1_m_74', 'HHON_TR26683_c5_g1_i1_m_54924', 'HHON_TR28332_c3_g1_i1_m_59669', 'HHON_TR3244_c3_g3_i1_m_8222', 'HHON_TR4147_c2_g1_i1_m_10062', 'HHON_TR10172_c1_g1_i1_m_21864', 'HHON_TR22417_c0_g2_i1_m_44272', 'HHON_TR8303_c0_g1_i1_m_18199', 'HHON_TR4342_c0_g1_i1_m_10312', 'HHON_TR7287_c5_g1_i1_m_15624', 'HHON_TR21144_c0_g1_i1_m_40335', 'HHON_TR8021_c1_g1_i1_m_17414', 'HHON_TR12538_c2_g1_i1_m_24350', 'HHON_TR13760_c3_g1_i1_m_26179', 'HHON_TR16633_c1_g2_i1_m_32605', 'HHON_TR24783_c2_g1_i1_m_50847', 'HHON_TR27396_c3_g1_i1_m_55937', 'HHON_TR3244_c3_g5_i1_m_8223', 'HHON_TR619_c5_g1_i1_m_296', 'HHON_TR10704_c0_g1_i1_m_22772', 'HHON_TR10172_c3_g2_i3_m_21867', 'HHON_TR10172_c4_g11_i1_m_21878', 'HHON_TR10172_c4_g7_i1_m_21872', 'HHON_TR4410_c1_g1_i1_m_10667', 'HHON_TR4410_c1_g1_i2_m_10668', 'HHON_TR18004_c0_g1_i1_m_34041', 'HHON_TR27396_c4_g1_i1_m_55939', 'HHON_TR10172_c4_g8_i3_m_21874', 'HHON_TR28332_c3_g11_i2_m_59671', 'HHON_TR12538_c3_g1_i2_m_24351', 'HHON_TR28332_c3_g13_i2_m_59674', 'HHON_TR10172_c4_g10_i2_m_21876', 'HHON_TR10172_c4_g10_i3_m_21877', 'HHON_TR28332_c2_g1_i1_m_59668', 'HHON_TR29947_c3_g2_i1_m_64089', 'HHON_TR27396_c6_g12_i1_m_55945', 'HHON_TR13570_c0_g1_i1_m_25812', 'HHON_TR13760_c5_g3_i1_m_26181', 'HHON_TR14167_c4_g3_i1_m_27282', 'HHON_TR14167_c4_g5_i1_m_27283', 'HHON_TR27707_c0_g3_i1_m_56528', 'HHON_TR14167_c4_g8_i1_m_27285', 'HHON_TR18475_c0_g3_i1_m_34691', 'HHON_TR1901_c0_g1_i1_m_1660', 'HHON_TR27447_c13_g7_i1_m_56136', 'HHON_TR986_c2_g2_i1_m_1060', 'HHON_TR986_c2_g9_i1_m_1061', 'HHON_TR28332_c3_g6_i1_m_59670', 'HHON_TR5648_c0_g1_i1_m_13193', 'HHON_TR18642_c0_g1_i1_m_35192', 'HHON_TR3471_c0_g2_i1_m_8722', 'HHON_TR16668_c10_g2_i3_m_32801', 'HHON_TR16668_c10_g2_i5_m_32803', 'HHON_TR31066_c0_g3_i1_m_67109', 'HHON_TR31066_c0_g7_i1_m_67110', 'HHON_TR16767_c1_g2_i2_m_32983', 'HHON_TR14219_c2_g1_i1_m_27735', 'HHON_TR10704_c1_g1_i1_m_22773', 'HHON_TR10704_c1_g1_i2_m_22774', 'HHON_TR19202_c1_g1_i1_m_36190', 'HHON_TR7287_c3_g1_i1_m_15614', 'HHON_TR16450_c1_g1_i1_m_31334', 'HHON_TR19202_c0_g1_i1_m_36189', 'HHON_TR30971_c0_g1_i1_m_66647', 'HHON_TR10172_c4_g4_i1_m_21870', 'HHON_TR619_c4_g3_i5_m_291', 'HHON_TR619_c4_g3_i3_m_288', 'HHON_TR619_c4_g3_i1_m_286', 'HHON_TR26637_c2_g1_i1_m_54640', 'HHON_TR619_c4_g3_i8_m_294', 'HHON_TR619_c4_g3_i2_m_287', 'HHON_TR25941_c0_g1_i1_m_51975', 'HHON_TR3244_c2_g1_i1_m_8218', 'HHON_TR16697_c0_g1_i1_m_32955', 'HHON_TR10361_c0_g1_i1_m_22729', 'HHON_TR7811_c0_g1_i1_m_16321', 'HHON_TR10314_c0_g1_i1_m_22722', 'HHON_TR22995_c1_g1_i1_m_45511', 'HHON_TR22995_c1_g1_i2_m_45513', 'HHON_TR24861_c0_g1_i1_m_51126', 'HHON_TR12543_c0_g10_i1_m_24363', 'HHON_TR12543_c0_g9_i1_m_24362', 'HHON_TR8126_c0_g1_i2_m_17936', 'HHON_TR11676_c0_g1_i1_m_23765', 'HHON_TR14256_c1_g2_i1_m_27946', 'HHON_TR14256_c0_g1_i1_m_27945', 'HHON_TR17507_c0_g1_i1_m_33268', 'HHON_TR13033_c5_g8_i1_m_25288', 'HHON_TR16428_c0_g1_i1_m_31328', 'HHON_TR26683_c5_g8_i1_m_54933', 'HHON_TR7144_c0_g1_i1_m_15502', 'HHON_TR1616_c0_g1_i1_m_1635', 'HHON_TR25660_c0_g1_i1_m_51860', 'HHON_TR22310_c1_g1_i1_m_43398', 'HHON_TR12538_c3_g2_i2_m_24352', 'HHON_TR12538_c3_g4_i1_m_24355', 'HHON_TR31066_c0_g8_i2_m_67112', 'HHON_TR14167_c3_g1_i1_m_27280', 'HHON_TR14167_c3_g1_i2_m_27281', 'HHON_TR18321_c0_g2_i1_m_34600', 'HHON_TR20491_c0_g1_i1_m_40042', 'HHON_TR22995_c0_g1_i2_m_45510', 'HHON_TR15632_c1_g10_i1_m_29594', 'HHON_TR3052_c3_g10_i2_m_6639', 'HHON_TR3052_c3_g11_i1_m_6638', 'HHON_TR3052_c3_g12_i7_m_6644', 'HHON_TR3052_c3_g12_i2_m_6641', 'HHON_TR3052_c3_g12_i3_m_6642', 'HHON_TR3052_c3_g12_i5_m_6643', 'HHON_TR9414_c0_g1_i1_m_20541', 'HHON_TR2533_c0_g1_i1_m_3609', 'HHON_TR24586_c0_g1_i1_m_49953', 'HHON_TR9591_c1_g1_i2_m_20683', 'HHON_TR13810_c0_g1_i1_m_26330', 'HHON_TR26683_c5_g4_i1_m_54926', 'HHON_TR12801_c1_g2_i1_m_24607', 'HHON_TR14256_c4_g16_i1_m_27954', 'HHON_TR14256_c4_g4_i1_m_27950', 'HHON_TR17551_c8_g1_i1_m_33347', 'HHON_TR4558_c3_g10_i1_m_11073', 'HHON_TR4558_c3_g5_i1_m_11071', 'HHON_TR4558_c3_g11_i2_m_11075', 'HHON_TR4558_c3_g1_i1_m_11068', 'HHON_TR4558_c3_g2_i1_m_11069', 'HHON_TR14256_c4_g8_i3_m_27952', 'HHON_TR19740_c5_g10_i1_m_38478', 'HHON_TR19740_c5_g1_i1_m_38468', 'HHON_TR2477_c0_g2_i1_m_3442', 'HHON_TR4558_c3_g11_i1_m_11074', 'HHON_TR26683_c5_g10_i1_m_54936', 'HHON_TR26683_c5_g5_i3_m_54929', 'HHON_TR4558_c2_g1_i2_m_11067', 'HHON_TR26683_c5_g5_i4_m_54932', 'HHON_TR4790_c0_g1_i1_m_12609', 'HHON_TR4790_c0_g2_i1_m_12610', 'HHON_TR26683_c5_g7_i1_m_54928', 'HHON_TR21144_c1_g1_i1_m_40337', 'HHON_TR16767_c0_g1_i1_m_32981', 'HHON_TR19090_c0_g1_i1_m_36180', 'HHON_TR16767_c1_g1_i1_m_32982', 'HHON_TR4795_c0_g2_i1_m_12611', 'HHON_TR20630_c0_g2_i1_m_40064', 'HHON_TR22995_c1_g2_i1_m_45512', 'HHON_TR26683_c5_g7_i2_m_54930', 'HHON_TR20214_c0_g1_i1_m_39316', 'HHON_TR3555_c0_g1_i1_m_8781', 'HHON_TR20764_c0_g1_i1_m_40086', 'HHON_TR220_c3_g1_i1_m_76', 'HHON_TR4558_c3_g8_i1_m_11072', 'HHON_TR220_c2_g1_i1_m_75', 'HHON_TR19720_c1_g2_i1_m_38310', 'HHON_TR4147_c5_g1_i3_m_10067', 'HHON_TR26637_c2_g3_i1_m_54642', 'HHON_TR26637_c2_g3_i2_m_54643', 'HHON_TR26637_c1_g1_i2_m_54638', 'HHON_TR26637_c1_g2_i1_m_54637', 'HHON_TR30753_c0_g1_i1_m_66625', 'HHON_TR12538_c3_g7_i3_m_24357', 'HHON_TR14256_c2_g1_i1_m_27947', 'HHON_TR14256_c2_g1_i2_m_27948', 'HHON_TR4558_c4_g1_i1_m_11077', 'HHON_TR23840_c3_g10_i2_m_47424', 'HHON_TR23840_c3_g4_i1_m_47423', 'HHON_TR26052_c0_g2_i1_m_52637', 'HHON_TR26052_c0_g4_i1_m_52638', 'HHON_TR26052_c0_g4_i2_m_52640', 'HHON_TR26052_c0_g6_i6_m_52646', 'HHON_TR26052_c0_g6_i9_m_52647', 'HHON_TR17143_c0_g1_i1_m_33028', 'HHON_TR2596_c0_g1_i2_m_3614', 'HHON_TR2596_c0_g1_i3_m_3615', 'HHON_TR27396_c6_g10_i1_m_55943', 'HHON_TR27396_c6_g9_i1_m_55941', 'HHON_TR4573_c0_g1_i2_m_11185', 'HHON_TR27396_c6_g14_i1_m_55946', 'HHON_TR27396_c6_g11_i1_m_55944', 'HHON_TR27396_c6_g4_i1_m_55940', 'HHON_TR31423_c2_g6_i1_m_68294', 'HHON_TR31423_c2_g6_i2_m_68295', 'HHON_TR31423_c2_g7_i1_m_68296', 'HHON_TR6014_c1_g7_i3_m_14720', 'HHON_TR6014_c1_g7_i6_m_14722', 'HHON_TR13005_c0_g1_i1_m_25044', 'HHON_TR29878_c0_g1_i1_m_63784', 'HHON_TR27762_c2_g1_i6_m_56915', 'HHON_TR29947_c6_g1_i2_m_64091', 'HHON_TR25835_c0_g1_i1_m_51924', 'HHON_TR25835_c0_g1_i2_m_51925', 'HHON_TR27475_c2_g1_i2_m_56241', 'HHON_TR27475_c2_g1_i3_m_56242', 'HHON_TR27475_c2_g1_i6_m_56245', 'HHON_TR29902_c3_g2_i4_m_63916', 'HHON_TR29902_c3_g8_i1_m_63927', 'HHON_TR10593_c0_g1_i1_m_22748', 'HHON_TR10593_c0_g2_i1_m_22749', 'HHON_TR3870_c0_g1_i1_m_9826', 'HHON_TR13411_c1_g1_i1_m_25777', 'HHON_TR26989_c6_g4_i1_m_55405', 'HHON_TR9848_c0_g1_i1_m_20758', 'HHON_TR11291_c0_g1_i1_m_22885', 'HHON_TR12044_c1_g1_i1_m_24055', 'HHON_TR12340_c2_g1_i1_m_24217', 'HHON_TR31268_c0_g2_i1_m_67716', 'HHON_TR13033_c4_g2_i2_m_25283', 'HHON_TR19728_c0_g1_i1_m_38410', 'HHON_TR12334_c1_g1_i1_m_24210', 'HHON_TR22506_c0_g1_i1_m_44282', 'HHON_TR26685_c1_g4_i1_m_54950', 'HHON_TR3719_c3_g1_i1_m_9670', 'HHON_TR3719_c4_g6_i1_m_9681', 'HHON_TR13033_c4_g1_i2_m_25281', 'HHON_TR18649_c0_g1_i1_m_35194', 'HHON_TR238_c1_g1_i2_m_102', 'HHON_TR13910_c0_g1_i1_m_26341', 'HHON_TR19728_c0_g1_i3_m_38411', 'HHON_TR13033_c3_g1_i1_m_25279', 'HHON_TR28485_c0_g1_i1_m_60266', 'HHON_TR23320_c4_g1_i3_m_46050', 'HHON_TR20925_c1_g1_i1_m_40155', 'HHON_TR9591_c0_g1_i1_m_20682', 'HHON_TR13033_c5_g7_i1_m_25286', 'HHON_TR19728_c0_g4_i2_m_38413', 'HHON_TR13033_c2_g1_i1_m_25278', 'HHON_TR13033_c5_g2_i1_m_25284', 'HHON_TR17551_c5_g1_i1_m_33346', 'HHON_TR2737_c0_g1_i1_m_3795', 'HHON_TR26989_c6_g8_i1_m_55412', 'HHON_TR14167_c1_g2_i1_m_27278', 'HHON_TR5229_c4_g2_i1_m_13067', 'HHON_TR16146_c0_g1_i1_m_31196', 'HHON_TR18578_c1_g1_i1_m_35185', 'HHON_TR14415_c0_g1_i1_m_28340', 'HHON_TR16668_c10_g1_i1_m_32798', 'HHON_TR16668_c9_g1_i1_m_32797', 'HHON_TR26989_c1_g1_i1_m_55398', 'HHON_TR20925_c0_g1_i1_m_40154', 'HHON_TR3244_c3_g2_i1_m_8221', 'HHON_TR12801_c0_g1_i1_m_24605', 'HHON_TR12801_c0_g2_i1_m_24606', 'HHON_TR23320_c4_g4_i1_m_46054', 'HHON_TR23320_c4_g5_i3_m_46058', 'HHON_TR26683_c1_g2_i1_m_54922', 'HHON_TR26989_c6_g5_i1_m_55406', 'HHON_TR26683_c4_g1_i1_m_54923', 'HHON_TR22073_c0_g1_i1_m_43372', 'HHON_TR24026_c1_g1_i2_m_48376', 'HHON_TR24386_c0_g2_i1_m_49911', 'HHON_TR27996_c0_g1_i1_m_57732', 'HHON_TR26989_c6_g7_i1_m_55411', 'HHON_TR25787_c2_g1_i1_m_51883', 'HHON_TR27447_c13_g3_i1_m_56133', 'HHON_TR13005_c1_g11_i1_m_25049', 'HHON_TR4578_c4_g5_i1_m_11222', 'HHON_TR29059_c0_g11_i1_m_61422', 'HHON_TR4578_c4_g3_i1_m_11221', 'HHON_TR13005_c1_g1_i1_m_25045', 'HHON_TR13005_c1_g2_i2_m_25046', 'HHON_TR3621_c1_g11_i1_m_9006', 'HHON_TR4578_c4_g1_i5_m_11219', 'HHON_TR4578_c4_g9_i1_m_11223', 'HHON_TR4460_c0_g1_i1_m_10944', 'HHON_TR18475_c0_g7_i1_m_34693', 'HHON_TR14039_c1_g1_i1_m_26416', 'HHON_TR12538_c3_g2_i3_m_24353', 'HHON_TR14283_c0_g1_i2_m_28166', 'HHON_TR4147_c4_g1_i1_m_10063', 'HHON_TR4147_c4_g1_i2_m_10065', 'HHON_TR4147_c4_g2_i2_m_10066', 'HHON_TR31423_c2_g8_i1_m_68297', 'HHON_TR10209_c0_g1_i1_m_22219', 'HHON_TR4147_c5_g2_i1_m_10068', 'HHON_TR10209_c0_g1_i3_m_22221', 'HHON_TR26637_c2_g4_i1_m_54644', 'HHON_TR4546_c3_g1_i1_m_11003', 'HHON_TR15429_c0_g1_i1_m_29392', 'HHON_TR13548_c0_g1_i1_m_25806', 'HHON_TR22947_c0_g1_i1_m_45422', 'HHON_TR22947_c0_g1_i6_m_45428', 'HHON_TR15632_c0_g1_i1_m_29592', 'HHON_TR2768_c1_g2_i3_m_4130', 'HHON_TR2768_c1_g2_i4_m_4131', 'HHON_TR27447_c13_g2_i1_m_56131', 'HHON_TR27269_c0_g2_i4_m_55900', 'HHON_TR27269_c0_g2_i7_m_55903', 'HHON_TR19742_c2_g3_i3_m_38491', 'HHON_TR19742_c2_g4_i1_m_38490', 'HHON_TR4558_c3_g14_i1_m_11076', 'HHON_TR3486_c0_g1_i2_m_8757', 'HHON_TR19742_c2_g3_i4_m_38492', 'HHON_TR19742_c2_g3_i5_m_38493', 'HHON_TR27719_c0_g1_i1_m_56587', 'HHON_TR27719_c0_g1_i5_m_56590', 'HHON_TR27719_c0_g1_i9_m_56594', 'HHON_TR29338_c0_g1_i2_m_63158', 'HHON_TR19742_c2_g3_i7_m_38495', 'HHON_TR3024_c0_g3_i3_m_6308', 'HHON_TR3024_c0_g3_i2_m_6307', 'HHON_TR30610_c1_g1_i1_m_66560', 'HHON_TR24229_c1_g2_i1_m_49278', 'HHON_TR30650_c2_g1_i1_m_66588', 'HHON_TR3024_c0_g2_i1_m_6303', 'HHON_TR30650_c1_g1_i1_m_66587', 'HHON_TR19742_c3_g2_i1_m_38500', 'HHON_TR24229_c6_g1_i1_m_49279', 'HHON_TR24229_c6_g7_i1_m_49281', 'HHON_TR25907_c1_g1_i1_m_51960', 'HHON_TR24229_c6_g2_i1_m_49280', 'HHON_TR20205_c0_g1_i1_m_39310', 'HHON_TR20205_c0_g2_i1_m_39311', 'HHON_TR21193_c2_g3_i1_m_40847', 'HHON_TR29386_c0_g1_i1_m_63537', 'HHON_TR29386_c0_g2_i4_m_63540', 'HHON_TR19742_c1_g1_i1_m_38488', 'HHON_TR19742_c3_g1_i4_m_38503', 'HHON_TR2975_c3_g1_i1_m_5922', 'HHON_TR2975_c3_g1_i2_m_5923', 'HHON_TR2975_c3_g1_i7_m_5928', 'HHON_TR2975_c3_g1_i8_m_5930', 'HHON_TR2975_c3_g2_i1_m_5929', 'HHON_TR19742_c3_g1_i2_m_38499', 'HHON_TR3502_c0_g1_i1_m_8766', 'HHON_TR2975_c3_g3_i1_m_5932', 'HHON_TR16714_c1_g1_i1_m_32964', 'HHON_TR19699_c2_g12_i1_m_38135', 'HHON_TR19699_c2_g6_i2_m_38117', 'HHON_TR19699_c2_g6_i4_m_38119', 'HHON_TR19699_c2_g6_i6_m_38121', 'HHON_TR19699_c2_g6_i8_m_38122', 'HHON_TR19699_c2_g6_i1_m_38116', 'HHON_TR19699_c2_g9_i6_m_38131', 'HHON_TR19699_c2_g9_i7_m_38132', 'HHON_TR19699_c2_g9_i4_m_38129', 'HHON_TR19699_c2_g9_i5_m_38130', 'HHON_TR8456_c0_g2_i1_m_18235', 'HHON_TR19699_c2_g9_i9_m_38134', 'HHON_TR19720_c1_g4_i1_m_38312', 'HHON_TR19720_c1_g8_i1_m_38318', 'HHON_TR19699_c2_g1_i1_m_38112', 'HHON_TR19720_c1_g6_i2_m_38315', 'HHON_TR19720_c1_g6_i3_m_38316', 'HHON_TR2675_c0_g1_i1_m_3650', 'HHON_TR9874_c0_g1_i1_m_20761', 'HHON_TR2675_c1_g1_i1_m_3652', 'HHON_TR2675_c1_g2_i1_m_3653', 'HHON_TR28574_c2_g3_i1_m_60809', 'HHON_TR28574_c2_g3_i5_m_60813', 'HHON_TR28574_c2_g4_i1_m_60810', 'HHON_TR23356_c1_g1_i1_m_46358', 'HHON_TR19699_c2_g8_i1_m_38124', 'HHON_TR19720_c2_g1_i1_m_38325', 'HHON_TR25907_c0_g1_i1_m_51959', 'HHON_TR19705_c0_g1_i1_m_38171', 'HHON_TR19705_c0_g1_i6_m_38177', 'HHON_TR19705_c0_g1_i8_m_38181', 'HHON_TR15022_c0_g1_i1_m_28468', 'HHON_TR2361_c0_g1_i1_m_2941', 'HHON_TR2361_c1_g1_i2_m_2943', 'HHON_TR28545_c1_g1_i2_m_60573', 'HHON_TR28545_c1_g2_i1_m_60574'): In-Paralogs
('SDOL_gnl_est_initial1030.p2', 'SDOL_gnl_est_initial10959.p1', 'SDOL_gnl_est_initial10960.p1', 'SDOL_gnl_est_initial27637.p1', 'SDOL_gnl_est_initial27627.p1', 'SDOL_gnl_est_initial27630.p1', 'SDOL_gnl_est_initial27628.p1', 'SDOL_gnl_est_initial27629.p1', 'SDOL_gnl_est_initial27632.p1', 'SDOL_gnl_est_initial27636.p1', 'SDOL_gnl_est_initial27634.p1', 'SDOL_gnl_est_initial26538.p1', 'SDOL_gnl_est_initial27318.p1', 'SDOL_gnl_est_initial26539.p1', 'SDOL_gnl_est_initial6728.p1', 'SDOL_gnl_est_initial27635.p1', 'SDOL_gnl_est_initial27631.p1', 'SDOL_gnl_est_initial20542.p1', 'SDOL_gnl_est_initial3175.p1', 'SDOL_gnl_est_initial21653.p1', 'SDOL_gnl_est_initial8063.p1', 'SDOL_gnl_est_initial26540.p2', 'SDOL_gnl_est_initial31500.p1', 'SDOL_gnl_est_initial1185.p1', 'SDOL_gnl_est_initial26209.p2', 'SDOL_gnl_est_initial26214.p1', 'SDOL_gnl_est_initial26269.p1', 'SDOL_gnl_est_initial26271.p1', 'SDOL_gnl_est_initial26270.p1', 'SDOL_gnl_est_initial15120.p1', 'SDOL_gnl_est_initial15121.p1', 'SDOL_gnl_est_initial27900.p1', 'SDOL_gnl_est_initial27901.p1', 'SDOL_gnl_est_initial26973.p1', 'SDOL_gnl_est_initial26975.p1', 'SDOL_gnl_est_initial26974.p1', 'SDOL_gnl_est_initial4726.p4', 'SDOL_gnl_est_initial26657.p1', 'SDOL_gnl_est_initial26972.p1', 'SDOL_gnl_est_initial27899.p1', 'SDOL_gnl_est_initial15236.p1', 'SDOL_gnl_est_initial26272.p1', 'SDOL_gnl_est_initial7413.p2', 'SDOL_gnl_est_initial13417.p1', 'SDOL_gnl_est_initial8146.p2', 'SDOL_gnl_est_initial27804.p1', 'SDOL_gnl_est_initial27806.p1', 'SDOL_gnl_est_initial27807.p2', 'SDOL_gnl_est_initial26207.p1', 'SDOL_gnl_est_initial27803.p1', 'SDOL_gnl_est_initial1365.p1', 'SDOL_gnl_est_initial26917.p1', 'SDOL_gnl_est_initial26919.p1', 'SDOL_gnl_est_initial7401.p2', 'SDOL_gnl_est_initial26206.p1', 'SDOL_gnl_est_initial29045.p1', 'SDOL_gnl_est_initial15466.p1', 'SDOL_gnl_est_initial26210.p1', 'SDOL_gnl_est_initial26212.p1', 'SDOL_gnl_est_initial26660.p1', 'SDOL_gnl_est_initial13415.p1', 'SDOL_gnl_est_initial29940.p1', 'SDOL_gnl_est_initial13416.p1', 'SDOL_gnl_est_initial19008.p1', 'SDOL_gnl_est_initial26653.p2', 'SDOL_gnl_est_initial26655.p1', 'SDOL_gnl_est_initial28587.p1', 'SDOL_gnl_est_initial26213.p1', 'SDOL_gnl_est_initial15235.p1', 'SDOL_gnl_est_initial31037.p1', 'SDOL_gnl_est_initial27322.p1', 'SDOL_gnl_est_initial28774.p1', 'SDOL_gnl_est_initial6776.p2', 'SDOL_gnl_est_initial3134.p1', 'SDOL_gnl_est_initial12938.p1', 'SDOL_gnl_est_initial40383.p3', 'SDOL_gnl_est_initial22481.p1', 'SDOL_gnl_est_initial42803.p2', 'SDOL_gnl_est_initial22482.p1', 'SDOL_gnl_est_initial22483.p3', 'SDOL_gnl_est_initial15077.p1', 'SDOL_gnl_est_initial2958.p1', 'SDOL_gnl_est_initial27344.p1', 'SDOL_gnl_est_initial27345.p1', 'SDOL_gnl_est_initial26498.p1', 'SDOL_gnl_est_initial26747.p1', 'SDOL_gnl_est_initial26748.p1', 'SDOL_gnl_est_initial3555.p1', 'SDOL_gnl_est_initial27769.p1', 'SDOL_gnl_est_initial1299.p1', 'SDOL_gnl_est_initial27846.p1', 'SDOL_gnl_est_initial27847.p1', 'SDOL_gnl_est_initial25908.p1', 'SDOL_gnl_est_initial32274.p2', 'SDOL_gnl_est_initial27849.p2', 'SDOL_gnl_est_initial27044.p1', 'SDOL_gnl_est_initial27853.p1', 'SDOL_gnl_est_initial28672.p2', 'SDOL_gnl_est_initial27898.p1', 'SDOL_gnl_est_initial27902.p1', 'SDOL_gnl_est_initial27045.p1', 'SDOL_gnl_est_initial3650.p1', 'SDOL_gnl_est_initial11182.p1', 'SDOL_gnl_est_initial27767.p1', 'SDOL_gnl_est_initial2742.p1', 'SDOL_gnl_est_initial26104.p1', 'SDOL_gnl_est_initial11180.p1', 'SDOL_gnl_est_initial11181.p1', 'SDOL_gnl_est_initial13551.p1', 'SDOL_gnl_est_initial47987.p1', 'SDOL_gnl_est_initial11179.p1', 'SDOL_gnl_est_initial25728.p1', 'SDOL_gnl_est_initial5603.p1', 'SDOL_gnl_est_initial43659.p1', 'SDOL_gnl_est_initial7128.p1', 'SDOL_gnl_est_initial1119.p1', 'SDOL_gnl_est_initial27315.p1', 'SDOL_gnl_est_initial27316.p2', 'SDOL_gnl_est_initial26035.p1', 'SDOL_gnl_est_initial26037.p2', 'SDOL_gnl_est_initial27343.p2', 'SDOL_gnl_est_initial43745.p1', 'SDOL_gnl_est_initial17668.p1', 'SDOL_gnl_est_initial25729.p1', 'SDOL_gnl_est_initial27340.p1', 'SDOL_gnl_est_initial26656.p1', 'SDOL_gnl_est_initial27046.p1', 'SDOL_gnl_est_initial26405.p1', 'SDOL_gnl_est_initial26408.p1', 'SDOL_gnl_est_initial26406.p1', 'SDOL_gnl_est_initial43837.p1', 'SDOL_gnl_est_initial26407.p1', 'SDOL_gnl_est_initial26409.p1', 'SDOL_gnl_est_initial43835.p1', 'SDOL_gnl_est_initial43838.p1', 'SDOL_gnl_est_initial43839.p2', 'SDOL_gnl_est_initial2442.p2', 'SDOL_gnl_est_initial18804.p1', 'SDOL_gnl_est_initial27042.p1', 'SDOL_gnl_est_initial10961.p1', 'SDOL_gnl_est_initial24986.p1', 'SDOL_gnl_est_initial24987.p1', 'SDOL_gnl_est_initial3857.p1', 'SDOL_gnl_est_initial14064.p1', 'SDOL_gnl_est_initial27893.p2', 'SDOL_gnl_est_initial27894.p1', 'SDOL_gnl_est_initial4108.p2', 'SDOL_gnl_est_initial24371.p1', 'SDOL_gnl_est_initial26087.p1', 'SDOL_gnl_est_initial27041.p2', 'SDOL_gnl_est_initial27341.p1', 'SDOL_gnl_est_initial27342.p1', 'SDOL_gnl_est_initial17669.p1', 'SDOL_gnl_est_initial27043.p1', 'SDOL_gnl_est_initial8225.p2', 'SDOL_gnl_est_initial4818.p1', 'SDOL_gnl_est_initial26086.p2'): In-Paralogs
('HHON_TR10068_c4_g1_i2_m_20869', 'HHON_TR10068_c5_g1_i1_m_20870', 'HHON_TR13161_c0_g1_i1_m_25695', 'HHON_TR10068_c5_g4_i1_m_20871', 'HHON_TR23729_c0_g1_i1_m_46536', 'HHON_TR10104_c0_g1_i1_m_21165', 'HHON_TR29135_c0_g1_i1_m_62144', 'HHON_TR12340_c1_g3_i1_m_24216', 'HHON_TR21485_c0_g1_i1_m_41512', 'HHON_TR21989_c0_g2_i2_m_42557', 'HHON_TR21989_c0_g2_i7_m_42562', 'HHON_TR2444_c0_g1_i1_m_3331', 'HHON_TR27898_c0_g1_i1_m_57272', 'HHON_TR14089_c1_g3_i2_m_26533', 'HHON_TR27939_c0_g1_i1_m_57536', 'HHON_TR14192_c0_g1_i1_m_27579', 'HHON_TR18230_c0_g1_i1_m_34572', 'HHON_TR21697_c1_g1_i1_m_41800', 'HHON_TR22001_c0_g1_i1_m_42667', 'HHON_TR21998_c0_g1_i1_m_42628', 'HHON_TR21998_c0_g1_i3_m_42630', 'HHON_TR28284_c0_g1_i1_m_59111', 'HHON_TR22645_c0_g1_i1_m_44641', 'HHON_TR28543_c0_g1_i1_m_60560', 'HHON_TR23985_c1_g2_i1_m_48344', 'HHON_TR14725_c0_g1_i1_m_28375', 'HHON_TR19725_c3_g2_i4_m_38402', 'HHON_TR8663_c0_g1_i1_m_18754', 'HHON_TR8663_c0_g1_i6_m_18759', 'HHON_TR25822_c0_g1_i1_m_51917', 'HHON_TR31268_c0_g4_i1_m_67721', 'HHON_TR31268_c0_g2_i4_m_67722', 'HHON_TR31268_c0_g3_i1_m_67717', 'HHON_TR31268_c0_g3_i3_m_67723', 'HHON_TR12809_c0_g1_i1_m_24614', 'HHON_TR19725_c2_g1_i1_m_38388', 'HHON_TR2675_c3_g1_i1_m_3654', 'HHON_TR18539_c0_g1_i2_m_35016', 'HHON_TR27215_c0_g1_i1_m_55679', 'HHON_TR25813_c1_g1_i1_m_51909', 'HHON_TR22811_c0_g1_i1_m_44886', 'HHON_TR9809_c1_g1_i1_m_20752', 'HHON_TR11087_c0_g2_i1_m_22833', 'HHON_TR27469_c7_g1_i1_m_56236', 'HHON_TR21989_c0_g1_i1_m_42552', 'HHON_TR31371_c0_g1_i2_m_68054', 'HHON_TR23320_c3_g1_i1_m_46044', 'HHON_TR28664_c0_g1_i1_m_61129', 'HHON_TR8380_c2_g1_i1_m_18218', 'HHON_TR31515_c0_g1_i1_m_68352', 'HHON_TR13033_c7_g1_i1_m_25289', 'HHON_TR14167_c0_g4_i1_m_27277', 'HHON_TR20824_c1_g1_i1_m_40092', 'HHON_TR13033_c4_g2_i1_m_25282', 'HHON_TR22602_c0_g1_i1_m_44377', 'HHON_TR22849_c3_g1_i1_m_45102', 'HHON_TR22849_c3_g1_i3_m_45104', 'HHON_TR30411_c0_g5_i1_m_65655', 'HHON_TR10068_c6_g1_i1_m_20874', 'HHON_TR10641_c0_g1_i1_m_22757', 'HHON_TR3719_c0_g1_i1_m_9666', 'HHON_TR3719_c4_g4_i1_m_9675', 'HHON_TR238_c3_g1_i1_m_108', 'HHON_TR238_c3_g1_i2_m_109', 'HHON_TR28664_c1_g1_i1_m_61130', 'HHON_TR3719_c4_g4_i1_m_9676', 'HHON_TR8380_c1_g1_i1_m_18217', 'HHON_TR238_c2_g1_i1_m_104', 'HHON_TR238_c2_g2_i1_m_106', 'HHON_TR8349_c0_g1_i1_m_18212', 'HHON_TR3719_c1_g1_i1_m_9667', 'HHON_TR3719_c4_g3_i1_m_9674', 'HHON_TR3719_c2_g1_i1_m_9669', 'HHON_TR3719_c4_g4_i3_m_9678', 'HHON_TR3719_c4_g1_i1_m_9671', 'HHON_TR3719_c4_g5_i1_m_9679', 'HHON_TR5894_c0_g1_i1_m_13961', 'HHON_TR238_c2_g3_i1_m_107', 'HHON_TR3719_c4_g2_i1_m_9672', 'HHON_TR13033_c4_g1_i1_m_25280', 'HHON_TR16668_c6_g1_i1_m_32794', 'HHON_TR9110_c0_g1_i1_m_19785', 'HHON_TR11500_c3_g4_i1_m_23128', 'HHON_TR27247_c0_g1_i3_m_55824', 'HHON_TR18711_c0_g1_i1_m_35207', 'HHON_TR26039_c0_g1_i1_m_52563', 'HHON_TR16668_c8_g1_i1_m_32795', 'HHON_TR9681_c0_g2_i1_m_20738', 'HHON_TR6599_c0_g1_i1_m_15424', 'HHON_TR14808_c0_g1_i1_m_28384', 'HHON_TR26989_c6_g6_i1_m_55408', 'HHON_TR26039_c5_g1_i1_m_52564', 'HHON_TR7295_c7_g2_i4_m_15668', 'HHON_TR17551_c10_g1_i1_m_33348', 'HHON_TR26989_c6_g9_i1_m_55413', 'HHON_TR721_c2_g1_i1_m_697', 'HHON_TR26989_c6_g3_i1_m_55403', 'HHON_TR5229_c5_g1_i1_m_13069', 'HHON_TR7295_c4_g3_i1_m_15656', 'HHON_TR8816_c0_g1_i1_m_19014', 'HHON_TR16067_c0_g1_i1_m_31191', 'HHON_TR5229_c4_g1_i1_m_13065', 'HHON_TR7295_c7_g2_i1_m_15662', 'HHON_TR7295_c7_g2_i3_m_15665', 'HHON_TR23320_c4_g6_i1_m_46059', 'HHON_TR23320_c4_g6_i2_m_46060', 'HHON_TR16712_c0_g1_i1_m_32963', 'HHON_TR23320_c0_g1_i1_m_46042', 'HHON_TR23320_c4_g7_i1_m_46063', 'HHON_TR26989_c4_g1_i1_m_55399', 'HHON_TR17282_c0_g2_i3_m_33080', 'HHON_TR27240_c0_g6_i1_m_55792', 'HHON_TR28918_c0_g1_i1_m_61166', 'HHON_TR7605_c0_g1_i1_m_16303', 'HHON_TR20692_c0_g1_i1_m_40074', 'HHON_TR29855_c0_g1_i1_m_63714', 'HHON_TR28820_c0_g1_i1_m_61150', 'HHON_TR13462_c0_g1_i1_m_25790', 'HHON_TR5665_c0_g1_i1_m_13195', 'HHON_TR146_c2_g1_i1_m_18', 'HHON_TR7008_c0_g1_i1_m_15490', 'HHON_TR20371_c3_g2_i6_m_39379', 'HHON_TR25942_c0_g1_i1_m_51976', 'HHON_TR19182_c0_g1_i1_m_36186', 'HHON_TR12859_c0_g1_i1_m_24637', 'HHON_TR20925_c2_g1_i1_m_40156', 'HHON_TR30119_c1_g1_i1_m_64190', 'HHON_TR7295_c0_g1_i1_m_15653', 'HHON_TR14648_c0_g1_i1_m_28367', 'HHON_TR21746_c2_g1_i10_m_42251', 'HHON_TR21746_c2_g1_i16_m_42253', 'HHON_TR16720_c0_g1_i1_m_32966', 'HHON_TR5229_c2_g1_i3_m_13061', 'HHON_TR28185_c0_g1_i1_m_58228', 'HHON_TR25911_c0_g2_i1_m_51964', 'HHON_TR2643_c0_g1_i1_m_3633', 'HHON_TR2643_c0_g2_i1_m_3634', 'HHON_TR17097_c0_g1_i1_m_33024', 'HHON_TR27680_c0_g2_i1_m_56404', 'HHON_TR2737_c1_g1_i1_m_3797', 'HHON_TR14274_c0_g1_i1_m_28121', 'HHON_TR5229_c5_g5_i2_m_13070', 'HHON_TR7295_c4_g1_i1_m_15654', 'HHON_TR10068_c5_g6_i1_m_20873', 'HHON_TR22316_c3_g11_i1_m_43443', 'HHON_TR20438_c0_g2_i3_m_39897', 'HHON_TR27240_c0_g2_i1_m_55790', 'HHON_TR27240_c0_g3_i2_m_55791', 'HHON_TR22316_c3_g2_i2_m_43441', 'HHON_TR29456_c0_g1_i1_m_63671', 'HHON_TR146_c1_g1_i1_m_17'): In-Paralogs
('MLEY_Mnemiopsis_leydi_ML000128a', 'MLEY_Mnemiopsis_leydi_ML003255a', 'MLEY_Mnemiopsis_leydi_ML01806a', 'MLEY_Mnemiopsis_leydi_ML031711a', 'MLEY_Mnemiopsis_leydi_ML049310a', 'MLEY_Mnemiopsis_leydi_ML07027a', 'MLEY_Mnemiopsis_leydi_ML150410a', 'MLEY_Mnemiopsis_leydi_ML078914a', 'MLEY_Mnemiopsis_leydi_ML25741a', 'MLEY_Mnemiopsis_leydi_ML083810a', 'MLEY_Mnemiopsis_leydi_ML21544a', 'MLEY_Mnemiopsis_leydi_ML061520a', 'MLEY_Mnemiopsis_leydi_ML06273a', 'MLEY_Mnemiopsis_leydi_ML04078a', 'MLEY_Mnemiopsis_leydi_ML114621a', 'MLEY_Mnemiopsis_leydi_ML149618a', 'MLEY_Mnemiopsis_leydi_ML17732a', 'MLEY_Mnemiopsis_leydi_ML205714a', 'MLEY_Mnemiopsis_leydi_ML115517a', 'MLEY_Mnemiopsis_leydi_ML115518a', 'MLEY_Mnemiopsis_leydi_ML084410a', 'MLEY_Mnemiopsis_leydi_ML01193a', 'MLEY_Mnemiopsis_leydi_ML06717a', 'MLEY_Mnemiopsis_leydi_ML214311a', 'MLEY_Mnemiopsis_leydi_ML451314a', 'MLEY_Mnemiopsis_leydi_ML221322a', 'MLEY_Mnemiopsis_leydi_ML00679a', 'MLEY_Mnemiopsis_leydi_ML06551a', 'MLEY_Mnemiopsis_leydi_ML28568a', 'MLEY_Mnemiopsis_leydi_ML003265a', 'MLEY_Mnemiopsis_leydi_ML148518a', 'MLEY_Mnemiopsis_leydi_ML16134a', 'MLEY_Mnemiopsis_leydi_ML16594a', 'MLEY_Mnemiopsis_leydi_ML198511a', 'MLEY_Mnemiopsis_leydi_ML07241a', 'MLEY_Mnemiopsis_leydi_ML00391a', 'MLEY_Mnemiopsis_leydi_ML218928a', 'MLEY_Mnemiopsis_leydi_ML128424a', 'MLEY_Mnemiopsis_leydi_ML202618a', 'MLEY_Mnemiopsis_leydi_ML004439a', 'MLEY_Mnemiopsis_leydi_ML32344a', 'MLEY_Mnemiopsis_leydi_ML00643a', 'MLEY_Mnemiopsis_leydi_ML016329a', 'MLEY_Mnemiopsis_leydi_ML160312a', 'MLEY_Mnemiopsis_leydi_ML12424a', 'MLEY_Mnemiopsis_leydi_ML20834a', 'MLEY_Mnemiopsis_leydi_ML20835a', 'MLEY_Mnemiopsis_leydi_ML409410a', 'MLEY_Mnemiopsis_leydi_ML14764a', 'MLEY_Mnemiopsis_leydi_ML006934a', 'MLEY_Mnemiopsis_leydi_ML120737b', 'MLEY_Mnemiopsis_leydi_ML029911a', 'MLEY_Mnemiopsis_leydi_ML029912a', 'MLEY_Mnemiopsis_leydi_ML03596a', 'MLEY_Mnemiopsis_leydi_ML083015a', 'MLEY_Mnemiopsis_leydi_ML069711a', 'MLEY_Mnemiopsis_leydi_ML32093a', 'MLEY_Mnemiopsis_leydi_ML00918a', 'MLEY_Mnemiopsis_leydi_ML038833a', 'MLEY_Mnemiopsis_leydi_ML015730a', 'MLEY_Mnemiopsis_leydi_ML050815a', 'MLEY_Mnemiopsis_leydi_ML06576a', 'MLEY_Mnemiopsis_leydi_ML03194a', 'MLEY_Mnemiopsis_leydi_ML36061a', 'MLEY_Mnemiopsis_leydi_ML36936a', 'MLEY_Mnemiopsis_leydi_ML03238a', 'MLEY_Mnemiopsis_leydi_ML141723a', 'MLEY_Mnemiopsis_leydi_ML141721a', 'MLEY_Mnemiopsis_leydi_ML226712a', 'MLEY_Mnemiopsis_leydi_ML102912a', 'MLEY_Mnemiopsis_leydi_ML160325a', 'MLEY_Mnemiopsis_leydi_ML148521a', 'MLEY_Mnemiopsis_leydi_ML042010a', 'MLEY_Mnemiopsis_leydi_ML07013a', 'MLEY_Mnemiopsis_leydi_ML05048a', 'MLEY_Mnemiopsis_leydi_ML07028a', 'MLEY_Mnemiopsis_leydi_ML078921a', 'MLEY_Mnemiopsis_leydi_ML20562a', 'MLEY_Mnemiopsis_leydi_ML01809a', 'MLEY_Mnemiopsis_leydi_ML045236a', 'MLEY_Mnemiopsis_leydi_ML039813a', 'MLEY_Mnemiopsis_leydi_ML07899a', 'MLEY_Mnemiopsis_leydi_ML078910a', 'MLEY_Mnemiopsis_leydi_ML15156a', 'MLEY_Mnemiopsis_leydi_ML17014a', 'MLEY_Mnemiopsis_leydi_ML040511a', 'MLEY_Mnemiopsis_leydi_ML218813a', 'MLEY_Mnemiopsis_leydi_ML07102a', 'MLEY_Mnemiopsis_leydi_ML40984a', 'MLEY_Mnemiopsis_leydi_ML07454a', 'MLEY_Mnemiopsis_leydi_ML45532a', 'MLEY_Mnemiopsis_leydi_ML005348a', 'MLEY_Mnemiopsis_leydi_ML065756a', 'MLEY_Mnemiopsis_leydi_ML065757a', 'MLEY_Mnemiopsis_leydi_ML01013a', 'MLEY_Mnemiopsis_leydi_ML311615a', 'MLEY_Mnemiopsis_leydi_ML07841a', 'MLEY_Mnemiopsis_leydi_ML064932a', 'MLEY_Mnemiopsis_leydi_ML351725a', 'MLEY_Mnemiopsis_leydi_ML097511a', 'MLEY_Mnemiopsis_leydi_ML030232a', 'MLEY_Mnemiopsis_leydi_ML34349a', 'MLEY_Mnemiopsis_leydi_ML059820a', 'MLEY_Mnemiopsis_leydi_ML09537a', 'MLEY_Mnemiopsis_leydi_ML16588a', 'MLEY_Mnemiopsis_leydi_ML16032a', 'MLEY_Mnemiopsis_leydi_ML16705a', 'MLEY_Mnemiopsis_leydi_ML07086a', 'MLEY_Mnemiopsis_leydi_ML259912a', 'MLEY_Mnemiopsis_leydi_ML097516a', 'MLEY_Mnemiopsis_leydi_ML138314a', 'MLEY_Mnemiopsis_leydi_ML00639a', 'MLEY_Mnemiopsis_leydi_ML08213a', 'MLEY_Mnemiopsis_leydi_ML11979a', 'MLEY_Mnemiopsis_leydi_ML28255a', 'MLEY_Mnemiopsis_leydi_ML00735a', 'MLEY_Mnemiopsis_leydi_ML010530a', 'MLEY_Mnemiopsis_leydi_ML019110a', 'MLEY_Mnemiopsis_leydi_ML06302a', 'MLEY_Mnemiopsis_leydi_ML040718a', 'MLEY_Mnemiopsis_leydi_ML06301a', 'MLEY_Mnemiopsis_leydi_ML065313a', 'MLEY_Mnemiopsis_leydi_ML22306a', 'MLEY_Mnemiopsis_leydi_ML077638a', 'MLEY_Mnemiopsis_leydi_ML12862a', 'MLEY_Mnemiopsis_leydi_ML279842a', 'MLEY_Mnemiopsis_leydi_ML21401a', 'MLEY_Mnemiopsis_leydi_ML368910a', 'MLEY_Mnemiopsis_leydi_ML190411a', 'MLEY_Mnemiopsis_leydi_ML16595a', 'MLEY_Mnemiopsis_leydi_ML065713a', 'MLEY_Mnemiopsis_leydi_ML095510a', 'MLEY_Mnemiopsis_leydi_ML327420a', 'MLEY_Mnemiopsis_leydi_ML02332a'): In-Paralogs
('STRO_XP_011676996.2', 'STRO_XP_780641.3'): In-Paralogs
```
- 1. First we count the in-paralogs by using count_inparalogs.py
```
from collections import Counter

# Read file
with open("inparalogs_after_out.txt", "r") as file:
    data = file.readlines()

# Initialize the list
sequences = []

# Loop through each line to get the seq ID
for line in data:
    # Get the part in single quotation
    start = line.find("(")
    end = line.find(")")
    
    if start != -1 and end != -1:
        # Extract the sequences and remove extra symbol
        seqs = line[start+1:end].replace("'", "").split(", ")
        # Add to the seq list
        sequences.extend(seqs)

# Count the number of occurence of each seq ID
sequence_count = Counter(sequences)

# Print out seq ID and times
for seq, count in sequence_count.items():
    print(f"{seq}: {count}")

# Export all ID into file
with open("inparalogs_4_n_index.txt", "w") as output_file:
    for seq in sequences:
        output_file.write(f"{seq}\n")

# Pring number of sequences
print(f"Total number of sequences: {len(sequences)}")
```
 Then the output should be a file called  inparalogs_4_n_index.txt and it looks like:
```
HHON_TR10104_c1_g1_i1_m_21166
HHON_TR22145_c0_g1_i1_m_43391
HHON_TR3035_c1_g3_i10_m_6505
HHON_TR25907_c2_g1_i1_m_51961
HHON_TR14192_c0_g2_i1_m_27580
HHON_TR28827_c0_g1_i1_m_61152
HHON_TR23333_c1_g1_i1_m_46165
```
- 2. Extract all the sequences above using extract_inparalogs_seq.py
```
from Bio import SeqIO
import os

# Path
index_file = "inparalogs_4_n_index.txt"
fasta_dir = "/user/work/qw23953/6_Compare_Software/8_Redo/4_all_OG_no_inparalogs/"
output_file = "inparalogs_4_n.fasta"

# Read the sequence list
with open(index_file, 'r') as f:
    target_sequences = set(line.strip() for line in f)

# Initialize the list
sequences_found = []
sequences_not_found = []

# Loop through all .fa file and extract seqs
for fasta_file in os.listdir(fasta_dir):
    if fasta_file.endswith(".fa"):
        file_path = os.path.join(fasta_dir, fasta_file)
        # Parse the FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            # If the sequence name is in the target sequence name list
            if record.id in target_sequences:
                sequences_found.append(record)
                target_sequences.remove(record.id)  # Removes the found sequence

# Save the extracted sequence to a new fasta file
with open(output_file, 'w') as output_handle:
    SeqIO.write(sequences_found, output_handle, 'fasta')

# Print the sequence name that was not found
sequences_not_found = list(target_sequences)
if sequences_not_found:
    print("Below sequence was not found:")
    for seq in sequences_not_found:
        print(seq)

# Print the total number of extracted sequences
print(f"Total number of extracted sequences: {len(sequences_found)}")
```

  Then we got an output of inparalogs_4_n.fasta, next step nwe need to map them with BUSCO sequences and keep them single copy, here we use 
  find_inparalogs_busco.py
  
```
from Bio import SeqIO
import os

# PATH
input_fasta = "inparalogs_4_n.fasta"
busco_dir = "/user/work/qw23953/6_Compare_Software/5_Sep_redo/busco_12_taxa"
output_file = "hits_match.txt"
to_delete_file = "inparalogs_to_delete.txt"

# Read 1736 sequences in inparalogs_4_n.fasta
query_sequences = {}
for record in SeqIO.parse(input_fasta, "fasta"):
    query_sequences[str(record.seq)] = record.id

# Initializes the matching record
matches = []
found_sequences = set()

# Iterate over all fasta files in the busco_12_taxa directory
for fasta_file in os.listdir(busco_dir):
    if fasta_file.endswith(".fasta"):
        file_path = os.path.join(busco_dir, fasta_file)
        for record in SeqIO.parse(file_path, "fasta"):
            seq_str = str(record.seq)
            if seq_str in query_sequences:
                # Record the matching result: busco file name and matching sequence name
                matches.append(f"{fasta_file}: {record.id}")
                found_sequences.add(seq_str)

# Generate the hits_match.txt file
with open(output_file, 'w') as f:
    for match in matches:
        f.write(match + "\n")

# Save the sequence of hits not found to inparalogs_to_delete.txt and count
not_found_count = 0
with open(to_delete_file, 'w') as f:
    for seq_str, seq_id in query_sequences.items():
        if seq_str not in found_sequences:
            f.write(f">{seq_id}\n{seq_str}\n")
            not_found_count += 1

# Count the number of sequences that find a match
found_count = len(found_sequences)
total_queries = len(query_sequences)

# Output statistics
print(f"In {total_queries} sequences, found {found_count} matchings ")
print(f"Matching results were save to {output_file} ")
print(f"Unmatching sequences were saved to{to_delete_file} ")
print(f"Number of unmatching sequences {not_found_count} ")
```
Then we got out oput of hits_match.txt and inparalogs_to_delete.txt

- 3. Delete all the inparalogs sequences above using delete_inparalogs.py
```
from Bio import SeqIO
import os

# Path
delete_fasta = "inparalogs_to_delete.txt"
target_dir = "/user/work/qw23953/6_Compare_Software/8_Redo/4_all_OG_no_inparalogs"

# Read the sequence to be deleted
sequences_to_delete = set()
for record in SeqIO.parse(delete_fasta, "fasta"):
    sequences_to_delete.add(str(record.seq))

# Count the number of deleted sequences
deleted_sequences_count = 0

# Loop through all .fa files in the directory
for fa_file in os.listdir(target_dir):
    if fa_file.endswith(".fa"):
        file_path = os.path.join(target_dir, fa_file)
        temp_file_path = file_path + ".temp"

        # Read the original file and writes to the temporary file, excluding sequences to be deleted
        with open(file_path, 'r') as infile, open(temp_file_path, 'w') as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                if str(record.seq) not in sequences_to_delete:
                    SeqIO.write(record, outfile, 'fasta')
                else:
                    deleted_sequences_count += 1

        # Replace original file
        os.replace(temp_file_path, file_path)

# Output
print(f"Number of deleted sequences is: {deleted_sequences_count}")
```
   - 4. Then, a lot of OrthoGroups become null file as they have no sequences remained, we move them to null_file folder using a script move_null_file.sh
   
```
#!/bin/bash

source_folder="/user/work/qw23953/6_Compare_Software/8_Redo/4_all_OG_no_inparalogs"
target_folder="${source_folder}/null_file"


mkdir -p "$target_folder"

# Find empty files and move to the destination folder
find "$source_folder" -type f -empty -exec mv {} "$target_folder" \;

echo "All empty files have been moved to $target_folder."
```


   
5. To now, for 13594 OGs, 3972 are out-paralog OGs, after removing in-paralog sequences, we lost 3993 OGs, so 5629 OGs are clean to use after filtering, they are in a folder called OF_4_N.
   
   + For now, we get the most relaxed level of single-copy Orthologues Groups.
   + Before we filtering, we also have OF_4_12 using BUSCO definition of single copy orthologs, the number of OF_4_12 is 1659 OGs,
   + We also have OrthoFinder definition of 365 Single_Copy_Orthologue_Sequences.

  **Using these diferent level of restrictions we make conparison with BUSCO output, here the BUSCO output was reformatted according to gene families.**
   

  main comparison script compare.py
```
import os

# Specify the path to the BUSCO and Orthofinder folders
busco_folder = '/user/work/qw23953/6_Compare_Software/4_Aug_redo/busco_12_taxa'
orthofinder_folder = '/user/work/qw23953/6_Compare_Software/8_Redo/5_all_filtered'

# Create output file
output_file_path = '/user/work/qw23953/6_Compare_Software/8_Redo/OF_all_results.txt'

# Open the output file and prepare to write the comparison results
with open(output_file_path, 'w') as output_file:
    output_file.write("busco_file_name\tSeq_ID\torthofinder_file_name\tSeq_ID\tlength of sites identitical\n")
    
    # Loop through each fasta file in the Orthofinder folder
    for orthofinder_file in os.listdir(orthofinder_folder):
        if orthofinder_file.endswith(".fa"):  # Make sure only process fasta file
            orthofinder_filepath = os.path.join(orthofinder_folder, orthofinder_file)
            
            # Extract the file name of the Orthofinder file as an identifier
            orthofinder_filename = orthofinder_file.split('.')[0]
            
            # Read the sequence information in the Orthofinder file
            orthofinder_sequences = {}
            with open(orthofinder_filepath, 'r') as orthofinder_handle:
                current_sequence = ""
                for line in orthofinder_handle:
                    if line.startswith(">"):
                        current_sequence = line.strip().split(" ")[0]  # Extract sequence ID
                        orthofinder_sequences[current_sequence] = ""
                    else:
                        orthofinder_sequences[current_sequence] += line.strip()
            
            # Loop through each fasta file in the BUSCO folder
            for busco_file in os.listdir(busco_folder):
                if busco_file.endswith(".fasta"):  # Only fasta file
                    busco_filepath = os.path.join(busco_folder, busco_file)
                    
                    # Extract the file name of the Orthofinder file as an identifier
                    busco_filename = busco_file.split('.')[0]
                    
                    # Read the sequence information in the BUSCO file
                    busco_sequences = {}
                    with open(busco_filepath, 'r') as busco_handle:
                        current_sequence = ""
                        for line in busco_handle:
                            if line.startswith(">"):
                                current_sequence = line.strip().split(" ")[0]  # Extract sequence ID
                                busco_sequences[current_sequence] = ""
                            else:
                                busco_sequences[current_sequence] += line.strip()
                    
                    # Check if each sequence in the Orthofinder file exists in each sequence in the BUSCO file
                    for orthofinder_id, orthofinder_seq in orthofinder_sequences.items():
                        for busco_id, busco_seq in busco_sequences.items():
                            if orthofinder_seq in busco_seq:
                                # Write the result to a txt file
                                output_file.write(f"{busco_filename}\t{busco_id}\t{orthofinder_filename}\t{orthofinder_id}\t{len(orthofinder_seq)}\n")

print("Comparison results have been written to", output_file_path)
```



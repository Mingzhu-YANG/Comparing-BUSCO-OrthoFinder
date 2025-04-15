import os
import shutil
from collections import Counter
from Bio import SeqIO

# ========== CONFIGURATION AREA ==========
base_dir = os.getcwd()
busco_dir = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/3_BUSCO_2025_processed/final"
outparalog_family_dir = os.path.join(base_dir, "outparalog_family")
backup_dir = os.path.join(base_dir, "with_inparalogs_copy")
null_dir = os.path.join(base_dir, "null_file")

# ========== Step 1: Merge outparalog family files ==========
os.system("ls *_outparalogs.txt > outparalog_family.txt")

# ========== Step 2: Move all related outparalog files ==========
os.makedirs(outparalog_family_dir, exist_ok=True)
with open("outparalog_family.txt", "r") as f:
    for line in f:
        line = line.strip()
        if line.endswith("_outparalogs.txt"):
            prefix = line.replace("_tree_outparalogs.txt", "")
            for suffix in ["_tree.txt", "_tree_outparalogs.txt", "_tree_inparalogs.txt", ".fa"]:
                filename = f"{prefix}{suffix}"
                src = os.path.join(base_dir, filename)
                dst = os.path.join(outparalog_family_dir, filename)
                if os.path.exists(src):
                    shutil.move(src, dst)
                    print(f"Outparalogs Moved: {filename}")
                else:
                    print(f"Inparalogs in Outparalos File not found: {filename}")

# ========== Step 3: Merge inparalog txt and extract sequences ==========
os.system("cat *_inparalogs.txt > inparalogs_after_out.txt")
sequences = []
with open("inparalogs_after_out.txt") as f:
    for line in f:
        if "(" in line and ")" in line:
            seqs = line.split("(")[-1].split(")")[0].replace("'", "").split(", ")
            sequences.extend(seqs)
Counter(sequences)
with open("inparalogs_4_n_index.txt", "w") as f:
    for s in sequences:
        f.write(f"{s}\n")

# ========== Step 4: Extract short ID ==========
with open("inparalogs_4_n_index.txt") as f, open("inparalogs_4_n_index_extracted.txt", "w") as out:
    for line in f:
        parts = line.strip().split("_")
        extracted = parts[2]
        if len(parts) > 3:
            extracted += "_" + "_".join(parts[3:])
        out.write(extracted + "\n")

# ========== Step 5: Extract matching sequences ==========
with open("inparalogs_4_n_index_extracted.txt") as f:
    target_ids = set(line.strip() for line in f)

matched_seqs = []
for file in os.listdir(base_dir):
    if file.endswith(".fa"):
        for record in SeqIO.parse(file, "fasta"):
            if record.id in target_ids:
                matched_seqs.append(record)
                target_ids.remove(record.id)

SeqIO.write(matched_seqs, "inparalogs_4_n.fasta", "fasta-2line")
print(f"✅ Extracted {len(matched_seqs)} sequences")
if target_ids:
    print("⚠️ Not found:")
    print("\n".join(target_ids))

# ========== Step 6: Compare with BUSCO database ==========
input_fasta = "inparalogs_4_n.fasta"
output_file = "hits_match.txt"
to_delete_file = "inparalogs_to_delete.txt"

query_sequences = {}
for record in SeqIO.parse(input_fasta, "fasta"):
    query_sequences[str(record.seq)] = record.id

matches = []
found_sequences = set()
for fasta_file in os.listdir(busco_dir):
    if fasta_file.endswith(".faa"):
        file_path = os.path.join(busco_dir, fasta_file)
        for record in SeqIO.parse(file_path, "fasta"):
            seq_str = str(record.seq)
            if seq_str in query_sequences:
                matches.append(f"{fasta_file}: {record.id}")
                found_sequences.add(seq_str)

with open(output_file, 'w') as f:
    for match in matches:
        f.write(match + "\n")

not_found_count = 0
with open(to_delete_file, 'w') as f:
    for seq_str, seq_id in query_sequences.items():
        if seq_str not in found_sequences:
            f.write(f">{seq_id}\n{seq_str}\n")
            not_found_count += 1

print(f"📊 Found {len(found_sequences)} matches out of {len(query_sequences)} sequences")
print(f"📁 Match results: {output_file}")
print(f"📁 To delete: {to_delete_file}, count = {not_found_count}")

# ========== Step 7: Backup .fa files before deletion ==========
os.makedirs(backup_dir, exist_ok=True)
for file in os.listdir(base_dir):
    if file.endswith(".fa"):
        shutil.copy(os.path.join(base_dir, file), os.path.join(backup_dir, file))
print(f"🗂️ Backup complete: all .fa files saved to {backup_dir}")

# ========== Step 8: Delete unmatched sequences ==========
seqs_to_delete = {str(r.seq) for r in SeqIO.parse("inparalogs_to_delete.txt", "fasta")}
deleted = 0
for file in os.listdir(base_dir):
    if file.endswith(".fa"):
        input_file = os.path.join(base_dir, file)
        temp_file = input_file + ".tmp"
        with open(input_file) as infile, open(temp_file, "w") as outfile:
            for rec in SeqIO.parse(infile, "fasta"):
                if str(rec.seq) not in seqs_to_delete:
                    SeqIO.write(rec, outfile, "fasta")
                else:
                    deleted += 1
        os.replace(temp_file, input_file)
print(f"🗑️ Deleted {deleted} unmatched sequences")

# ========== Step 9: Move empty files ==========
os.makedirs(null_dir, exist_ok=True)
for file in os.listdir(base_dir):
    if file.endswith(".fa") and os.path.getsize(os.path.join(base_dir, file)) == 0:
        shutil.move(os.path.join(base_dir, file), os.path.join(null_dir, file))
print(f"📂 All empty files moved to: {null_dir}")

# ========== Step 10: Organize final outputs and print summary ==========
final_dir = os.path.join(base_dir, "final_OGs_no_paralogs")
intermediate_dir = os.path.join(base_dir, "intermediate")

os.makedirs(final_dir, exist_ok=True)
os.makedirs(intermediate_dir, exist_ok=True)

# Move final .fa files
for file in os.listdir(base_dir):
    if file.endswith(".fa"):
        shutil.move(os.path.join(base_dir, file), os.path.join(final_dir, file))

# Move .txt and other intermediate files
for file in os.listdir(base_dir):
    if file.endswith(".txt") or file.endswith(".fasta") or file.endswith(".log"):
        shutil.move(os.path.join(base_dir, file), os.path.join(intermediate_dir, file))

# Count final .fa files and total sequences
fa_files = [f for f in os.listdir(final_dir) if f.endswith(".fa")]
file_count = len(fa_files)
seq_count = 0
for fa in fa_files:
    with open(os.path.join(final_dir, fa)) as f:
        seq_count += sum(1 for line in f if line.startswith(">"))

print(f"✅ Final cleaned OGs moved to: {final_dir}")
print(f"📦 Intermediate files moved to: {intermediate_dir}")
print(f"📊 Final summary:")
print(f"   - Number of final clean OGs: {file_count}")
print(f"   - Total number of sequences: {seq_count}")

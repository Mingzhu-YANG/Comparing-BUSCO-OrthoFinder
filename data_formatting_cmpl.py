import os
import re
import argparse
import logging
import shutil
from Bio import SeqIO
from collections import defaultdict

VERSION = "data_formatting_cmpl.py v1.1.0"

def setup_logger(base_dir):
    log_file = os.path.join(base_dir, "pipeline.log")
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger().addHandler(console)
    return log_file

def rename_sequences(base_dir):
    for folder in os.listdir(base_dir):
        if folder.endswith(".out"):
            folder_name = folder.replace(".fasta.out", "")
            faa_dir = os.path.join(
                base_dir, folder,
                "run_metazoa_odb12", "busco_sequences",
                "single_copy_busco_sequences"
            )

            if not os.path.exists(faa_dir):
                logging.warning(f"Directory not found: {faa_dir}")
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
                    logging.info(f"Renamed sequences in file: {faa_path}")
    logging.info("All faa sequences are renamed.")

def merge_sequences(base_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for folder in os.listdir(base_dir):
        if folder.endswith(".out"):
            source_dir = os.path.join(
                base_dir, folder,
                "run_metazoa_odb12/busco_sequences/single_copy_busco_sequences"
            )
            if not os.path.exists(source_dir):
                logging.warning(f"Source dir not found: {source_dir}, skipping.")
                continue

            output_filename = f"{folder.replace('.fasta.out', '')}.faa"
            output_file = os.path.join(output_dir, output_filename)

            with open(output_file, "w") as outfile:
                for faa_file in os.listdir(source_dir):
                    if faa_file.endswith(".faa"):
                        faa_path = os.path.join(source_dir, faa_file)
                        for record in SeqIO.parse(faa_path, "fasta"):
                            SeqIO.write(record, outfile, "fasta")
            logging.info(f"Merged into {output_file}")
    logging.info("All faa files merged.")

def group_by_gene_family(input_dir, output_dir):
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
                    logging.warning(f"Unexpected sequence format skipped: {record.id}")

    for family_id, sequences in gene_family_dict.items():
        output_file = os.path.join(output_dir, f"{family_id}at33208.faa")
        with open(output_file, "w") as f:
            SeqIO.write(sequences, f, "fasta")
        logging.info(f"Gene family {family_id} written to {output_file}")
    logging.info("Sequences grouped by gene family.")

def extract_bitscores(base_dir, output_dir, bitscore_file):
    os.makedirs(output_dir, exist_ok=True)
    with open(bitscore_file, "w") as final_out:
        for folder in os.listdir(base_dir):
            if folder.endswith(".out"):
                species_name = folder.replace(".fasta.out", "")
                results_dir = os.path.join(
                    base_dir, folder,
                    "run_metazoa_odb12", "hmmer_output", "initial_run_results"
                )
                output_file = os.path.join(output_dir, f"{species_name}_scores.txt")

                if not os.path.exists(results_dir):
                    logging.warning(f"Results dir not found: {results_dir}, skipping {species_name}")
                    continue

                with open(output_file, "w") as outfile:
                    for out_file in os.listdir(results_dir):
                        if out_file.endswith("at33208.out"):
                            out_path = os.path.join(results_dir, out_file)
                            with open(out_path) as f:
                                lines = f.readlines()
                                if len(lines) >= 4:
                                    parts = lines[3].split()
                                    if len(parts) >= 8:
                                        gene, col4, col8 = parts[0], parts[3], parts[7]
                                        line_out = f"{species_name} {gene} {col4} {col8}\n"
                                        outfile.write(line_out)
                                        final_out.write(line_out)
                logging.info(f"Bitscores extracted for {species_name}")
    logging.info(f"All species processed! Final results saved in {bitscore_file}")

# ===== 新增 3_rename.py 功能 =====
def read_bitscore(bitscore_file):
    bitscore_dict = {}
    with open(bitscore_file, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 4:
                species, gene_name, score = parts[0], parts[2], parts[3]
                key = (species, gene_name)
                bitscore_dict[key] = score
    return bitscore_dict

def rename_sequences_with_bitscore(fasta_directory, bitscore_file):
    bitscore_dict = read_bitscore(bitscore_file)
    for filename in os.listdir(fasta_directory):
        if filename.endswith(".faa"):
            fasta_file = os.path.join(fasta_directory, filename)
            modified_records = []
            for record in SeqIO.parse(fasta_file, "fasta"):
                match = re.match(r"(.+?)_(\d+)at33208", record.id)
                if match:
                    species, gene_name = match.group(1), match.group(2)
                    key = (species, gene_name)
                    if key in bitscore_dict:
                        score = bitscore_dict[key]
                        new_id = f"{record.id}_@-{score}"
                        record.id = new_id
                        record.description = ""
                    else:
                        logging.warning(f"Bitscore not found for {key}")
                    modified_records.append(record)
                else:
                    logging.warning(f"Regex pattern mismatch for {record.id}")

            
            SeqIO.write(modified_records, fasta_file, "fasta")
            logging.info(f"Updated FASTA in place: {fasta_file}")
# =================================

def main():
    parser = argparse.ArgumentParser(
        description="BUSCO sequence and bitscore processing pipeline (part 1)."
    )
    parser.add_argument("--base_dir", required=True,
                        help="Base directory containing BUSCO .out folders")
    parser.add_argument("--keep_intermediates", action="store_true",
                        help="Keep intermediate merged .faa files (default: False)")
    parser.add_argument("--keep_scores", action="store_true",
                        help="Keep scores directory after generating bitscore_all.txt (default: False)")
    args = parser.parse_args()

    base_dir = args.base_dir
    log_file = setup_logger(base_dir)

    intermediate_dir = os.path.join(base_dir, "intermediate")
    final_output_dir = os.path.join(base_dir, "1_final_60_672")
    scores_dir = os.path.join(base_dir, "2_get_scores")
    bitscore_file = os.path.join(base_dir, "bitscore_all.txt")
    ready_dir = os.path.join(base_dir, "ready_data")

    logging.info(f"=== data_formatting_part1.py {VERSION} started ===")
    logging.info(f"Parameters: base_dir={base_dir}, "
                 f"keep_intermediates={args.keep_intermediates}, "
                 f"keep_scores={args.keep_scores}")

    # Step 1: rename sequences
    rename_sequences(base_dir)

    # Step 2: merge sequences
    merge_sequences(base_dir, intermediate_dir)

    # Step 3: group sequences
    group_by_gene_family(intermediate_dir, final_output_dir)

    # Step 4: clean up intermediate
    if not args.keep_intermediates:
        for file in os.listdir(intermediate_dir):
            os.remove(os.path.join(intermediate_dir, file))
            logging.info(f"Deleted intermediate file: {file}")
        os.rmdir(intermediate_dir)
        logging.info("Intermediate files removed.")

    # Step 5: extract bitscores
    extract_bitscores(base_dir, scores_dir, bitscore_file)

    # Step 6: remove scores_dir if not kept
    if not args.keep_scores:
        for file in os.listdir(scores_dir):
            os.remove(os.path.join(scores_dir, file))
            logging.info(f"Deleted score file: {file}")
        os.rmdir(scores_dir)
        logging.info("Scores directory removed, only bitscore_all.txt kept.")

    # Step 7: copy results to ready_data
    os.makedirs(ready_dir, exist_ok=True)
    # copy 1_final_60_672
    final_ready = os.path.join(ready_dir, "1_final_60_672")
    if os.path.exists(final_ready):
        shutil.rmtree(final_ready)
    shutil.copytree(final_output_dir, final_ready)
    logging.info(f"Copied {final_output_dir} -> {final_ready}")

    # copy bitscore_all.txt
    bitscore_ready = os.path.join(ready_dir, "2_bitscore_all.txt")
    shutil.copy(bitscore_file, bitscore_ready)
    logging.info(f"Copied {bitscore_file} -> {bitscore_ready}")

    # Step 8: run rename with bitscore (3_rename.py)
    rename_sequences_with_bitscore(final_ready, bitscore_ready)

    logging.info(f"=== data_formatting_part1.py {VERSION} finished successfully ===")
    print(f"Pipeline finished. Log file written to {log_file}")

if __name__ == "__main__":
    main()

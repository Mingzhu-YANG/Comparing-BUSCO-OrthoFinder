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


***get_sequences.py***



### S2.1-Copy the Orthogroup_Sequences directory to the working directory (now we work on OrthoFinder output).

1. Using the bash script to group all OGs into different categories, for the Orthogroup_Sequences dir, we have 71249 OGs, we then group them into 2:
   
   - folder OF_1_3: 57459 OGs, using script **get_groups_1_3_seqs.sh**
     these are OGs with less than 3 sequences, cannot be used to build trees, so we have to take them out to folder OF_1_3_seqs.
   - folder OF_4_N: 13790 OGs, we separate them using script **find_single_copy.sh**
       - we further select the OGs with 4-12 sequences without duplication, there are 1659 OGs
       - so the rest 12131 OGs with duplications
       - we can also use **get_strict_single_copy.sh** to get single-copy orthogroups with all taxa present

**get_groups_1_3_seqs.sh**


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
 
### S2.2-topology analysis (identify paralogs, i.e. in-paralogs and out-paralogs)
1. identify_paralogs_2.py


2. filtering the paralogs
   
   
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



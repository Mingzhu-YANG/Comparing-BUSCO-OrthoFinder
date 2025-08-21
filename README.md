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


2. Copy treefiles using copy_treefiles.sh

```
#!/bin/bash

# Set the directories
src_dir=/home/qw23953/mingzhu/compare_BUSCO_OrthoFinder/5_BUSCO_OrthoFinder_2024/2_OrthoFinder/Results_Mar11/Gene_Trees
dst_dir=/home/qw23953/mingzhu/compare_BUSCO_OrthoFinder/5_BUSCO_OrthoFinder_2024/2_OrthoFinder/Orthogroup_Sequences/OF_4_N

# Loop through the txt files in the destination directory
for file in "$dst_dir"/*.txt; do
  # Create a new directory for each file
  mkdir -p "${file%.txt}"

  # Loop through the contents of the file
  while IFS= read -r line; do
    # Extract the OGxxx part
    og_name="${line%.fa}"

    # Copy the corresponding tree file
    cp "$src_dir/$og_name"_tree.txt "${file%.txt}/"
  done < "$file"
done

```


 
### S2.2-Topology analysis (identify paralogs, i.e. in-paralogs and out-paralogs)
1. identify_paralogs_2.py


2. filtering the paralogs
   
   
5. To now, for 13594 OGs, 3972 are out-paralog OGs, after removing in-paralog sequences, we lost 3993 OGs, so 5629 OGs are clean to use after filtering, they are in a folder called OF_4_N.
   
   + For now, we get the most relaxed level of single-copy Orthologues Groups.
   + Before we filtering, we also have OF_4_12 using BUSCO definition of single copy orthologs, the number of OF_4_12 is 1659 OGs,
   + We also have OrthoFinder definition of 365 Single_Copy_Orthologue_Sequences.
   
![BUSCO filtering criteria and overlap](/busco_pie_chart.png)

  ### S2.3-Using these different level of restrictions we make comparison with BUSCO output, here the BUSCO output was reformatted according to gene families.**
   

  main comparison script 
  compare.py


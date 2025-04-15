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
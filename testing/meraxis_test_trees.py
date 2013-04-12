#!/usr/bin/env python

"""Doc-string here..."""

import numpy as np
import matplotlib.pyplot as plt

__author__ = "Simon Mutch"
__date__   = "12 Apr 2013"

catalog_desc = [
    ("id_MBP", 'q'),
    ("M_vir", 'f8'),
    ("n_particles", 'i4'),
    ("position_COM", ('f4', 3)),
    ("position_MBP", ('f4', 3)),
    ("velocity_COM", ('f4', 3)),
    ("velocity_MBP", ('f4', 3)),
    ("R_vir", 'f4'),
    ("R_halo", 'f4'),
    ("R_max", 'f4'),
    ("V_max", 'f4'),
    ("sigma_v", 'f4'),
    ("spin", ('f4', 3),),
    ("q_triaxial", 'f4'),
    ("s_triaxial", 'f4'),
    ("shape_eigen_vectors", ('f4', (3,3))),
    ("padding", 'S8'),
]

group_desc = [
    ("id", 'i4'),
    ("tree_flags", 'i4'),
    ("desc_id", 'i4'),
    ("dummy", 'i4'),
    ("file_offset", 'i4'),
    ("file_index", 'i4'),
    ("n_subgroups", 'i4'),
]

subgroup_desc = [
    ("id", 'i4'),
    ("tree_flags", 'i4'),
    ("desc_id", 'i4'),
    ("dummy", 'i4'),
    ("file_offset", 'i4'),
    ("file_index", 'i4'),
]

trees_header_desc = [
  ("n_groups", 'i4'),
  ("n_subgroups", 'i4'),
  ("n_halos_max", 'i4'),
  ("n_trees_subgroup", 'i4'),
  ("n_trees_group", 'i4'),
]

# Create the trees files
header = np.zeros(1, dtype=trees_header_desc)
header['n_groups'] = 3
header['n_subgroups'] = 3
header['n_halos_max'] = 6
header['n_trees_subgroup'] = 3
header['n_trees_group'] = 3

group = np.zeros(3, dtype=group_desc)
group['id']=[0,1,2]
group['tree_flags']=[1,1,1]
group['desc_id']=[0,1,1]
group['dummy']=[0,0,0]
group['file_offset']=[1,1,1]
group['file_index']=[0,1,1]
group['n_subgroups']=[1,1,1]

subgroup = np.zeros(3, dtype=subgroup_desc)
subgroup['id']=[0,1,2]
subgroup['tree_flags']=[1,1,1]
subgroup['desc_id']=[0,1,2]
subgroup['dummy']=[0,0,0]
subgroup['file_offset']=[1,1,1]
subgroup['file_index']=[0,1,2]

with open('test_step001_scan001.trees_horizontal.0', "wb") as fout:
    for v in header:
        v.tofile(fout)
    for g, s in zip(group, subgroup):
        for v in g:
            v.tofile(fout)
        for v in s:
            v.tofile(fout)

header = np.zeros(1, dtype=trees_header_desc)
header['n_groups'] = 2
header['n_subgroups'] = 3
header['n_halos_max'] = 6
header['n_trees_subgroup'] = 3
header['n_trees_group'] = 2

group = np.zeros(2, dtype=group_desc)
group['id']=[0,1]
group['tree_flags']=[1,2]
group['desc_id']=[0,1]
group['dummy']=[0,0]
group['file_offset']=[1,1]
group['file_index']=[0,1]
group['n_subgroups']=[1,2]

subgroup = np.zeros(3, dtype=subgroup_desc)
subgroup['id']=[0,1,2]
subgroup['tree_flags']=[1,1,1]
subgroup['desc_id']=[0,1,1]
subgroup['dummy']=[0,0,0]
subgroup['file_offset']=[1,1,1]
subgroup['file_index']=[0,1,1]

with open('test_step001_scan001.trees_horizontal.1', "wb") as fout:
    for v in header:
        v.tofile(fout)
    for v in group[0]:
        v.tofile(fout)
    for v in subgroup[0]:
        v.tofile(fout)
    for v in group[1]:
        v.tofile(fout)
    for v in subgroup[1]:
        v.tofile(fout)
    for v in subgroup[2]:
        v.tofile(fout)

header = np.zeros(1, dtype=trees_header_desc)
header['n_groups'] = 2
header['n_subgroups'] = 2
header['n_halos_max'] = 6
header['n_trees_subgroup'] = 2
header['n_trees_group'] = 2

group = np.zeros(2, dtype=group_desc)
group['id']=[0,1]
group['tree_flags']=[1,1]
group['desc_id']=[0,1]
group['dummy']=[0,0]
group['file_offset']=[1,1]
group['file_index']=[0,1]
group['n_subgroups']=[1,1]

subgroup = np.zeros(2, dtype=subgroup_desc)
subgroup['id']=[0,1]
subgroup['tree_flags']=[1,2]
subgroup['desc_id']=[0,1]
subgroup['dummy']=[0,0]
subgroup['file_offset']=[1,1]
subgroup['file_index']=[0,1]

with open('test_step001_scan001.trees_horizontal.2', "wb") as fout:
    for v in header:
        v.tofile(fout)
    for g, s in zip(group, subgroup):
        for v in g:
            v.tofile(fout)
        for v in s:
            v.tofile(fout)



# Now deal with the catalog files

group = np.zeros(3, dtype=catalog_desc)
group["id_MBP"] = [1,500,1000]
group["M_vir"] = [5e9,5e9,5e9]
group["n_particles"] = [100,100,100]
group["R_vir"] = [2, 2, 2]
group["V_max"] = [100, 100, 100]
with open('subfind_000.catalog_group_properties', "wb") as fout:
    for g in group:
        g.tofile(fout)

subgroup = np.zeros(3, dtype=catalog_desc)
subgroup = group.copy()
with open('subfind_000.catalog_subgroup_properties', "wb") as fout:
    for s in subgroup:
        s.tofile(fout)

group = np.zeros(2, dtype=catalog_desc)
group["id_MBP"] = [1,500]
group["M_vir"] = [5e9,10e9]
group["n_particles"] = [100,200]
group["R_vir"] = [2, 3]
group["V_max"] = [100, 200]
with open('subfind_001.catalog_group_properties', "wb") as fout:
    for g in group:
        g.tofile(fout)

with open('subfind_001.catalog_subgroup_properties', "wb") as fout:
    for s in subgroup:
        s.tofile(fout)

with open('subfind_002.catalog_group_properties', "wb") as fout:
    for g in group:
        g.tofile(fout)

subgroup = np.zeros(3, dtype=catalog_desc)
subgroup["id_MBP"] = [1,500]
subgroup["M_vir"] = [5e9,10e9]
subgroup["n_particles"] = [100,200]
subgroup["R_vir"] = [2, 3]
subgroup["V_max"] = [100, 200]
with open('subfind_002.catalog_subgroup_properties', "wb") as fout:
    for s in subgroup:
        s.tofile(fout)


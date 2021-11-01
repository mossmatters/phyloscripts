import sys
from ete3 import Tree

if len(sys.argv) < 3:
    print("Usage: python reroot_trees.py treefile.tre outgroup.list > rerooted.tre")
    
# outgroup.list should be one sample name per line
outgroup_names = [x.rstrip() for x in open(sys.argv[2])]

for line in open(sys.argv[1]):
    t = Tree(line.rstrip())
    outgroups_in_tree = list(set(t.get_leaf_names()).intersection(set(outgroup_names)))
    if len(outgroups_in_tree) > 1:
        ancestor = t.get_common_ancestor(outgroups_in_tree)
        if ancestor == t:
            ingroups_in_tree = list(set(t.get_leaf_names()).difference(set(outgroups_in_tree)))
            ancestor = t.get_common_ancestor(ingroups_in_tree)
            t.set_outgroup(ancestor)
            print(t.write())
        else:
            t.set_outgroup(ancestor)
            print(t.write())
    elif len(outgroups_in_tree) == 1:
        t.set_outgroup(outgroups_in_tree[0])
        print(t.write())
    else:
        continue

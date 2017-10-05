#!/usr/bin/env python

helptext='''This script identfies branch length outliers on a phylogeny. An outlier is a 
branch with a length that exceeds a percentage of the maximum depth of the tree.

The input is a file containing one tree in newick format, and an optional file containing
a list of outgroup taxa (one per line). 

The output will be an ASCII depiction of each branch with a length exceeding the threshold
(default 25% for ingroups, 75% for outgroups). A PNG file will also be generated for the 
tree, with outgroup taxa in blue and branch length outliers in red.

Dependencies: 
Python > 2.7
ETE3 installed with all graphical dependencies
'''

#Given a tree, determine if there are long branch lengths

import sys, argparse, os
from ete3 import Tree,TreeStyle,TextFace,NodeStyle


# outgroups=set([x.rstrip() for x in open(sys.argv[2])])
# outgroups_in_tree = list(set(t.get_leaf_names()).intersection(set(outgroups)))
# ingroups_in_tree = list(set(t.get_leaf_names()).difference(set(outgroups)))
# 
# if len(outgroups_in_tree) > 1:
#     ancestor = t.get_common_ancestor(outgroups_in_tree)
#     try:
#         t.set_outgroup(ancestor)
#         #ingroup_monophyly = t.check_monophyly(ingroups_in_tree,"name")
#         #if not ingroup_monophyly[0]:
#         #    sys.stdout.write("Ingroup polyphyletic! for {}\n".format(sys.argv[1]))
#         #    print(ingroup_monophyly)
#     except:
#         sys.stdout.write("Ingroup not monophyletic for {}!\n".format(sys.argv[1]))
#         sys.exit(1)
#     #print t.write()
# elif len(outgroups_in_tree) == 1:
#     t.set_outgroup(outgroups_in_tree[0])
#     ancestor = t.get_leaves_by_name(outgroups_in_tree[0])[0]
#     #print t.write()
# else:
#     sys.stdout.write("no outgroups found for {}!\n".format(sys.argv[1]))
#     sys.exit(1)

#ingroup = ancestor.get_sisters()[0].detach()
#ingroup_depth = ingroup.get_farthest_node()[1]



#print(outgroups)
def get_bad_nodes(t,inlen,outlen,leaflen,outgroups=None):
    bad_nodes = []
    tree_depth = t.get_farthest_node()[1]
    for node in t.traverse():
        isOutgroup = True
        for leaf in node.get_leaves():
            if outgroups:
                if leaf.name not in outgroups:
                    isOutgroup = False
            else:
                isOutgroup=False
        if isOutgroup:
            if node.dist > tree_depth * outlen:
                print(node)
                bad_nodes.append(node)
        elif node.is_leaf():
            if node.dist > tree_depth * leaflen:
                print(node)
                bad_nodes.append(node)
        elif node.dist > tree_depth * inlen:
            print(node)
            bad_nodes.append(node)
    return bad_nodes        

def make_png(t,bad_nodes,png_name,outgroups=None):
    for n in t.traverse():
        n.img_style["size"]=0
        if n in bad_nodes:
            nstyle = NodeStyle()
            nstyle["hz_line_color"] = "red"
            nstyle["hz_line_width"] = 3
            n.set_style(nstyle)
        if n.is_leaf():
            if outgroups:
                if n.name in outgroups:
                    name_face = TextFace(n.name, fgcolor="blue")
                else:
                    name_face = TextFace(n.name, fgcolor="black")
            else:
                name_face = TextFace(n.name, fgcolor="black")    
            n.add_face(name_face,0,"branch-right")

    if len(bad_nodes) > 0:
        ts = TreeStyle()
        ts.show_leaf_name = False
        gene_name = png_name
        ts.title.add_face(TextFace(gene_name,fsize=15,bold=True),0)
        my_png = t.render(png_name,tree_style=ts)     

def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("treefile",help="File containing one tree in newick format")
    parser.add_argument("--outgroups",help="File containing list of outgroup taxa, one per line")
    parser.add_argument("--png",help="Name of png file, default is same as tree file name",default=None)
    parser.add_argument("--inlen",help="Percentage of max tree depth for ingroup outliers default = %(default)s",default=0.25,type=float)
    parser.add_argument("--outlen",help="Percentage of max tree depth for outgroup outliers default = %(default)s",default=0.75,type=float)
    parser.add_argument("--leaflen",help="Percentage of max tree depth for leaf outliers default = %(default)s",default=0.25,type=float)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    if os.path.isfile(args.treefile):
        t = Tree(args.treefile)
    else:
        print("Treefile {} not found!\n".format(args.treefile))
        sys.exit(1)    
    if args.png:
        if args.png.endswith(".png"):
            png_name = args.png
        else:
            png_name =  args.png + ".png"    
    else:
        png_name = os.path.basename(args.treefile).split(".")[0] + ".png"
    
    if args.outgroups:
        outgroups = set([x.rstrip() for x in open(args.outgroups)])
    else:
        outgroups = None
    
    bad_nodes = get_bad_nodes(t,args.inlen,args.outlen,args.leaflen,outgroups=outgroups)
    if len(bad_nodes) > 0:
        make_png(t,bad_nodes,png_name,outgroups=outgroups)
    #else:
        #print("No outliers found for {}\n".format(os.path.basename(args.treefile)))
        
    
    
    
    
    
    
    
    
if __name__ == "__main__":main()    
    
    



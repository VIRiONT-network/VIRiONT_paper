import os
os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree, TreeStyle, NodeStyle , TextFace
import sys

newick_tree=sys.argv[1]
output_pdf=sys.argv[2]
# Load a tree structure from a newick file.
#t = Tree("/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/TEST_RESULTS_MI/10_QC_ANALYSIS/TREE/IQtree_analysis.treefile")
t = Tree(newick_tree)

#print(t)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False
#ts.scale =  20
#ts.branch_vertical_margin = 1 
ts.mode = "c" # draw tree in circular mode
#ts.title.add_face(TextFace("Hello ETE", fsize=20), column=0)
# Draws nodes as small red spheres of diameter equal to 10 pixels
nstyle = NodeStyle()
nstyle["shape"] = "sphere"
nstyle["size"] = 10
nstyle["fgcolor"] = "darkred"
# Applies the same static style to all nodes in the tree. 
for n in t.traverse():
   n.set_style(nstyle)

#t.render("/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/TEST_RESULTS_MI/10_QC_ANALYSIS/TREE/mytree.pdf", w=183, units="mm", tree_style=ts)
t.render(output_pdf, w=183, units="mm", tree_style=ts)
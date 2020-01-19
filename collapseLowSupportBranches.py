import dendropy
import sys
import random

def collapseLowSupportBranches(treeFile):
    tree = dendropy.Tree.get(path=treeFile,schema='newick')
    for e in tree.postorder_edge_iter():
        if e.head_node.label != None:
            if int(e.head_node.label) < 50:
                e.collapse()
    fileSplit = treeFile.split('.')
    outfile = open(fileSplit[0] + '.collapsed.trees','w')
    outfile.write(tree.as_string(schema='newick'))
    outfile.close()

infile = open(sys.argv[1])
for line in infile:
    realLine = line
    while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\n':
        realLine = realLine[0:-1]
    collapseLowSupportBranches(realLine)

from Bio import Phylo
import dendropy
from cStringIO import StringIO
import sys
import random
'''USAGE:

python subTreeIterator.py trees.fofn Sp1ID Sp2ID SpNID > summaryFile.txt'''


def subTreeIterator(treeFile,speciesIDList):
    tree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
    treeString = tree.as_string('newick')
    fileSplit = treeFile.split('.')
    og = fileSplit[0]
    ogNum = -1
    #seqDict = buildSeqDict(og + '.aligned.fasta')
    totalLeaves = []
    outgroupLeaves = []
    ingroupList = speciesIDList[1:]
    outgroupID = speciesIDList[0]
    ingroupLeaves = []
    for taxon in tree.taxon_namespace:
        totalLeaves.append(taxon.label)
        if speciesIDList[0] in taxon.label:
            outgroupLeaves.append(taxon.label)
        else:
            ingroupLeaves.append(taxon.label)
    speciesLeafDict,numSpecies = checkLeaves(totalLeaves,speciesIDList[1:])
    trimLeaves = []
    if numSpecies == len(speciesIDList[1:]):
        subTreeIteration = True
    else:
        subTreeIteration = False
    nodeNum = 0
    nodeDict = {}
    nodeList = []
    for node in tree:
        subTree = node.leaf_nodes()
        currLeaves = []
        for leaf in subTree:
            currLeaves.append(leaf.taxon.label)
        nodeDict[nodeNum] = currLeaves
        nodeList.append(nodeNum)
        nodeNum += 1
    if len(outgroupLeaves) > 0:
        speciesComplete = []
        minusOne = []
        for node in nodeList:
            currLeaves = nodeDict[node]
            currIngroup,ingroupSpecies = checkLeaves(currLeaves,ingroupList)
            currOutgroup,outgroupSpecies = checkLeaves(currLeaves,[outgroupID])
            if ingroupSpecies == len(ingroupList) and outgroupSpecies == 0:
                speciesComplete.append(node)
            elif ingroupSpecies == len(ingroupList)-1 and outgroupSpecies == 1:
                minusOne.append(node)
        nonOverlappingNodes = [] 
        while len(nonOverlappingNodes) < len(speciesComplete):
            minLength = False 
            minNode = False
            for node in speciesComplete:
                if node not in nonOverlappingNodes:
                    if minLength == False or len(nodeDict[node]) < minLength:
                        minLength = len(nodeDict[node])
                        minNode = node
            removeList = []
            for node in speciesComplete:
                if node != minNode:
                    for taxon in nodeDict[node]:
                        if taxon in nodeDict[minNode] and node not in removeList:
                            removeList.append(node)
            for item in removeList:
                speciesComplete.remove(item)
            nonOverlappingNodes.append(minNode)
        totalReportingLeaves = []
        for node in nonOverlappingNodes:
            currNodeLeaves = nodeDict[node]
            for leaf in currNodeLeaves:
                totalReportingLeaves.append(leaf)
        nonOverlappingMinusOne = [] 
        while len(nonOverlappingMinusOne) < len(minusOne):
            minLength = False 
            minNode = False
            for node in minusOne:
                if node not in nonOverlappingMinusOne:
                    if minLength == False or len(nodeDict[node]) < minLength:
                        minLength = len(nodeDict[node])
                        minNode = node
            removeList = []
            for node in minusOne:
                if node != minNode:
                    for taxon in nodeDict[node]:
                        if taxon in nodeDict[minNode] and node not in removeList:
                            removeList.append(node)
                for taxon in nodeDict[node]:
                    if taxon in totalReportingLeaves:
                        removeList.append(node)
            for item in removeList:
                minusOne.remove(item)
            nonOverlappingMinusOne.append(minNode)
        for minTree in nonOverlappingNodes: #This section (lines 76-82) writes the summary information of each minimal subtree to the standard out
            ogNum += 1
            topology = monophyletic(nodeDict[minTree],treeFile)
            sys.stdout.write(og + '\t' + str(ogNum) + '\t'+ topology + '\t')
            minimalSubTree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
            minimalSubTree.retain_taxa_with_labels(nodeDict[minTree]+outgroupLeaves,update_bipartitions=True,suppress_unifurcations=False)
            minRetainLeaves = nodeDict[minTree]+outgroupLeaves
            minLeafLength = len(minRetainLeaves)
            minNodeID = 0
            currNodeNum = 0
            for minNode in minimalSubTree:
                minSubTree = minNode.leaf_nodes()
                minCurrLeaves = []
                for leaf in minSubTree:
                    minCurrLeaves.append(leaf.taxon.label)
                minCurrSpecies,minSpeciesNum = checkLeaves(minCurrLeaves,speciesIDList)
                if len(minCurrLeaves) < minLeafLength and minSpeciesNum == len(speciesIDList):
                    minLeafLength = len(minCurrLeaves)
                    minNodeID = currNodeNum
                    minRetainLeaves = minCurrLeaves
                currNodeNum += 1
            minimalSubTree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
            minimalSubTree.retain_taxa_with_labels(minRetainLeaves,update_bipartitions=True,suppress_unifurcations=False)
            minimalSubTreeString = minimalSubTree.as_string(schema='newick')
            currSpeciesLeafDict,currSpeciesNum = checkLeaves(minRetainLeaves,speciesIDList)
            trimLeaves += nodeDict[minTree]
            for currSpecies in speciesIDList:
                sys.stdout.write(str(currSpeciesLeafDict[currSpecies]) + '\t')
            currSeqNum = 0
            for item in currSpeciesLeafDict:
                currSeqNum += currSpeciesLeafDict[item]
            sys.stdout.write(str(len(minRetainLeaves)) + '\t' + str(currSpeciesNum) + '\n')
            treeOutFile = open(og + '.' + str(ogNum) + '.pruned.trees','w') #This section lines(84-89) writes each minimally inclusive subtree to a newick file with a .trees suffix
            treeOutFile.write(minimalSubTreeString)
            treeOutFile.close()
        for minTree in nonOverlappingMinusOne: #This section (lines 76-82) writes the summary information of each minimal subtree to the standard out
            ogNum += 1
            topology = monophyletic(nodeDict[minTree],treeFile)
            sys.stdout.write(og + '\t' + str(ogNum) + '\t'+ topology + '\t')
            minimalSubTree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
            minimalSubTree.retain_taxa_with_labels(nodeDict[minTree],update_bipartitions=True,suppress_unifurcations=False)
            minRetainLeaves = nodeDict[minTree]
            minimalSubTreeString = minimalSubTree.as_string(schema='newick')
            currSpeciesLeafDict,currSpeciesNum = checkLeaves(minRetainLeaves,speciesIDList)
            trimLeaves += nodeDict[minTree]
            for currSpecies in speciesIDList:
                sys.stdout.write(str(currSpeciesLeafDict[currSpecies]) + '\t')
            currSeqNum = 0
            for item in currSpeciesLeafDict:
                currSeqNum += currSpeciesLeafDict[item]
            sys.stdout.write(str(len(minRetainLeaves)) + '\t' + str(currSpeciesNum) + '\n')
            treeOutFile = open(og + '.' + str(ogNum) + '.pruned.trees','w') #This section lines(84-89) writes each minimally inclusive subtree to a newick file with a .trees suffix
            treeOutFile.write(minimalSubTreeString)
            treeOutFile.close()
        try: #This section (lines 96-112) attemps to prune the taxa present in minimally inclusive subtrees out of the main tree. If it fails (e.g., fewer than 3 taxa remain outside the minimally inclusive subtree), the while loop will break sending the program to the final step. If the pruning is successful, the while loop criterion (i.e., all species being present) will be recalculated and reevaluated for the new tree, which is now missing the minimally inclusive subtree
            tree.prune_taxa_with_labels(trimLeaves,update_bipartitions=True,suppress_unifurcations=False)
            treeString = tree.as_string('newick')
            tree = dendropy.Tree.get(data=treeString,schema='newick',preserve_underscores=True)
            totalLeaves = []
            trimLeaves = []
            for taxon in tree.taxon_namespace:
                totalLeaves.append(taxon.label)
            speciesLeafDict,numSpecies = checkLeaves(totalLeaves,speciesIDList)
        except:
            for seq in trimLeaves:
                if seq in totalLeaves:
                    totalLeaves.remove(seq)
            speciesLeafDict,numSpecies = checkLeaves(totalLeaves,speciesIDList)
    else:
        speciesComplete = []
        minusOne = []
        for node in nodeList:
            currLeaves = nodeDict[node]
            currIngroup,ingroupSpecies = checkLeaves(currLeaves,ingroupList)
            currOutgroup,outgroupSpecies = checkLeaves(currLeaves,[outgroupID])
            if ingroupSpecies == len(ingroupList) and outgroupSpecies == 0:
                speciesComplete.append(node)
            elif ingroupSpecies == len(ingroupList)-1 and outgroupSpecies == 1:
                minusOne.append(node)
        nonOverlappingNodes = [] 
        while len(nonOverlappingNodes) < len(speciesComplete):
            minLength = False 
            minNode = False
            for node in speciesComplete:
                if node not in nonOverlappingNodes:
                    if minLength == False or len(nodeDict[node]) < minLength:
                        minLength = len(nodeDict[node])
                        minNode = node
            removeList = []
            for node in speciesComplete:
                if node != minNode:
                    for taxon in nodeDict[node]:
                        if taxon in nodeDict[minNode] and node not in removeList:
                            removeList.append(node)
            for item in removeList:
                speciesComplete.remove(item)
            nonOverlappingNodes.append(minNode)
        for minTree in nonOverlappingNodes: #This section (lines 76-82) writes the summary information of each minimal subtree to the standard out
            ogNum += 1
            topology = monophyletic(nodeDict[minTree],treeFile)
            sys.stdout.write(og + '\t' + str(ogNum) + '\t'+ topology + '\t')
            minimalSubTree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
            minimalSubTree.retain_taxa_with_labels(nodeDict[minTree],update_bipartitions=True,suppress_unifurcations=False)
            minRetainLeaves = nodeDict[minTree]
            minimalSubTreeString = minimalSubTree.as_string(schema='newick')
            currSpeciesLeafDict,currSpeciesNum = checkLeaves(minRetainLeaves,speciesIDList)
            trimLeaves += nodeDict[minTree]
            for currSpecies in speciesIDList:
                sys.stdout.write(str(currSpeciesLeafDict[currSpecies]) + '\t')
            currSeqNum = 0
            for item in currSpeciesLeafDict:
                currSeqNum += currSpeciesLeafDict[item]
            sys.stdout.write(str(len(minRetainLeaves)) + '\t' + str(currSpeciesNum) + '\n')
            treeOutFile = open(og + '.' + str(ogNum) + '.pruned.trees','w') #This section lines(84-89) writes each minimally inclusive subtree to a newick file with a .trees suffix
            treeOutFile.write(minimalSubTreeString)
            treeOutFile.close()
        try: #This section (lines 96-112) attemps to prune the taxa present in minimally inclusive subtrees out of the main tree. If it fails (e.g., fewer than 3 taxa remain outside the minimally inclusive subtree), the while loop will break sending the program to the final step. If the pruning is successful, the while loop criterion (i.e., all species being present) will be recalculated and reevaluated for the new tree, which is now missing the minimally inclusive subtree
            tree.prune_taxa_with_labels(trimLeaves,update_bipartitions=True,suppress_unifurcations=False)
            treeString = tree.as_string('newick')
            tree = dendropy.Tree.get(data=treeString,schema='newick',preserve_underscores=True)
            totalLeaves = []
            trimLeaves = []
            for taxon in tree.taxon_namespace:
                totalLeaves.append(taxon.label)
            speciesLeafDict,numSpecies = checkLeaves(totalLeaves,speciesIDList)
        except:
            for seq in trimLeaves:
                if seq in totalLeaves:
                    totalLeaves.remove(seq)
            speciesLeafDict,numSpecies = checkLeaves(totalLeaves,speciesIDList)
    currSpeciesLeafDict,currSpeciesNum = checkLeaves(totalLeaves,speciesIDList) #This section (lines 115-138) evaluates the remaining tree to see if there are more than two taxa present. If so, the program will write the remaining tree to a newick file (suffix .trees), the sequence information to a corresponding .fasta file, and the summary information to the standard output.
    if len(totalLeaves) > 2:
        ogNum += 1
        topology = monophyletic(totalLeaves,treeFile)
        sys.stdout.write(og + '\t' + str(ogNum) + '\t' + topology + '\t')
        for currSpecies in speciesIDList:
            sys.stdout.write(str(currSpeciesLeafDict[currSpecies]) + '\t')
        sys.stdout.write(str(len(totalLeaves)) + '\t' + str(currSpeciesNum) + '\n')
        remainingTree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
        remainingTree.retain_taxa_with_labels(totalLeaves,update_bipartitions=True,suppress_unifurcations=False)
        remainingTreeString = remainingTree.as_string('newick')
        treeOutFile = open(og + '.' + str(ogNum) + '.pruned.trees','w')
        treeOutFile.write(remainingTreeString)
        treeOutFile.close()
        outfile = open(og + '.' + str(ogNum) + '.pruned.fasta','w')
        '''for taxa in totalLeaves:
            outfile.write('>' + taxa + '\n' + seqDict['>' + taxa] + '\n')
        outfile.close()'''
    elif len(totalLeaves) > 1:
        ogNum += 1
        '''outfile = open(og + '.' + str(ogNum) + '.pruned.fasta','w')
        for taxa in totalLeaves:
            outfile.write('>' + taxa + '\n' + seqDict['>' + taxa] + '\n')
        outfile.close()'''
    
def checkLeaves(leaves,speciesIDList):
    speciesLeafDict = {}
    for speciesID in speciesIDList:
        speciesLeafDict[speciesID] = 0
    for leaf in leaves:
        currSpecies = False
        for speciesID in speciesIDList:
            if speciesID in leaf:
                currSpecies = speciesID
                break
        if currSpecies != False:
            speciesLeafDict[currSpecies] += 1
    numSpecies = 0
    for species in speciesIDList:
        if speciesLeafDict[species] > 0:
            numSpecies += 1
    return speciesLeafDict,numSpecies

def midpointRoot(treeFile):
    tree = dendropy.Tree.get(path = treeFile,schema='newick',preserve_underscores=True)
    tree.reroot_at_midpoint(update_bipartitions=True,suppress_unifurcations=False)
    rootedTree = tree.as_string(schema='newick')
    return rootedTree

def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            seqName = seqName.replace('_',' ')
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            seqName = seqName.replace('_',' ')
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict

def monophyletic(listOfTaxa,treeFile):
    bigTree = dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
    nodeDict = {}
    nodeList = []
    nodeNum = 0
    totalLeaves = []
    listOfTaxa.sort()
    for taxon in bigTree.taxon_namespace:
        totalLeaves.append(taxon.label)
    for node in bigTree: 
        subTree = node.leaf_nodes()
        currLeaves = []
        for leaf in subTree:
            currLeaves.append(leaf.taxon.label)
        currLeaves.sort()
        nodeDict[nodeNum] = currLeaves
        nodeList.append(nodeNum)
        nodeNum += 1
        if listOfTaxa == currLeaves:
            return 'Monophyletic'
    return 'Paraphyletic'

infile = open(sys.argv[1])
sys.stdout.write('OGID\tsubTreeNumber\tMonophyletic?\t')
for currSpecies in sys.argv[2:]:
    sys.stdout.write(currSpecies + '\t')
sys.stdout.write('numSequences\tnumSpecies\n')
for line in infile:
    realLine = line
    while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
        realLine = realLine[0:-1]
    subTreeIterator(realLine,sys.argv[2:])

import sys
import re
import logging
import igraph as ig

logger = logging.getLogger(__name__)

def read_node_file(filename):
    nodes = {}

    with open(filename) as fh:
        for line in fh.readlines():
            fields = line.rstrip().split()
            nodes[fields[0]] = 1

    return nodes

def read_graph_from_gfa(gfafile, args):

    allnodes = {}
    nodelines = {}
    edgelines = []
    nodelengths = {}
    with open(gfafile, "r") as gfh:
        for line in gfh.readlines():
            fields = line.rstrip().split()
            if fields[0] == 'S':
                nodename = fields[1]
                lengthmatch = re.search('(?<=LN:i:)\d+$', fields[3])
                if lengthmatch:
                    nodelengths[nodename] = int(lengthmatch.group(0))
                nodelines[nodename] = line
            elif fields[0] == 'L':
                edgelines.append(line)
                allnodes[fields[1]] = 1
                allnodes[fields[3]] = 1

    #construct graph:
    n_vertices = len(allnodes.keys())
    node_indices = {}
    nodes = []
    edges = []
    # directions will be ">" for "-" nodes and "<" for "+" nodes
    directionlist = []
    index = 0
    for node in allnodes.keys():
        node_indices[node] = index
        nodes.append(node)

        #print("Node " + node + " has length " + str(parednodelengths[node]))
        index = index + 1
    for line in edgelines:
        fields = line.rstrip().split()
        node1 = fields[1]
        dir1 = fields[2]
        node2 = fields[3]
        dir2 = fields[4]
        edges.append([node_indices[node1], node_indices[node2]])
        directionlist.append(dir1 + dir2)

    assemblygraph = ig.Graph(n=n_vertices, edges=edges, directed=True)
    assemblygraph.vs["name"] = nodes
    assemblygraph.es["direction"] = directionlist

    return assemblygraph

def find_hom_het_nodes(graph, args):

    # look for nodes with two incoming or two outgoing edges:
    nodelist = graph.vs
    nodenamelist = nodelist["name"]
    edgelist = graph.es
    for node in graph.vs:
        attribdict = node.attributes()
        nodename = attribdict["name"]
        #print(nodename)
        nodeedges = node.out_edges()
        numedges = len(nodeedges)
        #print("Node " + nodename + " has " + str(numedges) + " outgoing edges:")
        numbackedges = 0
        numforwardedges = 0
        backedges = []
        forwardedges = []
        for edge in nodeedges:
            node1 = nodenamelist[edge.tuple[0]]
            node2 = nodenamelist[edge.tuple[1]]
            node1dir = edge.attributes()["direction"][0]
            node2dir = edge.attributes()["direction"][1]
            if node1dir == "-":
                numbackedges = numbackedges + 1
                backedges.append(edge)
            elif node1dir == "+":
                numforwardedges = numforwardedges + 1
                forwardedges.append(edge)
            #print(str(edge.tuple) + " " +  node1 + " " + node2 + " " + edge.attributes()["direction"])
        if numbackedges == 2:
            # trace back to a hom node if possible:
            # make a list of nodes extending from the 5' end of the unitig, as well as nodes at ends of edges at ends of single edges from the opposite ends:
            backhetnodes = []
            backhetsecondnodes = []
            for backedge in backedges:
                node2dir = backedge.attributes()["direction"][1]
                backhetnodes.append(nodenamelist[backedge.tuple[1]])
                doublebackhetedges = nodelist[backedge.tuple[1]].out_edges()
                extendingdoublebackhetedges = []
                for hetedge in doublebackhetedges:
                    if hetedge.attributes()["direction"][0] == node2dir:
                        extendingdoublebackhetedges.append(hetedge)
                numdoublebackhetedges = len(extendingdoublebackhetedges)
                if numdoublebackhetedges == 1:
                    backhetsecondnodes.append(nodenamelist[extendingdoublebackhetedges[0].tuple[1]])
            if len(backhetsecondnodes)==2 and backhetsecondnodes[0]==backhetsecondnodes[1]:
                print("Found backward hets from " + node1 + " named " + backhetnodes[0] + " and " + backhetnodes[1] + " leading to " + backhetsecondnodes[0])
        if numforwardedges == 2:
            # trace fwd to a hom node if possible:
            # make a list of nodes extending from the 3' end of the unitig, as well as nodes at ends of edges at ends of single edges from the opposite ends:
            fwdhetnodes = []
            fwdhetsecondnodes = []
            for fwdedge in forwardedges:
                node2dir = fwdedge.attributes()["direction"][1]
                fwdhetnodes.append(nodenamelist[fwdedge.tuple[1]])
                doublefwdhetedges = nodelist[fwdedge.tuple[1]].out_edges()
                extendingdoublefwdhetedges = []
                for hetedge in doublefwdhetedges:
                    if hetedge.attributes()["direction"][0] == node2dir:
                        extendingdoublefwdhetedges.append(hetedge)
                numdoublefwdhetedges = len(extendingdoublefwdhetedges)
                if numdoublefwdhetedges == 1:
                    fwdhetsecondnodes.append(nodenamelist[extendingdoublefwdhetedges[0].tuple[1]])
            if len(fwdhetsecondnodes)==2 and fwdhetsecondnodes[0]==fwdhetsecondnodes[1]:
                print("Found forward hets from " + node1 + " named " + fwdhetnodes[0] + " and " + fwdhetnodes[1] + " leading to " + fwdhetsecondnodes[0])

def compare_compressed_fasta_to_node_seqs(fastafile:str, nodegfa:str, args):
    minalignlength = args.minalignlength
    if args.prefix is not None:
        assemblyprefix = args.prefix
        nodefasta = args.prefix + ".nodes.fasta"
        compressedfasta = args.prefix + ".hpc.fasta"
        chainfile = args.prefix + ".hpc.chain"
    else:
        compressedfasta = re.sub("\.fa.*$", ".hpc.fasta", fastafile)
        compressedfasta = re.sub(".*/", "", compressedfasta)
        assemblyprefix = re.sub("\.fa.*$", "", fastafile)
        assemblyprefix = re.sub(".*/", "", assemblyprefix)
        nodefasta = re.sub("\.gfa.*$", ".nodes.fasta", nodegfa)
        nodefasta = re.sub(".*/", "", nodefasta)
        chainfile = re.sub("\.fa.*$", ".hpc.chain", fastafile)
        chainfile = re.sub(".*/", "", chainfile)

    paffile = re.sub("\.fasta$", "_vs_", nodefasta)
    paffile = paffile + assemblyprefix + ".paf"
    hpcbed = re.sub(".paf$", ".sort.bed", paffile)
    decompbed = re.sub(".bed", ".decompressed.bed", hpcbed)
    decompunmapped = re.sub(".bed", ".decompressed.unmapped", hpcbed)

    if os.path.isfile(nodefasta):
        print("Using pre-existing node fasta file " + nodefasta)
    else:
        print("Writing node fasta file " + nodefasta)
        write_gfa_fasta(nodegfa, nodefasta)

    if os.path.isfile(paffile):
        print("Using pre-existing alignment file " + paffile)
    else:
        print("Running minimap2 to generate " + paffile)
        commandwords = ['minimap2', '-k', '28', '-w', '100', '-x', 'lr:hq', '-t', '2', '-r', '10,10', '-N', '1', nodefasta, compressedfasta, '-o', paffile]
        subprocess.run(commandwords)

    bedlines = {} # dictionary to be sorted before printing
    with open(paffile, "r") as pfh:
        for line in pfh:
            fields = line.split()
            compchrom = fields[0]
            compchromlength = int(fields[1])
            compstart = int(fields[2])
            compend = int(fields[3])
            strand = fields[4]
            node = fields[5]
            if strand == "+":
                node = node + "<"
            else:
                node = node + ">"
            nodelength = int(fields[6])
            nodealigned = int(fields[8]) - int(fields[7])
            fracaligned = int(100*nodealigned/nodelength)
            if compend - compstart < minalignlength:
                continue
            compstartstring = str(compstart).zfill(12)
            compendstring = str(compend).zfill(12)
            bedlines[compchrom + ":" + compstartstring + ":" + compendstring] = compchrom + "\t" + str(compstart) + "\t" + str(compend) + "\t" + node + "\t" + str(fracaligned) + "\n"

    with open(hpcbed, "w") as bfh:
        for posstring in sorted(bedlines.keys()):
            bfh.write(bedlines[posstring])

    print("Running liftover to decompress coords in " + hpcbed)
    commandwords = ['liftOver', hpcbed, chainfile, decompbed, decompunmapped]
    subprocess.run(commandwords)

def write_gfa_fasta(nodegfa:str, nodefasta:str):

    with open(nodefasta, "w") as nfh:
        with open(nodegfa, "r") as gfh:
            for line in gfh:
                fields = line.split()
                #print(fields[0])
                if fields[0]=="S":
                    nodename = fields[1]
                    nodeseq = fields[2]
                    nfh.write(">" + nodename + "\n" + nodeseq + "\n")

    return(0)


import toyplot

def create_sequences(f):
    # Open file that was outputted by Lighter
    file = open(f, 'r')
    data = file.readlines()

    # Create a list to hold all of the sequences from the file
    sequences = []

    # Iterate through data, adding the corrected sequences to the list
    for i in range(1,len(data),4):
        sequences.append(data[i][:-1])

    return sequences



def create_kmers(sequences, k):
    # Create a dictionary to hold the kmers generated
    kmers = {}

    # Iterate through the list of sequences
    for sequence in sequences:
        
        # Create kmers from sequence
        for i in range(len(sequence)):
            kmer = sequence[i:i+k]

            #Check to see if we need to loop around to the front of the sequence
            if len(kmer) != k:
                kmer += sequence[:(k-len(kmer))]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1

    # Return the list of kmers
    return kmers

def create_debruijn(kmers):
    # Create a set to hold the edge values:
    edges = set()

    # Iterate through all pairs of kmers
    for k1mer in kmers.keys():
        for k2mer in kmers.keys():

            # Check to make sure kmers are distinct
            if k1mer != k2mer:

                # If kmers are distinct and the first k-1 values of k1mer are equal to the last k-1 values of k2mer, add an edge to the set connecting the first k-1 values of both.
                if k1mer[:-1] == k2mer[1:]:
                    edges.add((k2mer[:-1], k1mer[:-1]))
                if k1mer[1:] == k2mer[:-1]:
                    edges.add((k1mer[:-1], k2mer[:-1]))
    return edges

# Function to plot debruijn_graph
# Plotting function taken from Eaton Lab: https://eaton-lab.org/slides/genomics/answers/nb-10.2-de-Bruijn.html
# Originally derived from https://www.nature.com/articles/nbt.2023
# Citations will be provided in README

def plot_debruijn_graph(edges, width=1500, height=1500):
    # returns a toyplot graph from an input of edges
    graph = toyplot.graph(
        [i[0] for i in edges],
        [i[1] for i in edges],
        width=width,
        height=height,
        tmarker=">", 
        vsize=25,
        vstyle={"stroke": "black", "stroke-width": 2, "fill": "none"},
        vlstyle={"font-size": "11px"},
        estyle={"stroke": "black", "stroke-width": 2},
        layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
    return graph

#### To create edges for usage in part 5 if necessary.
# seq = create_sequences('../data/output.cor.fq')
# k = create_kmers(seq, 4)
# edges = create_debruijn(k)

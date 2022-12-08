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


# Different, recursive approach, Unfinished, not yet working
def visit(kmers, edge_list, current, contig_list, contig):
    next_found = False

    if current[3] == 0:  # If the current edge has no remaining traversals...

        # End remaining contig
        contig_list.append(contig)
        contig = []

        for edge in edge_list:
            if edge[2] > 0:
                edge[2] = edge[2] - 1
                visit(edge)

    else:  # Else look for an edge whose start k-1mer is equal to the current's end k-1mer
        for edge in edge_list:
            if current[1] == edge[0] and current[2] > 0 and edge[2] > 0:  # Also check for remaining reads

                next_found = True

                current = edge
                # This will ultimately construct the sequence backwards, we'll reverse this later
                contig.append(current)

                index = edge_list.index(current)
                edge_list[index][2] = edge_list[index][2] - 1  # Decrement remaining edge traversals

                visit(kmers, edge_list, current, contig_list, contig)  # Make next recursive call

    if not next_found:
        contig_list.append(contig)
        contig = []

        current = edge_list[0]
        index = edge_list.index(current)
        contig.append(current)
        edge_list[index][2] = edge_list[index][2] - 1

        visit(kmers, edge_list, current, contig_list, contig)



def create_circuit(kmers, edges):
    contig = []
    contig_list = []  # If one continuous sequence cannot be recreated, a list of contigs is instead created
    edge_list = []

    # Create an edge list that contains [start_node, end_node, times read]
    for edge in edges:
        count = kmers.get(edge[0] + edge[1][-1])
        edge_list.append([edge[0], edge[1], count])

    # Make initial call to recursive visit function
    visit(kmers, edge_list, edge_list[0], contig_list, contig)

    contig_list.reverse()

    return contig_list


## Finished, but doesn't work correctly.
# def create_circuit(kmers, edges):
#     edge_list = []
#     contig = []
#     contig_list = []
#     for edge in edges:  # Create an edge list that contains start_node, end_node, times_travelled
#         count = kmers.get(edge[0] + edge[1][-1])
#         edge_list.append([edge[0], edge[1], count])
#
#     print(edge_list)
#
#     current = edge_list[0]
#     index = edge_list.index(current)
#     contig.append(current)
#     edge_list[index][2] = edge_list[index][2] - 1
#
#     i = 0
#
#     while len(edge_list) > 0:
#
#         found_match = False
#
#         print(edge_list)
#
#         i = i + 1
#
#         for edge in edge_list:
#             if current[1] == edge[0] and current[2] > 0 and edge[2] > 0:
#                 current = edge
#                 index = edge_list.index(current)
#                 contig.append(current)
#                 edge_list[index][2] = edge_list[index][2] - 1
#                 found_match = True
#                 break
#
#         if not found_match:
#             contig_list.append(contig)
#             contig = []
#
#             current = edge_list[0]
#             index = edge_list.index(current)
#             contig.append(current)
#             edge_list[index][2] = edge_list[index][2] - 1
#
#     return contig_list


## Creates a circuit, Correct and working, but does not handle duplicates
# def create_circuit(edges):
#     edge_list = []
#     contig = []
#     contig_list = []
#     for edge in edges:
#         edge_list.append(edge)
#
#     current = edge_list[0]
#     contig.append(current)
#     edge_list.pop(0)
#
#     while edge_list != []:
#
#         found_match = False
#
#         for edge in edge_list:
#             if current[1] == edge[0]:
#                 current = edge
#                 contig.append(current)
#                 edge_list.remove(edge)
#                 found_match = True
#                 break
#
#         if not found_match:
#             contig_list.append(contig)
#             contig = []
#
#             current = edge_list[0]
#             contig.append(current)
#             edge_list.pop(0)
#
#     return contig_list


seq = create_sequences('../data/output.cor.fq')
k = 4
kmers = create_kmers(seq, k)
edges = create_debruijn(kmers)

contigs = create_circuit(kmers, edges)

# create_circuit(kmers)

print(kmers)

print(contigs)

# for contig in contigs:
#     print("")
#     for x, y in contig:
#         print(x, end='')

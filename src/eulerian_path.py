import toyplot


def create_sequences(f):
    # Open file that was outputted by Lighter
    file = open(f, 'r')
    data = file.readlines()

    # Create a list to hold all of the sequences from the file
    sequences = []

    # Iterate through data, adding the corrected sequences to the list
    for i in range(1, len(data), 4):
        sequences.append(data[i][:-1])

    return sequences


def create_kmers(sequences, k):
    # Create a dictionary to hold the kmers generated
    kmers = {}

    # Iterate through the list of sequences
    for sequence in sequences:

        # Create kmers from sequence
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]

            # Count how many times we see each kmer
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1

    # Return the list of kmers
    return kmers


# def create_debruijn(kmers):
#     # Create a set to hold the edge values:
#     edges = set()
#
#     # Iterate through all pairs of kmers
#     for k1mer in kmers.keys():
#         for k2mer in kmers.keys():
#
#             # Check to make sure kmers are distinct
#             if k1mer != k2mer:
#
#                 # If kmers are distinct and the first k-1 values of k1mer are equal to the last k-1 values of k2mer, add an edge to the set connecting the first k-1 values of both.
#                 if k1mer[:-1] == k2mer[1:]:
#                     edges.add((k2mer[:-1], k1mer[:-1]))
#
#                 if k1mer[1:] == k2mer[:-1]:
#                     edges.add((k1mer[:-1], k2mer[:-1]))
#     return edges


def create_debruijn(kmers):
    # Create a set to hold the edge values:
    edges = set()

    # Iterate through all pairs of kmers, split them into k-1mers
    for kmer in kmers.keys():
        edges.add((kmer[:-1], kmer[1:]))
    return edges


def create_debruijn_with_multiplicity(kmers):
    # Create a set to hold the edge values:
    edges = set()

    # Iterate through all pairs of kmers, split them into k-1mers, and record the edge's multiplicity
    for kmer in kmers.items():
        edges.add((kmer[0][:-1], kmer[0][1:], kmer[1]))
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

# Creates a path from edges
def create_path(edges):
    edge_list = []
    contig = []
    contig_list = []
    count = 0

    for edge in edges:  # Create an edge list that we can change and loop through
        new_edge = [edge[0], edge[1], edge[2]]
        edge_list.append(new_edge)
        count += edge[2]  # Keep track of the total edge count so that we can stop once all have been travelled
    # max_count = count

    current = edge_list[0]  # Set starting edge

    current[2] -= 1  # Decrement starting edge's remaining traversals
    count -= 1  # Decrement total remaining edge traversals

    contig.append(current[:-1])

    while count > 0:

        # print("Creating Contigs:", '{:.2f}%'.format((max_count - count) / max_count * 100), end='\r')

        found_match = False

        for edge in edge_list:

            # If first k-1mer of current equals the 2nd k-1mer of edge and edge has remaining traversals left...
            if current[1] == edge[0] and edge[2] > 0:
                current = edge  # Step to new edge
                current[2] -= 1  # Decrement remaining traversals
                count -= 1  # Decrement total remaining edge traversals

                contig.append(current[:-1])  # Add edge to contig

                found_match = True
                break

        # If we did not find a match, then we need to start a new contig.
        if not found_match:
            contig_list.append(contig)  # Add current contig to the contig list
            contig = []  # Reset contig

            for edge in edge_list:  # Find new starting edge that has traversals remaining
                if edge[2] > 0:
                    current = edge

                    current[2] -= 1
                    count -= 1

                    contig.append(current[:-1])  # Add starting edge

                    break

    # print("Creating Contigs: 100.00%")

    contig_list.append(contig)  # Add final contig to contig list

    # Combine as many contigs as possible...
    # print("\nConnecting Contigs...")

    assembled_contig = True

    while assembled_contig:

        # print(f'Contigs Remaining: {len(contig_list)}', end='\r')
        assembled_contig = False

        for contig1 in contig_list:
            for contig2 in contig_list:
                if contig1 != contig2 and contig1[len(contig1) - 1][1] == contig2[0][0]:  # If contigs can be joined...
                    for edge in contig2:  # Append contig2 to contig1
                        contig1.append(edge)
                    contig_list.remove(contig2)  # Remove contig2

                    assembled_contig = True
                    break

            # Break out of second for loop
            if assembled_contig:
                break

    # print(f'Contigs Remaining: {len(contig_list)}')

    return contig_list


def find_eulerian_path(sequence_count, k, output_path, label):

    sequences = create_sequences('../data/output.cor.fq')
    kmers = create_kmers(sequences[:sequence_count], k)
    edges = create_debruijn_with_multiplicity(kmers)
    contigs = create_path(edges)

    # Find longest contig
    contigs.sort(key=len, reverse=True)
    longest_contig = contigs[0]

    # print(f'Length of Longest Contig: {len(longest_contig)}')

    # Write longest contig to contig.fna
    f = open(output_path, "w")
    f.write(f'>{label}\n')
    f.write(longest_contig[0][0])  # Write first node
    f.write(longest_contig[0][1][-1])  # Write last character of 2nd node

    for i in range(1, len(longest_contig)):  # Write the rest of the characters
        f.write(longest_contig[i][1][-1])

    f.write("\n")
    f.close()

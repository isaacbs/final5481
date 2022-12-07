# final5481
Final Project regarding shotgun sequencing of the SARS-CoV2 spike protein

## 1.) Using Trimmomatic to clean data:
Used the following command to clean the data: 
"java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 sars_spike_protein_raw_reads.fastq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:5:15:12 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:50"

# Output: 
Input Reads: 5000 Surviving: 4828 (96.56%) Dropped: 172 (3.44%)
TrimmomaticSE: Completed successfully

This output can be seen in data/output.fq
# Parameters:
The options perform the following:
- Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:5:15:12)
- Remove leading low quality or N bases (below quality 3) (LEADING:3)
- Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
- Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 25 (SLIDINGWINDOW:4:25)
- Drop reads below the 50 bases long (MINLEN:50)


## 2.) Using Lighter to correct data:
The following command was used and the following output was created.
./lighter -r ../University\ of\ Minnesota/University-of-Minnesota-Projects/Computational\ Techniques\ for\ Genomics/final5481/data/output.fq -k 17 4000 0.03 -t 8

# Output:
[2022-12-06 13:00:36] =============Start====================
[2022-12-06 13:00:36] Bad quality threshold is "E"
[2022-12-06 13:00:37] Finish sampling kmers
[2022-12-06 13:00:37] Bloom filter A's false positive rate: 0.002184
[2022-12-06 13:00:38] Finish storing trusted kmers
[2022-12-06 13:00:39] Finish error correction
Processed 4828 reads:
        3794 are error-free
        Corrected 1170 bases(1.131528 corrections for reads with errors)
        Trimmed 0 reads with average trimmed bases 0.000000
        Discard 0 reads

This output can be seen in data/output.cor.fq
# Parameters
    -r seq_file: seq_file is the path to the sequence file
    -k kmer_length genome_size alpha
    -t number of threads
These parameters were chosen because the sars_spike_protein has a length of 3822 which was rounded up to 4000. I chose a K-mer value of 17, and in order to calculate our alpha value I found our total number of reads to be 4828 after trimming using Trimmomatic. The length of the longest read was 151 bp so I used (151*4828)/4000 to find our average coverage. Which turned out to be ~180. From there I calculated the alpha using the formula:
alpha = 7/C, where C is the average coverage. This yielded a value of 0.03. For t, which represents the number of threads, I used a value of 8 just because I knew I had at least that many threads on my computer.


## 3 and 4.) Break the corrected sequencing errors into a list of k-mers for some small k and Create a DeBruijn Graph
The code for these two portions of the project can be accessed either through the jupyter notebook called "Parts 3 and 4.ipynb"
or through the file p3.py.

The notebook has more detailed representations and shows the graphs for the different levels of k and different percentages of the sequence data used.

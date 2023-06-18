#!/usr/bin/env python3

import argparse
import pandas as pd

# imports vsearch clustering file 
# outputs a three column TSV file listing, for unique sequence, the name, count, and relative abundance
# outputs a text file listing the name of each unique sequence

def make_feature_table(input_file, output_basename, min_count, min_abundance):
    df = pd.read_csv(args.input_file, sep = "\t", header = None)

    clusters = df[df[0] == "C"][[8, 2]]
    clusters.columns = ["Cluster", "Count"]

    clusters = clusters[clusters["Count"] >= args.min_count]
    clusters["RelAbund"] = 100 * clusters["Count"] / clusters["Count"].sum()
    clusters = clusters[clusters["RelAbund"] >= args.min_abundance]

    clustersoutput = args.output_basename + '_clusters.txt'
    clusters.to_csv(clustersoutput, sep="\t", index = False, header = False)

    readnameoutput = args.output_basename + '_readnames.txt'
    with open(readnameoutput, 'w') as f:
        for text in clusters['Cluster'].tolist():
            f.write(text + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Make a feature table from vsearch culstering output")
    parser.add_argument("-i", "--input_file", dest = 'input_file', required = True, help = "Input file")
    parser.add_argument("-o", "--output_basename", dest = 'output_basename', required = True, help = "Basebane (prefix) of output files")
    parser.add_argument("-m", "--min_count", dest = 'min_count', type = int, default = 2, help = "Minimum time each unique sequence must be detected to be included in relative abundance estimation. Sequences detected less frequently are discarded. Default value is 2 to discard singletons.")
    parser.add_argument("-a", "--min_abundance", dest = 'min_abundance', type = float, default = 0.1, help = "Minimum relative abundance of a sequence to be reported in feature table. Default is 0.1")
    args = parser.parse_args()

    make_feature_table(args.input_file, args.output_basename, args.min_count, args.min_abundance)

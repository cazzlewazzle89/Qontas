import argparse
import pandas as pd

def make_feature_table(input_file, output_file, min_count, min_abundance):
    df = pd.read_csv(args.input_file, sep = "\t", header = None)

    clusters = df[df[0] == "C"][["V9", "V3"]]
    clusters.columns = ["Cluster", "Count"]

    clusters = clusters[clusters["Count"] >= args.min_count]
    clusters["RelAbund"] = 100 * clusters["Count"] / clusters["Count"].sum()
    clusters = clusters[clusters["RelAbund"] >= args.min_abundance]

    clusters.to_csv(args.output_file, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", dest = 'input_file', required = True, help = "Input file")
    parser.add_argument("-o", "--output_file", dest = 'output_file', required = True, help = "Output file")
    parser.add_argument("-m", "--min_count", dest = 'min_count', type = int, default = 2, help = "Minimum time each unique sequence must be detected to be included in relative abundance estimation. Sequences detected less frequently are discarded. Default value is 2 to discard singletons.")
    parser.add_argument("-a", "--min_abundance", dest = 'min_abundance', type = float, default = 0.1, help = "Minimum relative abundance (%) of a sequence to be reported in feature table. Default is 0.1")
    args = parser.parse_args()
    make_feature_table(args.input_file, args.output_file, args.min_count, args.min_abundance)

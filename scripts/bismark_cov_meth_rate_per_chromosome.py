#' ---
#' title:  bismark_cov_meth_rate_per_chromosome.py
#' author: Sebastian Mueller (sebm_at_posteo.de)
#' date:   2019-03-28
#' details: script that takes bismark.cov files and calculates % methylation
#' for each chromosome (stdout). In particalar useful to assess conversion rates, e.g. for
#' plants add chloroplast or lambda phage to reference and calculate conversion rate for it
#' using this script.
#' ---

import argparse
import csv
parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input file (*.bismark.cov)")
# parser.add_argument("output", help="Output file (csv)")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
args = parser.parse_args()
# print(args.input)
# print(args.output)

methall, unmethall = 0, 0

# create dictonary with chromosomes as keys and added up meth counts as value
count_dict_meth = {}
count_dict_unmeth = {}

with open(args.input, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
    for row in spamreader:
        # print(', '.join(row))
        chrom, meth, unmeth = row[0], row[4], row[5]
        count_dict_meth[chrom]   = count_dict_meth.get(chrom, 0)   + int(meth)
        count_dict_unmeth[chrom] = count_dict_unmeth.get(chrom, 0) + int(unmeth)

print("Chromosome\tMethylated\tUnmethylated\t%Methylated")
for key in count_dict_meth:
    methall   = count_dict_meth[key]
    unmethall = count_dict_unmeth[key]
    print("{}\t{}\t{}\t{:.2f}".format(key, methall, unmethall, methall / (methall+unmethall) * 100))

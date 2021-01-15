# -*- coding: utf-8 -*-


# The input are the splice junction files from STAR will look like such:
# column 1:  chromosome
# column 2:  first base of the intron (1-based)
# column 3:  last base of the intron (1-based)
# column 4:  strand (0:  undefined, 1:  +, 2:  -)
# column 5:  intron  motif:  0:  non-canonical;  1:  GT/AG,  2:  CT/AC,  3:  GC/AG,  4:  CT/GC,  5:AT/AC, 6:  GT/AT
# column 6:  0:  unannotated, 1:  annotated (only if splice junctions database is used)
# column 7:  number of uniquely mapping reads crossing the junction
# column 8:  number of multi-mapping reads crossing the junction
# column 9:  maximum spliced alignment overhang
# input looks like.....
# 0      1        2     3   4   5  6   7     8
# chr1    11672    12009    1    1    1    0    2    67
# chr1    12228    12612    1    1    1    0    1    31
# chr1    14830    14969    2    2    1    75    162    71
#############################################
# The output will look like such so that score is the number of unique mappers
# and the name is the original 1-based star junction and if it was annotated
#############################################
# chr1    11671    12008    chr1:11672-12009|1    0    +
# chr1    12227    12611    chr1:12228-12612|1    0    +
# chr1    14829    14968    chr1:14830-14969|1    75    -

import argparse, sys, os.path
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input", help="File with the regions in bed format")
    parser.add_argument("-o","--output", help="Name of the output bed you want")
    parser.add_argument("-m","--motifFilter", help=" filter out splice junctions with \n\t\t non-canonical motifs",action="store_true")
    parser.add_argument("-n","--name", help="File with the regions in bed format",action="store_true")


    args = parser.parse_args()

    infile = args.input
    outfile = args.output
    motifON = args.motifFilter
    nameBase = args.name


    if infile is not None and outfile is not None:
        run(infile, outfile, motifON,nameBase)
    else:
        usage()


def usage():
    print("\nUsage: python splicejunction2bed [options] <mandatory>")
    print("Used for turning STAR's SJ.out.tab into bed files")
    print("Options:")
    print("\t-h, --help:\n\t\t show this help message and exit")
    print("\t-n, --name:\n\t\t change name column to the name of the file")

    print("Mandatory:")
    print("\t-i, --input:\n\t\t File with the regions in bed format")
    print("\t-o, --output:\n\t\t Name of the gtf file output file. Directory where the file will be created should exist!")


def run(infile, outfile, motifON, nameBase):

    inf  = open(infile, 'r')
    outf = open(outfile,'w')


    for line in inf:
        linea_split = line.split(',')

        chrom = linea_split[0]
        # if the intromotif filter is on
        # column 5:  intron  motif:  0:  non-canonical;  1:  GT/AG,  2:  CT/AC,  3:  GC/AG,  4:  CT/GC,  5:AT/AC, 6:  GT/AT
        if (str(linea_split[0]) == "seqnames"):
            name = str(linea_split[8]).rstrip()
            continue
        else:
            # convert from one based to zero based with -1
            ini_pos = str(int(linea_split[1]) - 1)
            fin_pos = str(int(linea_split[2]) + 1)

            strand = linea_split[3]
            # score is the PSI
            score = str(linea_split[8]).rstrip()
            annotation = str(linea_split[4]).rstrip()
            score = score + "_" + annotation
            # chr1    11671    12008    chr1_11672_12009:1    0    +
            # chr1    12227    12611    chr1_12228_12612:1    0    +
            # chr1    14829    14968    chr1_14830_14969:1    75    -
            write_line = [chrom, ini_pos, fin_pos, name, score, strand]
            write_line = "\t".join(write_line)

            outf.write(write_line + "\n")



    inf.close()
    outf.close()


if __name__ == "__main__":
    main()

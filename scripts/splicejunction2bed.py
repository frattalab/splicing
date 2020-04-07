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
# chr1	11672	12009	1	1	1	0	2	67
# chr1	12228	12612	1	1	1	0	1	31
# chr1	14830	14969	2	2	1	75	162	71
#############################################
# The output will look like such so that score is the number of unique mappers
# and the name is the original 1-based star junction and if it was annotated
#############################################
# chr1	11671	12008	chr1:11672-12009|1	0	+
# chr1	12227	12611	chr1:12228-12612|1	0	+
# chr1	14829	14968	chr1:14830-14969|1	75	-

import getopt, sys, os.path


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:m:", ["help", "input=","output=","motifFilter="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None
    outfile = None
    motifON = False

    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            if os.path.isfile(arg):
                infile = arg
        elif opt in ("-o", "--output"):
            outfile = arg
        elif opt in ("-m", "--motifFilter"):
            motifON = True
        else:
            print(opt)
            assert False, "Unhandled option"

    if infile is not None and outfile is not None:
        run(infile, outfile, motifON)
    else:
        usage()


def usage():
    print("\nUsage: python splicejunction2bed [options] <mandatory>")
    print("Used for turning STAR's SJ.out.tab into bed files")
    print("Options:")
    print("\t-h, --help:\n\t\t show this help message and exit")
    print("\t-mf, --motifFilter:\n\t\t filter out splice junctions with \n\t\t non-canonical motifs")
    print("Mandatory:")
    print("\t-i, --input:\n\t\t File with the regions in bed format")
    print("\t-o, --output:\n\t\t Name of the gtf file output file. Directory where the file will be created should exist!")


def run(infile, outfile, motifON):

    inf  = open(infile, 'r')
    outf = open(outfile,'w')

    for line in inf:
        linea_split = line.split()

        chrom = linea_split[0]
        # if the intromotif filter is on
        # column 5:  intron  motif:  0:  non-canonical;  1:  GT/AG,  2:  CT/AC,  3:  GC/AG,  4:  CT/GC,  5:AT/AC, 6:  GT/AT
        if (motifON) & (str(linea_split[4]) == "0"):
            continue
        else:
            # convert from one based to zero based with -1
            ini_pos = str(int(linea_split[1]) - 1)
            fin_pos = str(int(linea_split[2]) + 1)
            # strand needs to be converted to -/+/*
            strand = linea_split[3]
            # # column 4:  strand (0:  undefined, 1:  +, 2:  -)
            if(strand == "1"):
                strand = "+"
            elif(strand == "2"):
                strand = "-"
            elif(strand == "0"):
                strand = "*"

            # name here will be the original coords separated by underscore
            name = str(linea_split[0]) + ":" + str(linea_split[1]) + "-" + str(linea_split[2])

            # and then : and the annotation
            name = name + "|" + str(linea_split[5])
            # score is the number of uniquemappers
            score = str(linea_split[6])
            # chr1	11671	12008	chr1_11672_12009:1	0	+
            # chr1	12227	12611	chr1_12228_12612:1	0	+
            # chr1	14829	14968	chr1_14830_14969:1	75	-
            write_line = [chrom, ini_pos, fin_pos, name, score, strand]
            write_line = "\t".join(write_line)

            outf.write(write_line + "\n")



    inf.close()
    outf.close()


if __name__ == "__main__":
    main()

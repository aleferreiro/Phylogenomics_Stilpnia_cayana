#!/usr/bin/env python3.10

import numpy as np


def make(in_loci, out_vcf, names_file, full_vcf=True):
    mindepth = 10 # A fixed value
    version = 0.1 # A fixed value
    outfile  =  open(out_vcf, 'w')
    loci = open(in_loci).read().split("|")[:-1]

    with open(names_file) as n:
        names = n.read().split()
    names.sort()

    outfile.write(VCF_HEADER.format(version))
    outfile.write("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO    ","FORMAT"]+list(names)))
    outfile.write("\n")

    vcflist = []
    for locusnumber in range(len(loci)):
#        import pdb; pdb.set_trace()
        samps = [i.split()[0] for i in loci[locusnumber].strip().split("\n")[:-1]]
        loc = np.array([tuple(i.split()[-1]) for i in loci[locusnumber].strip().split("\n")[:-1]])
        NS = str(len(loc))
        DP = str(mindepth)
        for base in range(len(loc.T)):
            col = []
            site = list(loc.T[base])
            site = list("".join(site).replace("-","").replace("N",""))
            if site:
                for bb in site:
                    if bb in list("RKYSWM"):
                        col += unstruct(bb)[0]
                        col += unstruct(bb)[1]
                    else:
                        col += bb
                REF = most_common([i for i in col if i not in list("-RKYSWMN")])
                ALT = set([i for i in col if (i in list("ATGC-N")) and (i!=REF)])
                if not ALT:
                    if not full_vcf:
                        # If only writing snps then skip this site
                        continue
                    else:
                        ALT = "."
                     
                GENO = [REF]+list(ALT)
                GENOS = []
                for samp in names:
                    if samp in samps:
                        idx = samps.index(samp)
                        f = unstruct(loc.T[base][idx])
                        if ('-' in f) or ('N' in f):
                            GENOS.append("./.")
                        else:
                            GENOS.append(str(GENO.index(f[0]))+"/"+str(GENO.index(f[1])))
                    else:
                        GENOS.append("./.")
                vcflist.append("\t".join([str(locusnumber+1), str(base+1), '.', REF, ",".join(ALT), "20", "PASS",
                                          ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS))
        if not locusnumber % 1000:
            outfile.write( "\n".join(vcflist)+"\n" )
            vcflist = []
                                              

    outfile.write( "\n".join(vcflist) )
    outfile.close()


def most_common(L):
    return max(L, key=L.count)


def unstruct(amb):
    amb = amb.upper()
    " returns bases from ambiguity code"
    D = {"R":["G","A"],
         "K":["G","T"],
         "S":["G","C"],
         "Y":["T","C"],
         "W":["T","A"],
         "M":["C","A"],
         "A":["A","A"],
         "T":["T","T"],
         "G":["G","G"],
         "C":["C","C"],
         "N":["N","N"],
         "-":["-","-"]}
    return D.get(amb)


VCF_HEADER = """##fileformat=VCFv4.1
##source=pyRAD.v.{}
##reference=common_allele_at_each_locus
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
"""

if __name__ == "__main__":
    import sys
    in_loci = sys.argv[1]
    out_vcf = sys.argv[2]
    names_file = sys.argv[3]
    try:
        full_vcf = eval(sys.argv[4])
    except:
        full_vcf = True
    make(in_loci, out_vcf, names_file, full_vcf=full_vcf)

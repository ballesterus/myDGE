#!/usr/bin/env python3

from sys import argv 
import re

def make_dict_gff(gff):
    R= {}
    with open(gff, 'r') as G:
        for line in G:
            if not line.startswith('#'):
                f=line.split("\t")[8]
                gene= re.findall("gene=(.+?);", f)
                prod= re.findall("product=(.+?)[;\n]", f)
                if len(gene) > 0:
                    if  gene in list(R.keys()):
                        if prod not in R[gene[0]] and len(prod) > 0:
                            R[gene].append(prod[0])
                    elif gene not in list(R.keys()) and len(prod) > 0:
                        R[gene[0]]=prod
    return R


def add_annot(afile, gffdict):
    oname = afile.split('.')[0] + "_wanno.csv"
    o=open(oname, 'w')
    with open(afile, 'r') as F:
        h=F.readline().split(',')
        h.insert(1, "Annot")
        o.write(','.join(h))
        for l in F:
            l = l.replace('"', '')
            f=l.split(",")
            q=re.sub("gene-", "", f[0])
            try:
                q=re.sub("gene-", "", f[0])
                a=";".join(gffdict[q])
                f.insert(1, a)
            except:
                print("WARNING %s, not in GFF" %q)
                a = "UNKNOWN"
                f.insert(1, a)
            o.write(",".join(f))
    o.close()


def main(gff, target):
    D=make_dict_gff(gff)
    add_annot(target, D)



if __main_

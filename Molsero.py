#!/usr/bin/env python3

import argparse
from SeqUtils import Fsa
import io
#import gzip
import subprocess as sp
import os

devnull=open(os.devnull,"w")

parser = argparse.ArgumentParser(description='')
parser.add_argument("contigs",type=argparse.FileType())
parser.add_argument("-v","--verbose",action="store_true",default=False,help="Print errors to StdErr")
parser.add_argument("-p","--primers",type=argparse.FileType(), required=True)
args=parser.parse_args()

primerfile=args.primers.name
primers=Fsa("".join(args.primers.readlines()))

contigfile=args.contigs.name
idxh=sp.Popen(("bwa","index",contigfile),stderr=devnull)
idxh.wait()

samOutput=sp.Popen(("bwa","mem","-B","5","-p","-T","15","-k","13",contigfile,primerfile),stdout=sp.PIPE,stderr=devnull).communicate()[0]
samOutput=samOutput.decode().split("\n")
matches=dict()
for line in samOutput:
    if len(line)==0:
        continue
    if line.startswith("@"):
        if args.verbose:
            print(line)
        continue
    field=line.split("\t")
    if(not int(field[1]) & 4):
        matches[field[0]]=field

bands=set()
for k in sorted(matches.keys()):
    if args.verbose:
        print("\t".join(matches[k]))
    else:
        if k.endswith("-F"):
            bands.add(k.rstrip("-F"))
            print("{}\t{}".format(k.rstrip("-F"),abs(int(matches[k][8]))))



if bands== set(('lmo0737','prs')):
    print('IIa')
elif bands== set(('ORF2819','prs')):
    print('IIb')
elif bands== set(('lmo1118','lmo0737','prs')):
    print('IIc')
elif bands== set(('ORF2110','ORF2819','prs')):
    print('IVb')
elif bands== set(('prs')):
    print('L')
else:
    print("Eureka!!")

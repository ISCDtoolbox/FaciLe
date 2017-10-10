#!/usr/bin/python
# -*- coding: utf-8 -*

import os
import sys
import numpy as np
import argparse

def parse():
    parser = argparse.ArgumentParser(description="converts CSV from Lydie to true CSV")
    parser.add_argument("-i", "--inputCSV", help="input CSV", type=str, required=True)
    parser.add_argument("-o", "--outputCSV", help="output CSV", type=str, required=True)
    return parser.parse_args()
def checkArgs(args):
    if not os.path.isfile(args.inputCSV):
        print args.inputCSV + "is not a valid directory"
        sys.exit()
    args.inputCSV  = os.path.abspath(args.inputCSV)
    args.outputCSV = os.path.abspath(args.outputCSV)

if __name__=="__main__":

    args = parse()
    checkArgs(args)
    indsToDelete = [1,4,5,11,13,15]

    with open(args.inputCSV) as f:
        LINES = f.readlines()
        NEWLINES = []

        for i,l in enumerate(LINES[1:]):
            if i%3 == 0:
                NEWLINES.append(l.strip().split(","))
            if i%3 == 1:
                NEWLINES[-1].extend(l.strip().split(","))
            else:
                pass

        for i in range(len(NEWLINES[0]))[::-1]:
            hasSomething = False
            for n in NEWLINES:
                if n[i]!="":
                    hasSomething = True
                    break
            if hasSomething:
                pass
            else:
                NEWLINES = [ [x for j,x in enumerate(n) if j!=i] for n in NEWLINES]

        #NEWLINES = [ [n for j,n in enumerate(l) if j not in indsToDelete] for n in NEWLINES ]

        with open(args.outputCSV, "w") as ff:
            for l in NEWLINES:
                ff.write(",".join(l) + "\n")

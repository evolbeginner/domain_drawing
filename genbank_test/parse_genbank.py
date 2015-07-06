#! /bin/env python

from Bio import SeqIO
import getopt
import sys
import re

##################################################################
input = None
output = None

try:
    opts,args = getopt.getopt(
        sys.argv[1:],
        "hi:o:",
        ["in=","out=",]
    )
except getopt.GetoptError:
    print "illegal error!"
    show_help()

for opt, value in opts:
    if opt == "-i" or opt == "--in":
        input = value
    if opt == "-o" or opt == "--out":
        output = value

##################################################################
for record in SeqIO.parse(open(input, "rU"), "genbank"):
    for feature in record.features:
        #print feature
        if feature.type == "Region":
            if 'db_xref' in feature.qualifiers and re.search('CDD',feature.qualifiers['db_xref'][0]):
                print feature
 


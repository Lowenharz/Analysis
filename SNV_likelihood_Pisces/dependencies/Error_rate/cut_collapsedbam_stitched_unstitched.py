#/illumina/thirdparty/python/python-2.7.5/bin/python
'''
@author: tjiang1
'''

import os
import glob
import re
import argparse
import pysam
import sys
from array import array

def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--Bam", type=str, required=True,
                        help="Stitched bam to calculate noise AF for")
    return parser

parser = build_parser()
args = parser.parse_args()
input = args.Bam
sampleID = re.sub(r".stitched.bam$", "", os.path.basename(input))

bam = pysam.AlignmentFile(input)
bam_stitched = pysam.AlignmentFile(os.path.join(os.path.dirname(input), sampleID + "_st.bam"), "wb", template=bam)
bam_unstitched = pysam.AlignmentFile(os.path.join(os.path.dirname(input), sampleID + "_ust.bam"), "wb", template=bam)

for s in bam:
    try:
        XD = s.get_tag("XD")
    except:
        bam_unstitched.write(s) #output unstitched reads
        continue
    Stitch_num = 0
    left_num = 0
    right_num = 0
    if re.search(r'\d+S', XD)!= None:
        Stitch_num = int(re.findall(r'\d+',re.findall(r'\d+S', XD)[0])[0])
    if re.search(r'\d+[FR].+S', XD)!= None:       #######XD:Z:1F67S4F
        left_num = int(re.findall(r'\d+',re.findall(r'\d+[FR].+S', XD)[0])[0])
    if re.search(r'S+\d+[FR]', XD)!=None:
        right_num = int(re.findall(r'\d+',re.findall(r'S+\d+[FR]', XD)[0])[0])
    if left_num+Stitch_num+right_num == len(s.query_qualities):    #cigar string with D need to be fix later
        s_st=s.__copy__()
        qual=list(s_st.qual)
        qual[0 : (left_num)] = ['!'] * left_num
        qual[( left_num + Stitch_num ) : ] = ['!'] * right_num
        s_st.qual = ''.join(qual)
        s_ust=s.__copy__()
        qual=list(s_ust.qual)
        qual[left_num : ( left_num + Stitch_num )] = ['!'] * Stitch_num
        s_ust.qual = ''.join(qual)
        bam_stitched.write(s_st)
        bam_unstitched.write(s_ust)
        

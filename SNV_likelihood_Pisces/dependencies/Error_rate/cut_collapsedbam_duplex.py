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
sampleID = re.sub(r".bam$", "", os.path.basename(input))
bam = pysam.AlignmentFile(input)
bam_duplex = pysam.AlignmentFile(os.path.join(os.path.dirname(input), sampleID + "_duplex.bam"), "wb", template=bam)
bam_simplex_forward = pysam.AlignmentFile(os.path.join(os.path.dirname(input), sampleID + "_simplex_forward.bam"), "wb", template=bam)
bam_simplex_reverse = pysam.AlignmentFile(os.path.join(os.path.dirname(input), sampleID + "_simplex_reverse.bam"), "wb", template=bam)
bam_simplex_other = pysam.AlignmentFile(os.path.join(os.path.dirname(input), sampleID + "_simplex_other.bam"), "wb", template=bam)

for s in bam:
    flag = s.flag
    XV = s.get_tag("XV")
    XW = s.get_tag("XW")
    if (int(XW)>0):
	bam_duplex.write(s)
    if (int(XW)==0):
	if (flag==99 or flag==147):
        	bam_simplex_forward.write(s)
	elif (flag==83 or flag==163):
		bam_simplex_reverse.write(s)
	else:
		bam_simplex_other.write(s)

"""
Script for calling snv with likelihood model
Input collapsed bam folder or indiviual bam
"""

import os
import sys
import re
import pyflow
import argparse
import multiprocessing
import snv_pyflow_workflow

def build_arg_parser():
    """Builds the argument parser"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-b", "--bams", nargs="*", required=True, type=str)
    parser.add_argument("-m", "--manifest", required=True, type=str)
    parser.add_argument("-o", "--outfolder", required=False, type=str)
    parser.add_argument("-e", '--email', required=True, type=str)
#    parser.add_argument("-pdir", '--pyflowDir', required=True, type=str, help=argparse.SUPPRESS)
    return parser


def main():
    """Runs the script"""
    parser = build_arg_parser()
    args = parser.parse_args()

    inputDir = os.path.dirname(args.bams[0])
    pyflowDir = os.path.join(inputDir,'pyflow')
    errorFile = os.path.join(pyflowDir, "snv_likelihood.log")
    warningFile = os.path.join(pyflowDir, "snv_likelihood.log")
    wf = snv_pyflow_workflow.SNV_calling(
        args.bams,
        args.manifest,
    )
    wf.run(mode='sge', mailTo=args.email,dataDirRoot=pyflowDir,
           errorLogFile=errorFile, warningLogFile=warningFile)


if __name__ == "__main__":
    main()

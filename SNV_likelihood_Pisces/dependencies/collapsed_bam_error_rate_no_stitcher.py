"""
Calculate noiseAF independent of Pisces output. This script will create several output files,
including a stitched and sorted bam file, a *.error.txt that contains a table of per-locus
counts/coverage, and a *.noiseaf.json file which contains summary values, including the desired
noiseAF value. Note that all files will be created in the directory as the input bam file.
Methodology + original script: tjiang1
Modified and integrated into Bevmo: kwu1
"""
import os
import sys
import re
import time
import argparse
import subprocess
import json
import gzip

def median(values):
    """Computes the median of the given values"""
    sorted_values = list(sorted(values))
    if len(sorted_values) == 1:
        return sorted_values[0]
    center = (len(sorted_values) + 1) / 2 - 1
    # Even
    if len(sorted_values) % 2 == 0:
        return float(sorted_values[center] + sorted_values[center + 1]) / 2
    # Odd
    else:
        return sorted_values[center]

def mad(values):
    """Computes the median absolute deviation of the given values"""
    c = 0.67448975019608171
    MTC = median(values)
    values = [abs(x - MTC) for x in values]
    sorted_values = list(sorted(values))
    mad = float(median(sorted_values))/c
    return mad

def run_samtools_mpileup(stitched_bam, manifest, genome_fa, samtools_binary, min_base_q=0):
    """Runs samtools mpileup on the given bam. Returns the stdout output"""
    samplename = re.sub(r".stitched.bam$", "", os.path.basename(stitched_bam))
    argdict = {
        "samtools": samtools_binary,
        "genome": genome_fa,
        "manifest": manifest,
        "bam": stitched_bam,
        "samplename": samplename,
        "min_base_q": min_base_q
    }
    # Refernece command:
    # # samtools mpileup -d 80000 -f /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa -l $manifest -Q min_base_q -B ${sample}.stitched.bam |python $noiseAF $sample    
    pileup_command = "{samtools} mpileup -d 80000 -f {genome} -Q {min_base_q} -l {manifest} -B {bam}".format(**argdict)
    pileup_handle = subprocess.Popen(pileup_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pileup_handle.communicate()
    # Return the standard output since it contains the actual output
    stdout_tokenized = stdout.split("\n") # Split by newline delimiter
    return stdout_tokenized

def calculate_coverage_metrics(pileup_lines, bedFile):
    """
    Given the output of mpiileup splitted by newline chars, calculate the coverage metric
    Note that this pileup should be the output of mpileup WITHOUT base quality cutoff
    """
    Depth_list = {}
    with open(bedFile, 'r') as f_manifest:
        for line in f_manifest:
            l = line.rstrip().split("\t")
            chr = l[0]
            start = int(l[1])
            end = int(l[2])
            for i in range (start+1,end+1):
                key = chr + ':' + str(i) 
                Depth_list[key] = 0
    
    Count_1000X = 0
    Count_1500X = 0
    Count_2000X = 0
    for line in pileup_lines:
        l = line.rstrip().split("\t")
        try:
            chr = l[0]
            pos = l[1]
            DP = l[3]
        except:
            continue
        key = chr + ':' + str(pos)
        Depth_list[key] = int(DP)
        if int(DP) >= 1000:
            Count_1000X += 1
        if int(DP) >= 1500:
            Count_1500X += 1
        if int(DP) >= 2000:
            Count_2000X += 1
    Depth = [x[1] for x in Depth_list.items()]
    MTC = median(Depth)
    try:
        CV_MTC = mad(Depth)/MTC
    except ZeroDivisionError:
        CV_MTC = 0.0
    PCT_1000X = float(Count_1000X)/len(Depth_list)*100
    PCT_1500X = float(Count_1500X)/len(Depth_list)*100
    PCT_2000X = float(Count_2000X)/len(Depth_list)*100
    Coverage_metrics=dict(MTC=MTC, CV_MTC=CV_MTC, PCT_1000X=PCT_1000X, PCT_1500X=PCT_1500X,PCT_2000X = PCT_2000X)
    return Coverage_metrics

def calculate_noise_af(pileup_lines, af_cutoff, error_black_list, outfile):
    """
    Given the output of mpiileup splitted by newline chars, calculate the noise AF
    Note that this should be the output of mpileup WITH base quality cutoff of 65
    """
    ####Define black list to exclude for error calculation
    black_list=[]
    if error_black_list != 'none':
        f_black_list = gzip.open(error_black_list, 'r')
        for line in f_black_list:
            l=line.rstrip().split("\t")
            black_list.append(l[0]+'_'+l[1])
        
        black_list=set(black_list)

    total_AD=0
    total_DP=0
    with open(outfile, 'w') as f_out:
        header=['chr','pos','ref','DP','nA','nC','nG','nT','nIndel','Included']
        f_out.write('\t'.join(header)+'\n')
        for line in pileup_lines:
            # Original logic here by Tingting
            l=line.rstrip().split("\t")
            if (len(l)<5):
                continue
            chromosome=l[0]
            pos=l[1]
            ref=l[2].upper()
            DP=l[3]
            seq=l[4]
            if chromosome+'_'+pos in black_list:
                included=False
            else:
                included=True
            #remove ^
            seq=re.sub('\\^.','',seq)
            ##remove $
            seq=re.sub('\\$','',seq)
            ##remove indel##
            while (re.search('[0-9]+',seq)!=None):
                match_pos=re.search('[0-9]+',seq).span()
                number=int(seq[match_pos[0]:match_pos[1]])
                s=list(seq)
                s[match_pos[0]:(match_pos[1]+number)]=''
                seq=''.join(s)
            ##remove *
            seq=re.sub('\\*','',seq)

            seq=re.sub(',',ref,seq)
            seq=re.sub('\\.',ref,seq)
            seq=seq.upper()
            nA=seq.count('A')
            nC=seq.count('C')
            nG=seq.count('G')
            nT=seq.count('T')
            nIndel=seq.count('+')+seq.count('-')
            RD=seq.count(ref)
            AD=nA+nC+nG+nT+nIndel-RD
            AD2=min(AD,RD)
            try:
                AF = AD2/float(AD+RD)
            except ZeroDivisionError:
                AF = 0.0
            output=[chromosome,pos,ref,DP,str(nA),str(nC),str(nG),str(nT),str(nIndel),str(included)]
            f_out.write('\t'.join(output)+'\n')

            if AF < af_cutoff and included == True:
                total_AD=total_AD+AD2
                total_DP=total_DP+AD+RD
    try:
        VAF=float(total_AD)/total_DP
    except ZeroDivisionError:
        VAF = -1.0 # Set VAF to be -1 when total_DP==0
    noiseAF=dict(NoiseAF=VAF, total_AD=total_AD, total_DP=total_DP)
    return noiseAF


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--Bam", type=str, required=True,
                        help="Stitched bam to calculate noise AF for")
    parser.add_argument("--Manifest", type=str, required=True,
                        help="Manifest interval file")
    parser.add_argument("--GenomeFasta", type=str, required=True,
                        help="Genome fasta file")
    parser.add_argument("--Samtools", type=str, required=True,
                        help="Path to samtools binary")
    parser.add_argument("--af_cutoff", type=float, default=0.05,
                        help="af cutoff for error calculation")
    parser.add_argument("--error_black_list", type=str, default=os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources", "1000G_blacklist_error_rate.vcf.gz"),
                        help="Blacklist for error rate calculation")
    return parser


def main():
    """Run the script"""
    parser = build_parser()
    args = parser.parse_args()

    if not os.path.isfile(args.Bam):
        raise RuntimeError("Cannot find stitched bam %s" % args.Bam)

    if not os.path.isfile(args.Manifest):
        raise RuntimeError("Cannot find manifest file: %s" % args.Manifest)
    
    if not os.path.isfile(args.GenomeFasta):
        raise RuntimeError("Cannot find genome fasta: %s" % args.GenomeFasta)

    if not os.path.isfile(args.error_black_list):
        raise RuntimeError("Cannot find error blacklist file: %s" % args.error_black_list)

    # Calculate the pileups
    pileup_without_base_qual_cutoff = run_samtools_mpileup(args.Bam, args.Manifest, args.GenomeFasta, args.Samtools)
    pileup_with_base_qual_cutoff = run_samtools_mpileup(args.Bam, args.Manifest, args.GenomeFasta, args.Samtools, min_base_q=65)
    # Calculate metrics from pileups
    samplename = re.sub(r".bam$", "", args.Bam)
    coverage_metric_dict = calculate_coverage_metrics(pileup_without_base_qual_cutoff,args.Manifest)
    output_errortable = os.path.join(os.path.dirname(args.Bam), samplename + ".error.txt")
    noise_af_dict = calculate_noise_af(pileup_with_base_qual_cutoff, args.af_cutoff, args.error_black_list, output_errortable)
    metric_dict = coverage_metric_dict.copy()
    metric_dict.update(noise_af_dict)

    output_json = os.path.join(os.path.dirname(args.Bam), samplename + ".noiseAF.json")
    with open(output_json, 'w') as f_out:
        json.dump(metric_dict, f_out, indent=4, sort_keys=True)

if __name__ == "__main__":
    main()

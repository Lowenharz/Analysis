import os
import sys
import multiprocessing
import pyflow
import re
import gzip
import strandErrorMode


def check_workflow_name(x):
    """Fixes workflow names such that they don't cause pyflow to error out"""
    return re.sub(r"(?!([A-Za-z0-9_-])).", "_", x)


class StitcherWorkflow(pyflow.WorkflowRunner):
    """
    Run stithcer on ReCo bam
    """
    def __init__(self, input_bam, manifest, output_folder, nthreads=1):
        self.input_bam = input_bam
        self.output_folder = output_folder
        self.threads = nthreads
        self.manifest = manifest
        self.samplename = os.path.basename(input_bam).split(".")[0]
        self.output_bam = os.path.join(
            self.output_folder,
            self.samplename+'.stitched.bam'
        )
        self.dotnet = '/illumina/thirdparty/dotnet/latest/dotnet'
        # Find the binaries
        self.stitcher_binary = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Stitcher",
            "Stitcher.dll"
        )
        assert os.path.isfile(self.stitcher_binary), "Cannot find stitcher binary at %s" % self.stitcher_binary

        self.noiseAF = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "collapsed_bam_error_rate_no_stitcher.py"
        )

        # Build the commands for stitcher
        cmd = "{0} {1}".format(self.dotnet, self.stitcher_binary)
        cmd += " -Bam {0}".format(self.input_bam)
        cmd += " -OutFolder {0}".format(self.output_folder)
        self.stitcher_commands = cmd
        print self.stitcher_commands
        # Build the commands for bam sort
        self.output_bam_sort_tmp = self.output_bam+'_tmp'
        cmd = "samtools sort {0} -o {1} ;".format(self.output_bam, self.output_bam_sort_tmp)
        cmd += "mv {0} {1} ".format(self.output_bam_sort_tmp, self.output_bam)
        self.samtools_sort_commands = cmd
        print self.samtools_sort_commands
        # Build the commands for bam index
        cmd = "samtools index {0}".format(self.output_bam)
        self.samtools_index_commands = cmd
        print self.samtools_index_commands
        # Build the commands for noiseAF
        self.noise_cmd = "/home/kwu1/software/python/bin/python %s" % self.noiseAF
        self.noise_cmd += " --Bam %s" % self.output_bam
        self.noise_cmd += " --Manifest %s" % self.manifest
        self.noise_cmd += " --Samtools /illumina/thirdparty/samtools/latest/bin/samtools"
        self.noise_cmd += " --GenomeFasta /illumina/sync/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
        self.noise_cmd += " --error_black_list /illumina/scratch/K2-I/Users/tjiang1/JIRA/ONBI-979/Code/resource/1000G_blacklist_error_rate.vcf.gz"
        print self.noise_cmd


    def workflow(self):
        # First run stitcher
        stitcher_label = self.addTask(
            check_workflow_name("varcall_stitcher_{0}".format(self.samplename)),
            self.stitcher_commands,
            nCores=1,  # Stitcher doesn't yet support multithreading in non-beta
            memMb=12*1024  # Allocate 12 GB RAM for stitcher
            )
        samtools_sort_label = self.addTask(
            check_workflow_name("sort_%s" % self.output_bam),
            self.samtools_sort_commands,
            nCores=4,
            memMb=16*1024,
            dependencies=stitcher_label
            )
        samtools_index_label = self.addTask(
            check_workflow_name("index_%s" % self.output_bam),
            self.samtools_index_commands,
            nCores=1,
            memMb=2*1024,
            dependencies=samtools_sort_label
            )
        self.addTask(
            check_workflow_name("noiseAF_%s" % self.output_bam),
            self.noise_cmd,
            nCores=1,
            memMb=60*1024,
            dependencies=samtools_sort_label
            )


class PiscesWorkflow(pyflow.WorkflowRunner):
    """
    Run Pisces and PEPE on stitched bam
    """
    def __init__(self, input_folder, manifest, output_folder, nthreads=4):
        self.input_folder = input_folder
        self.samplename = os.path.basename(self.input_folder)
        self.manifest = manifest
        self.output_folder = output_folder
        self.threads = nthreads
        self.dotnet = '/illumina/thirdparty/dotnet/latest/dotnet'
        self.mono = '/illumina/sync/software/unofficial/Isas/packages/mono-4.0.2/bin/mono'
        # Find the binaries
        self.pisces_binary = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Pisces",
            "Pisces.dll"
        )
        assert os.path.isfile(self.pisces_binary), "Cannot find Pisces binary at %s" % self.pisces_binary

        self.pepe_binary = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Pepe",
            "Pepe.exe"
        )
        assert os.path.isfile(self.pepe_binary), "Cannot find PePe binary at %s" % self.pepe_binary
        # Find the resource file
        self.pepe_baseline = "/illumina/scratch/K2-I/software/lychee_for_bevmo/qscorebaselines/grail_panel_qscore_baseline.txt"
        self.pepe_priors = ["/illumina/scratch/K2-I/Users/tjiang1/JIRA/ONBI-979/Code/resource/grail_panel_1000GenePrior.vcf",
                            "/illumina/scratch/K2-I/Users/tjiang1/JIRA/ONBI-979/Code/resource/grail_panel_cosmic_prior.vcf"]

        # Build the command for pisces
        self.pisces_cmd = "{0} {1}".format(self.dotnet, self.pisces_binary)
        self.pisces_cmd += " -BAMFolder '{0}'".format(self.input_folder)  # Takes stitched bams as input
        self.pisces_cmd += " -I '{0}'".format(self.manifest)  # Manifest file
        self.pisces_cmd += " -G '/illumina/sync/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta'"
        self.pisces_cmd += " -MinVF '0.0001'"
        self.pisces_cmd += " -SSFilter 'False'"
        self.pisces_cmd += " -MinBQ '65'"
        self.pisces_cmd += " -MaxVQ '100'"
        self.pisces_cmd += " -MinDepthFilter '500'"
        self.pisces_cmd += " -MinVQ '0'"
        self.pisces_cmd += " -VQFilter '20'"
        self.pisces_cmd += " -gVCF 'True'"
        self.pisces_cmd += " -ReportNoCalls 'True'"
        self.pisces_cmd += " -CallMNVs 'False'"
        self.pisces_cmd += " -MaxMNVLength '15'"
        self.pisces_cmd += " -MaxGapBetweenMNV '10'"
        self.pisces_cmd += " -OutFolder '{0}'".format(self.output_folder)
        self.pisces_cmd += " -t '{0}'".format(self.threads)
        self.pisces_cmd += " -Collapse 'False'"
        self.pisces_cmd += " -RepeatFilter '8'"
        self.pisces_cmd += " -ld '500'"
        self.pisces_cmd += " -NifyDisagreements 'True'"
        self.pisces_cmd += " -ReportRcCounts 'True'"
        self.pisces_cmd += " -ReportTsCounts 'True'"
        print self.pisces_cmd
        # Create the pepe command
        self.pepe_cmd = "{0} {1}".format(self.mono, self.pepe_binary)
        self.pepe_cmd += " -vcffolder {0}".format(self.output_folder)
        self.pepe_cmd += " -baselinefile {0}".format(self.pepe_baseline)
        self.pepe_cmd += " -mode lowfrequency"
        self.pepe_cmd += " -priorFiles {0}".format(",".join(self.pepe_priors))
        self.pepe_cmd += " -baselineMethod wilcox"
        print self.pepe_cmd
    def workflow(self):
        # pisces + pepe
        Pisces_label = self.addTask(
            check_workflow_name("varcall_pisces_{0}".format(self.samplename)),
            self.pisces_cmd,
            nCores=8,  
            memMb=60*1024  
            )
        PePe_label = self.addTask(
            check_workflow_name("varcall_pepe_{0}".format(self.samplename)),
            self.pepe_cmd,
            nCores=8,  
            memMb=60*1024,
            dependencies=Pisces_label
            )
        


class SNV_calling(pyflow.WorkflowRunner):
    """
    Given bam files
    """
    def __init__(self, bamfiles, manifest_bed):
        self.bamfiles = bamfiles
        self.bamDir = os.path.dirname(self.bamfiles[0])
        self.manifest = manifest_bed

    def workflow(self):
        observed_samplenames = set()
        ####stitcher#####
        stitcher_labels = []
        for bamfile in self.bamfiles:
            assert os.path.isfile(bamfile), "Cannot find input bam: %s" % bamfile
            # Parse samplename and ensure uniqueness
            samplename = os.path.basename(bamfile).split(".")[0]
            assert samplename not in observed_samplenames, "Duplicated samplename: %s" % samplename
            observed_samplenames.add(samplename)

            Stitched_output_folder = os.path.join(
                self.bamDir,
                'StitchedBam'
            )
            Stitcher_wf = StitcherWorkflow(bamfile, self.manifest, Stitched_output_folder)
            Stitcher_label = self.addWorkflowTask(
                check_workflow_name("dispatch_Stitcher_%s" % samplename),
                Stitcher_wf,
            )
            stitcher_labels.append(Stitcher_label)

        # variant call
        Pisces_output_folder = os.path.join(
            self.bamDir,
            'Variants'
        )
        Pisces_wf = PiscesWorkflow(Stitched_output_folder, self.manifest, Pisces_output_folder)
        Pisces_label = self.addWorkflowTask(
        check_workflow_name("dispatch_Pisces_%s" % samplename),
            Pisces_wf,
            dependencies=stitcher_labels
        )


            

if __name__ == "__main__":
    # wf = RecoTranslocationWorkflow(
    #     sys.argv[1],
    #     "foo"
    # )
    # wf.run(mode='sge')

    # wf.run(mode='sge')

    pass
    

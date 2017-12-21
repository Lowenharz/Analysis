import os.path
import sys
import re

# Specify the pyflow work directory

from pyflow import WorkflowRunner
def check_workflow_name(x):
    """Fixes workflow names such that they don't cause pyflow to error out"""
    return re.sub(r"(?!([A-Za-z0-9_-])).", "_", x)




def ensureDir(d):
    """
    make directory if it already doesn't exists
    """

    if os.path.exists(d):
        if not os.path.isdir(d):
            raise Exception("Can't create directory: %s" %(d))
    else:
        os.makedirs(d)


class strandSpecificErrorMode(WorkflowRunner):
    """
    Takes a collapsed BAM as an input, splits it by strand, then stitch them and finally
    separate them by support status. Finally calcualte the error rate in strand splitted 
    and support separated files.

    params:
    inputBamDir
    """

    def __init__(self, collapsedBAM, manifest):
        self.bam = collapsedBAM
        self.sample = os.path.splitext(os.path.split(self.bam)[1])[0]
        self.path = os.path.split(collapsedBAM)[0]
        self.manifest = manifest
        self.Dir =  os.path.dirname(self.bam)
        self.stitchedDir = os.path.join(
            self.Dir,
            'StitchedBam'
        )
        ensureDir(self.stitchedDir)

    def workflow(self):
        splitterWF = StrandSplitter(self.bam)
        splitterLabel = self.addWorkflowTask(check_workflow_name("split_" + self.sample), splitterWF)
        
        stitcherWF0 = StrandStitcher(self.path+'/'+self.sample + ".bam")
        stitcherLabel0 = self.addWorkflowTask(check_workflow_name("stitch_0_" + self.sample), stitcherWF0, dependencies =  splitterLabel)
        
        stitcherWF1 = StrandStitcher(self.path+'/'+self.sample + "_duplex.bam")
        stitcherLabel1 = self.addWorkflowTask(check_workflow_name("stitch_1_" + self.sample), stitcherWF1, dependencies =  splitterLabel)
        stitcherWF2 = StrandStitcher(self.path+'/'+self.sample + "_simplex_forward.bam")
        stitcherLabel2 = self.addWorkflowTask(check_workflow_name("stitch_2_" + self.sample), stitcherWF2, dependencies =  splitterLabel)
        stitcherWF3 = StrandStitcher(self.path+'/'+self.sample + "_simplex_reverse.bam")
        stitcherLabel3 = self.addWorkflowTask(check_workflow_name("stitch_3_" + self.sample), stitcherWF3, dependencies =  splitterLabel)
        

        stitchSeparatorWF1 = stitchSeparator(self.path+'/'+self.sample + "_duplex.stitched.bam")
        stitchSeparatorLabel1 = self.addWorkflowTask(check_workflow_name("stitchSeparate_1_" + self.sample), stitchSeparatorWF1, dependencies = stitcherLabel1)
        stitchSeparatorWF2 = stitchSeparator(self.path+'/'+self.sample + "_simplex_forward.stitched.bam")
        stitchSeparatorLabel2 = self.addWorkflowTask(check_workflow_name("stitchSeparate_2_" + self.sample), stitchSeparatorWF2, dependencies = stitcherLabel2)
        stitchSeparatorWF3 = stitchSeparator(self.path+'/'+self.sample + "_simplex_reverse.stitched.bam")
        stitchSeparatorLabel3 = self.addWorkflowTask(check_workflow_name("stitchSeparate_3_" + self.sample), stitchSeparatorWF3, dependencies = stitcherLabel3)
        
        errorLociIdentifierWF0 = ErrorLociIdentifier(self.path+'/'+self.sample + ".stitched.bam", self.manifest)
        errorLociIDLabel0 = self.addWorkflowTask(check_workflow_name("errorIdentification_0_" + self.sample), errorLociIdentifierWF0, dependencies = stitcherLabel0)
        errorLociIdentifierWF1 = ErrorLociIdentifier(self.path+'/'+self.sample + "_duplex_ust.bam",self.manifest)
        errorLociIDLabel1 = self.addWorkflowTask(check_workflow_name("errorIdentification_1_" + self.sample), errorLociIdentifierWF1, dependencies = stitchSeparatorLabel1)
        errorLociIdentifierWF2 = ErrorLociIdentifier(self.path+'/'+self.sample + "_duplex_st.bam",self.manifest)
        errorLociIDLabel2 = self.addWorkflowTask(check_workflow_name("errorIdentification_2_" + self.sample), errorLociIdentifierWF2, dependencies = stitchSeparatorLabel1)
        errorLociIdentifierWF3 = ErrorLociIdentifier(self.path+'/'+self.sample + "_simplex_forward_ust.bam",self.manifest)
        errorLociIDLabel3 = self.addWorkflowTask(check_workflow_name("errorIdentification_3_" + self.sample), errorLociIdentifierWF3, dependencies = stitchSeparatorLabel2)
        errorLociIdentifierWF4 = ErrorLociIdentifier(self.path+'/'+self.sample + "_simplex_forward_st.bam",self.manifest)
        errorLociIDLabel4 = self.addWorkflowTask(check_workflow_name("errorIdentification_4_" + self.sample), errorLociIdentifierWF4, dependencies = stitchSeparatorLabel2)
        errorLociIdentifierWF5 = ErrorLociIdentifier(self.path+'/'+self.sample + "_simplex_reverse_ust.bam",self.manifest)
        errorLociIDLabel5 = self.addWorkflowTask(check_workflow_name("errorIdentification_5_" + self.sample), errorLociIdentifierWF5, dependencies = stitchSeparatorLabel3)
        errorLociIdentifierWF6 = ErrorLociIdentifier(self.path+'/'+self.sample + "_simplex_reverse_st.bam",self.manifest)
        errorLociIDLabel6 = self.addWorkflowTask(check_workflow_name("errorIdentification_6_" + self.sample), errorLociIdentifierWF6, dependencies = stitchSeparatorLabel3)

        errorCalculatorWF = ErrorRateCalculator(self.path+'/'+self.sample+".bam")
        errorCalculatorLabel=self.addWorkflowTask(check_workflow_name("errorRateCalc_" + self.sample), errorCalculatorWF, dependencies = set([errorLociIDLabel0, errorLociIDLabel1, errorLociIDLabel2, errorLociIDLabel3, errorLociIDLabel4, errorLociIDLabel5, errorLociIDLabel6]))
        ##mv all stitched bam to stitched folder##
        self.mv_cmd = "mv {0}/{1}*stitched.bam* {2}".format(self.path,self.sample,self.stitchedDir)
        self.addTask(
            check_workflow_name('cleanstitchedbam_'+self.sample),
            self.mv_cmd,
            isForceLocal=True,
            dependencies = errorCalculatorLabel)


        


class StrandSplitter(WorkflowRunner):
    """
    Given a collpsed BAM file split the file into F1R2 and F2R1 strand
    """
    
    def __init__(self, bam):
        """
        params: 
        BAMfile
        """
        
        self.bamFile = bam
        self.sample = os.path.splitext(os.path.split(self.bamFile)[1])[0]
        self.python = '/home/kwu1/software/python/bin/python'
        self.splitterScript = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Error_rate",
            "cut_collapsedbam_duplex.py"
        )
        self.strandSplitCmd = "{0} {1} --Bam {2}".format(self.python, self.splitterScript, self.bamFile)

    def workflow(self):
        bamDir = os.path.dirname(self.bamFile)
        ensureDir(bamDir)
                                                                                
        (bamPrefix, bamExt) = os.path.splitext(self.bamFile)

        if bamExt != ".bam":
            raise Exception("bamFile argument must end in '.bam'. bamFile is: %s" % (bamFile))
        if bamPrefix=="":
            raise Exception("bamFile argument must have a prefix before the '.bam' extension")

        self.addTask(
            check_workflow_name('strandSplit_'+os.path.basename(self.bamFile)),
            self.strandSplitCmd,
            nCores=1,
            memMb=32*1024)


class StrandStitcher(WorkflowRunner):
    """
    Given a strand separated file stitches them together
    """
    
    def __init__(self, bam):
        """
        params: 
        strandedBAMfile
        outputfolder
        """
        
        self.bamFile = bam
        self.outputFolder = os.path.split(bam)[0]
        self.dotnet = '/illumina/thirdparty/dotnet/latest/dotnet'
        self.stitcher_binary = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Stitcher",
            "Stitcher.dll"
        )
        self.sample = os.path.splitext(os.path.split(self.bamFile)[1])[0]
        self.inputBAM = self.outputFolder+'/'+self.sample+'.stitched.bam'
        self.outBAM = self.outputFolder+'/'+self.sample+'.sorted.stitched.bam'
        cmd = "{0} {1}".format(self.dotnet, self.stitcher_binary)
        cmd += " -Bam {0}".format(self.bamFile)
        cmd += " -OutFolder {0}".format(self.outputFolder)
        self.strandStitchCmd1 = cmd
        self.strandStitchCmd2 = "samtools sort -o {0} {1};".format(self.outBAM,self.inputBAM)
        self.strandStitchCmd2 += "mv {0} {1};".format(self.outBAM,self.inputBAM)
        self.strandStitchCmd2 += "samtools index {0}".format(self.inputBAM)       
    def workflow(self): 
    
        stitchCommand = self.addTask(
            check_workflow_name('strandStitch_' + os.path.basename(self.bamFile)),
            self.strandStitchCmd1,
            nCores=1,
            memMb=32*1024)
        sortCommand = self.addTask(
            check_workflow_name('strandStitchSort_' + os.path.basename(self.bamFile)),
            self.strandStitchCmd2,
            dependencies = stitchCommand,
            nCores=1,
            memMb=32*1024)


class stitchSeparator(WorkflowRunner):
    """
    Given a strand separated file stitches them together
    """
    
    def __init__(self, bam):
        """
        params: 
        strandedBAMfile
        """
        
        self.bamFile = bam
        self.stitchSeparator = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Error_rate",
            "cut_collapsedbam_stitched_unstitched.py"
        )
        self.python = '/home/kwu1/software/python/bin/python'
        self.stitchSeparatorCmd = "{0} {1} --Bam {2}".format(self.python,self.stitchSeparator,self.bamFile)

    def workflow(self): 
        self.addTask(
            check_workflow_name('stitchSeparator_' + os.path.basename(self.bamFile)),
            self.stitchSeparatorCmd,
            nCores=1,
            memMb=32*1024)


class ErrorLociIdentifier(WorkflowRunner):
    """
    Given a bam file separated by stitch status/support status/strand status generate the error.txt 
    """
    
    def __init__(self, bam, manifest):
        self.stratBAM = bam
        self.manifest = manifest
        self.python = '/home/kwu1/software/python/bin/python'
        self.noiseAF = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Error_rate",
            "collapsed_bam_error_rate_no_stitcher.py"
        )
        self.blacklist = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Error_rate",
            "1000G_blacklist_error_rate.vcf.gz"
        )
        # Build the commands for noiseAF
        self.noise_cmd = "{0} {1}".format(self.python, self.noiseAF)
        self.noise_cmd += " --Bam %s" % self.stratBAM
        self.noise_cmd += " --Manifest %s" % self.manifest
        self.noise_cmd += " --Samtools /illumina/thirdparty/samtools/latest/bin/samtools"
        self.noise_cmd += " --GenomeFasta /illumina/sync/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
        self.noise_cmd += " --error_black_list {0}".format(self.blacklist)
    def workflow(self):
    
        self.addTask(
            check_workflow_name('errorRateID_' + os.path.basename(self.stratBAM)),
            self.noise_cmd,
            nCores=1,
            memMb=32*1024)



class ErrorRateCalculator(WorkflowRunner):
    """
    Given a strand splitted and stitch separated BAM files generate the 
    error file
    """
                     
    def __init__(self, bam):
        """
        params: 
        BAMfile
        sampleDirectory
        """
        self.sampleName = os.path.splitext(os.path.split(bam)[1])[0]
        self.sampleDirectory = os.path.split(bam)[0]
        self.errorRateCalculator = os.path.join(
            os.path.dirname(__file__),
            "dependencies",
            "Error_rate",
            "errorRateCalculation.R"
        )
        self.ErrorCalcCmd = "Rscript --vanilla {0} -s {1} -d {2}".format(self.errorRateCalculator,self.sampleName,self.sampleDirectory)

    def workflow(self):
        self.addTask(
            check_workflow_name('errorRateCalc_' + os.path.basename(self.sampleName)),
            self.ErrorCalcCmd,
            nCores=1,
            memMb=32*1024)

#def main():
#    """Runs the workflow"""
#    master_workflow = strandSpecificErrorMode()
#    master_workflow.run(dataDirRoot = "strandErrorLog",
#                        nCore = 2,
#                        memeMb = 4*1024
#                        mailTo = "dkangeyan@illumina.com",
#                        mode = "sge")

if __name__ == "__main__":
    pass

__author__ = 'benkjohnson'

import check_dependencies_linux
import optparse
import sys


cd = check_dependencies_linux.CheckDependencies()


optParser = optparse.OptionParser(
  usage="python SPARTA.py [options]",

  description="Simple Program for Automated reference-based bacterial RNA-seq Transcriptome Analysis (SPARTA)",

  epilog="Written by Benjamin K. Johnson (john3434@msu.edu), Michigan State University Department of " +
  "Microbiology and Molecular Genetics. (c) 2015")

#Future functionality
# optParser.add_option("--SE", help="Single-end read input. Default input choice is single-end if nothing is specified",
#                     action="store_true", default="True", dest="seqtype")
# optParser.add_option("--PE", help="Paired-end read input. Must have the exact same file name and end with _F for the forward read and _R for the reverse read",
#                     action="store_false", default="False", dest="seqtype")
optParser.add_option("--cleanup", help="Clean up the intermediate files to save space. Default action is to retain the intermediate files.", action="store_true", dest="cleanup")
optParser.add_option("--verbose", help="Display more output for each step of the analysis.", action="store_true", dest="verbose")
optParser.add_option("--noninteractive", help="Non-interactive mode. This is for running SPARTA without any user input. Assumes data is on the desktop. If this"
                    " option is specified, you must fill out the configuration file (ConfigFile.txt) with the appropriate experimental conditions in the SPARTA folder.", action="store_true", dest="noninteractive")
optParser.add_option("--threads", help="Define the number of threads that SPARTA should run with. This will enable some speed-up on multi-processor machines. As a generality, define the number of threads as the same number of cores in your computer. Default is 2.",
                       action="store", type="int", default=2, dest="threads")
optParser.add_option("--procs", help="Assign how many processors to use. This is most useful with --noninteractive and on a high performance computing environment. Usage: --procs=<integer value of processors>",
                       action="store", type="int", default=1, dest="procs")


trim = optparse.OptionGroup(optParser, 'Trimmomatic options', "The order the options will be run are: ILLUMINACLIP, LEADING, TRAILING, SLIDINGWINDOW, MINLEN")
trim.add_option("--clip", help="ILLUMINACLIP options. MiSeq & HiSeq usually TruSeq3.fa; GAII usually TruSeq2.fa. Default is ILLUMINACLIP:TruSeq3-SE.fa:2:30:10. Usage: --clip=<adapterseqs>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>",
                  action="store", default="TruSeq3-SE.fa:2:30:10", dest="illuminaclip")
trim.add_option("--lead", help="Set the minimun quality required to keep a base. Default is LEADING=3. Usage: --lead=<quality>",
                  action="store", type="int", default=3, dest="leading")
trim.add_option("--trail", help="Set the minimum quality required to keep a base. Default is TRAILING=3. Usage: --trail=<quality>",
                  action="store", type="int", default=3, dest="trailing")
trim.add_option("--slidewin", help="SLIDINGWINDOW options. Default is SLIDINGWINDOW:4:15. Usage: --slidewin=<window_size>:<required_quality>",
                  action="store", default="4:15", dest="slidingwindow")
trim.add_option("--minlentrim", help="Set the minimum read length to keep in base pairs. Default is 36. Usage: --minlentrim=<readlength>",
                  action="store", type="int", default=36, dest="minlentrim")

bowtie = optparse.OptionGroup(optParser, 'Bowtie options')
bowtie.add_option("--mismatch", help="Output alignments with at most a defined number of mismatches. Usage: --mismatch=<integer_value>",
                    action="store", type="int", dest="mismatch")
bowtie.add_option("--otherbowtieoptions", help="Bowtie has so many options that it is not worth listing them here. Go to http://bowtie-bio.sourceforge.net/manual.shtml#command-line for the manual and all available options. Usage: --otherbowtieoptions='all options inputed as a string (note the quotes!)'",
                    action="store", type="str", dest="otherbowtieoptions")

htseq = optparse.OptionGroup(optParser, 'HTSeq options')
htseq.add_option("--stranded", help="Stranded options: yes, no, reverse. Default is --stranded=reverse. Usage: --stranded=yes/no/reverse",
                   action="store", default="reverse", dest="stranded")
htseq.add_option("--order", help="Order options: name, pos. Usage: --order=name/pos.", action="store", type="str", default="name",  dest="order")
htseq.add_option("--minqual", help="Skip all reads with quality lower than the given value. Default is --minqual=10. Usage: --minqual=<value>",
                   action="store", type="int", default=10, dest="minqual")
htseq.add_option("--type", help="The feature type (3rd column in GTF file) to be used. Default is --type=exon (suitable for RNA-seq analysis)",
                    action="store", type="str", default="exon", dest="type")
htseq.add_option("--idattr", help="Feature ID from the GTF file to identify counts in the output table Default is --idattr=gene_id. Usage: --idattr=<id attribute>",
                   action="store", default="gene_id", dest="idattr")
htseq.add_option("--mode", help="Mode to handle reads overlapping more than one feature. Default is --mode=union. Usage: --mode=union/intersection-strict/intersection-nonempty",
                   action="store", default="union", dest="mode")

optParser.add_option_group(trim)
optParser.add_option_group(bowtie)
optParser.add_option_group(htseq)
(options, args) = optParser.parse_args()

# args = parser.parse_args()

if options.noninteractive:
    cond_lst, data_path = cd.parseConfigFile(options)
else:
    #Welcome the user to the software and check dependencies
    print "Welcome to SPARTA!"
    print "Let's make sure we have everything we need to get started..."
    print "Now checking dependencies...\n"

#Check for Java, R, and NumPy

javacheck = cd.checkjava()
Rcheck = cd.checkR()
numpycheck = cd.checknumpy(options)
# matplotlibcheck = cd.checkmatplotlib()

#If NumPy can't be found, SPARTA will attempt to download and install it

if numpycheck == False:
    cd.installdependencies()
    answer = cd.getanswerstate()
    if answer:
        if numpycheck == False:
            cd.getNumPy(options)
            cd.installNumPy(options)
        # elif numpycheck == False and matplotlibcheck == True:
        #     cd.getNumPy()
        #     cd.installNumPy()
        # elif numpycheck == True and matplotlibcheck == False:
        #     cd.getmatplotlib()
        #     cd.installmatplotlib()

print "Everything appears to be fine. Moving on.\n"

import qc_analysis
import mapping_and_counting
import differential_expression

qc = qc_analysis.QC_analysis()
mac = mapping_and_counting.Mapping_and_Counting()
de = differential_expression.DifferentialExpression()

if options.noninteractive:
    #Create a folder called RNAseq_Data with which to store all of the data analysis

    subfolderpath = qc.create_folder()

    #Read in the data and sort out the genome feature file and genome sequence files

    qc.findreferencefiles(data_path)

    #Trimmomatic

    qc.trimmomatic(data_path, subfolderpath, options)

    #FastQC

    qc.fastqc(data_path, subfolderpath, options)

    #Bowtie

    mac.bowtie(data_path, subfolderpath, options)

    #HTSeq

    mac.htseq(subfolderpath, options)

    #edgeR

    de.de_analysis_noninteractive(subfolderpath, cond_lst)

else:
    #Create a folder called RNAseq_Data with which to store all of the data analysis

    subfolderpath = qc.create_folder()

    #Find the RNAseq data folder location

    rawdatapath = qc.finddata()

    #Read in the data and sort out the genome feature file and genome sequence files

    qc.findreferencefiles(rawdatapath)

    #Trimmomatic

    qc.trimmomatic(rawdatapath, subfolderpath, options)

    #FastQC

    qc.fastqc(rawdatapath, subfolderpath, options)

    #Bowtie

    mac.bowtie(rawdatapath, subfolderpath, options)

    #HTSeq

    mac.htseq(subfolderpath, options)

    #edgeR

    de.de_analysis(subfolderpath)













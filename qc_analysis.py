__author__ = 'benkjohnson'

import os
import subprocess
import datetime
import check_dependencies_linux


class QC_analysis(object):
    def __init__(self):
        cd = check_dependencies_linux.CheckDependencies()
        self._mydesktoppath = cd.getdesktoppath()

    def finddata(self):
        """Attempt to locate the data folder on the user's desktop"""

        data_location_question = str(raw_input("Is the RNAseq data in a folder on the Desktop? (Y or N):"))

        if data_location_question.upper() == "N" or data_location_question.upper() == "NO":
            self._data_path = str(raw_input("Please input the full path to the folder containing your data:"))
            while not os.path.isdir(self._data_path):
                new_userdefined_location = str(raw_input("Sorry, we could not find that folder, try typing the path again or Q to quit:"))
                if new_userdefined_location.upper() == "Q":
                    quit()
                self._data_path = new_userdefined_location
        else:
            foldername = str(raw_input("What is the name of the folder on the Desktop containing the RNAseq data?:"))
            self._data_path = os.path.join(self._mydesktoppath, foldername)
            while not os.path.isdir(self._data_path):
                new_input = str(raw_input("Sorry, we could not find that folder on the desktop, try typing the folder name again or Q to quit:"))
                if new_input.upper() == "Q":
                    quit()
                self._data_path = os.path.join(self._mydesktoppath, new_input)
        return self._data_path

    def create_folder(self):
        """It creates a default folder 'RNAseq_Data' on the desktop if it
        doesn't already exist. And then creates a subfolder in RNAseq_Data with the date
        to store all of the data from the run.
        """
        print "Creating a folder with which all the generated data analysis will be placed."
        print "Default subfolder location will be in 'RNAseq_Data' which is located on the Desktop"
        default_folder_loc = os.path.join(self._mydesktoppath, "RNAseq_Data")
        now = datetime.datetime.now()
        #all data folder will go into HTS_Data, if it doesn't already exist from previous runs it wil be created
        if not os.path.isdir(default_folder_loc):
            os.mkdir(default_folder_loc)
        dateofrun = now.strftime("%Y-%m-%d")
        defined_full_path = os.path.join(default_folder_loc, dateofrun)  #path to the subfolder inside RNAseq_Data
        #if the specified folder is a valid directory name then it will create the folder and exit the function
        if not os.path.isdir(defined_full_path):
            os.mkdir(defined_full_path)
            os.mkdir(os.path.join(defined_full_path, "QC"))
            os.mkdir(os.path.join(defined_full_path, "Bowtie"))
            os.mkdir(os.path.join(defined_full_path, "HTSeq"))
            os.mkdir(os.path.join(defined_full_path, "DEanalysis"))
        else:
            taken = dateofrun
            num = 1
            dateofrun += "_" + str(num)
            defined_full_path = os.path.join(default_folder_loc, dateofrun)
            while os.path.isdir(defined_full_path):
                num += 1
                dateofrun = taken + "_" + str(num)
                defined_full_path = os.path.join(default_folder_loc, dateofrun)
            print "'" + taken + "' was already a directory in RNAseq_Data, '" + dateofrun + "' was used instead for this run."
            os.mkdir(defined_full_path)
            os.mkdir(os.path.join(defined_full_path, "QC"))
            os.mkdir(os.path.join(defined_full_path, "Bowtie"))
            os.mkdir(os.path.join(defined_full_path, "HTSeq"))
            os.mkdir(os.path.join(defined_full_path, "DEanalysis"))
        return defined_full_path

    def findreferencefiles(self, inputdata):
        """Find the reference .gtf and fasta file
        located in the raw data folder."""

        genomefeaturefile = None
        referencegenome = None
        for file in os.listdir(inputdata):
            extension = os.path.splitext(file)[1]
            if extension == ".gff" or extension == ".gtf":
                genomefeaturefile = os.path.join(inputdata, file)
            elif extension == ".fa" or extension == ".fasta":
                referencegenome = os.path.join(inputdata, file)
        if genomefeaturefile == None or referencegenome == None:
            print "You need to include a genome feature file (.gff or .gtf) and genome reference file (.fasta or .fa) in the raw data folder"
            print "Quitting"
            quit()
        return genomefeaturefile, referencegenome

    def trimmomatic(self, datalocation, analysislocation, options):
        """Run trimmomatic for SE reads and add the file prefix
        'trim' to the file name."""

        cd = check_dependencies_linux.CheckDependencies()
        os.chdir(os.path.join(cd.getSPARTAdir(options), "QC_analysis"))
        if not os.path.lexists(os.path.join(cd.getSPARTAdir(options), "QC_analysis", "Trimmomatic-0.33")):
            #This will be a problem for Windows. Just distribute with unzipped binaries?
            subprocess.call(["unzip", "Trimmomatic-0.33.zip"], stdout=open(os.devnull, 'wb'))
        os.chdir(os.path.join(cd.getpwd(), "Trimmomatic-0.33"))
        for file in os.listdir(datalocation):
            extension = file.split(".")[1]
            if extension == "fastq" or extension == "fq":
                subprocess.Popen("java -jar trimmomatic-0.33.jar SE -threads {threads} ".format(threads=options.threads) + os.path.join(datalocation, file) + " " + os.path.join(analysislocation, "QC", "trimmed" + file) + " ILLUMINACLIP:" + cd.getpwd() + "/adapters/{illuminaclip} LEADING:{leading} TRAILING:{trailing} SLIDINGWINDOW:{slidingwindow} MINLEN:{minlen}".format(illuminaclip=options.illuminaclip, leading=options.leading, trailing=options.trailing, slidingwindow=options.slidingwindow, minlen=options.minlentrim), shell=True).wait()

        return

    def fastqc(self, datalocation, analysislocation, options):
        """Run FastQC for trimmed data files."""

        cd = check_dependencies_linux.CheckDependencies()
        os.chdir(os.path.join(cd.getSPARTAdir(options), "QC_analysis"))
        if not os.path.lexists(cd.getSPARTAdir(options) + "/QC_analysis/FastQC"):
            subprocess.call(["unzip", "fastqc_v0.11.3.zip"], stdout=open(os.devnull, 'wb'))
        os.chdir(os.path.join(cd.getpwd(), "FastQC"))
        subprocess.call("chmod 755 fastqc", shell=True)
        print "FastQC is assessing your data set for overall quality"
        for file in os.listdir(datalocation):
            extension = file.split(".")[1]
            if extension == "fastq" or extension == "fq":
                if options.verbose:
                    subprocess.Popen("./fastqc " + os.path.join(analysislocation, "QC", "trimmed" + file), shell=True).wait()
                else:
                    subprocess.Popen("./fastqc --quiet " + os.path.join(analysislocation, "QC", "trimmed" + file), shell=True).wait()
        return
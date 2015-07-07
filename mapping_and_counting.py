__author__ = 'benkjohnson'


import os
import subprocess
import glob
from shutil import copy
import check_dependencies_linux
import qc_analysis

class Mapping_and_Counting(object):
    def __init__(self):

        return

    def bowtie(self, datalocation, analysislocation, options):
        """Run Bowtie for SE reads less than 50 bp in length.
        Will add the ability to run Bowtie2 for PE and SE with
        reads greater than 50 bp."""

        cd = check_dependencies_linux.CheckDependencies()
        qc = qc_analysis.QC_analysis()
        gff, genref = qc.findreferencefiles(datalocation)
        copy(genref, os.path.join(analysislocation, 'Bowtie'))
        genrefname = genref.split("/")[-1]
        copy(gff, os.path.join(analysislocation, 'HTSeq'))
        # subprocess.Popen("cp " + genref + " " + analysislocation + "/Bowtie", shell=True).wait()
        # subprocess.Popen("cp " + gff + " " + analysislocation + "/HTSeq", shell=True).wait()
        os.chdir(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting"))
        # os.chdir(cd.getSPARTAdir() + "/Mapping_and_counting")
        if not os.path.lexists(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting", "bowtie-1.1.1")):
            #This will be a problem for Windows users. Distribute with unzipped binaries?
            subprocess.call(["unzip", "bowtie-1.1.1-linux-x86_64.zip"], stdout=open(os.devnull, 'wb'))
        os.chdir(os.path.join(cd.getpwd(), "bowtie-1.1.1"))
        for file in os.listdir(os.path.join(analysislocation, "QC")):
            extension = file.split(".")[-1]
            if extension == "gz":
                subprocess.Popen("gunzip -c " + os.path.join(analysislocation, "QC", file) + " > " + os.path.join(analysislocation, "Bowtie", os.path.splitext(file)[0]), shell=True).wait()
            else:
                copy(os.path.join(analysislocation, "QC", file), os.path.join(analysislocation, "Bowtie"))

        if options.cleanup:
            for file in os.listdir(os.path.join(analysislocation, "QC")):
                extension = file.split(".")[-1]
                if extension == "gz" or extension in ["fq, fastq"]:
                    subprocess.Popen("rm " + os.path.join(analysislocation, "QC", "{file}".format(file=file)), shell=True).wait()

        print "Building the Bowtie index from the reference genome"
        if options.verbose:
            subprocess.Popen("./bowtie-build -f " + os.path.join(analysislocation, "Bowtie", genrefname) + " " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0], shell=True).wait()
        else:
            subprocess.Popen("./bowtie-build -q -f " + os.path.join(analysislocation, "Bowtie", genrefname) + " " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0], shell=True).wait()
        allebwtfiles = glob.glob("*.ebwt")[:]
        for ebwtfile in allebwtfiles:
            copy(ebwtfile, os.path.join(analysislocation, "Bowtie"))
            # subprocess.Popen("cp " + ebwtfile + " " + analysislocation + "/Bowtie/", shell=True).wait()
        print "Mapping reads to the reference genome with Bowtie"
        for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
            extension = os.path.splitext(file)[1]
            if extension == ".fq" or extension == ".fastq":
                fname = os.path.splitext(file)[0]
                strippedfile = fname[len('trimmed'):]
                if options.verbose:
                    subprocess.Popen("./bowtie -S --threads {threads} -p {procs} ".format(threads=options.threads, procs=options.procs) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                elif options.mismatch:
                    subprocess.Popen("./bowtie -S --threads {threads} -p {procs} -v {mismatch} ".format(threads=options.threads, procs=options.procs, mismatch=options.mismatch) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                elif options.otherbowtieoptions:
                    subprocess.Popen("./bowtie -S --threads {threads} -p {procs} {otherbowtieoptions} ".format(threads=options.threads, procs=options.procs, otherbowtieoptions=options.otherbowtieoptions) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                elif options.mismatch and options.otherbowtieoptions:
                    subprocess.Popen("./bowtie -S --threads {threads} -p {procs} -v {mismatch} {otherbowtieoptions} ".format(threads=options.threads, procs=options.procs, mismatch=options.mismatch, otherbowtieoptions=options.otherbowtieoptions) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                else:
                    subprocess.Popen("./bowtie -S --quiet --threads {threads} ".format(threads=options.threads) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
        return

    def htseq(self, analysislocation, options):
        """Run htseq-count to count gene features post-Bowtie mapping"""

        cd = check_dependencies_linux.CheckDependencies()
        os.chdir(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting"))
        if not os.path.lexists(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting", "HTSeq-0.6.1")):
            subprocess.Popen("tar -zxf HTSeq-0.6.1.tar.gz", stdout=open(os.devnull, 'w'), shell=True).wait()
        # htseqcheck = cd.checkhtseq()
        # if htseqcheck == False:
        os.chdir(os.path.join(cd.getpwd(), "HTSeq-0.6.1"))
        subprocess.Popen("python setup.py build install --user", shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb')).wait()
        os.chdir(os.path.join(cd.getpwd(), "build", "scripts-2.7"))
        gff = glob.glob(os.path.join(analysislocation, "HTSeq") + "/*.g*")[0]
        print "Counting gene features with HTSeq"
        for mapfile in os.listdir(os.path.join(analysislocation, "Bowtie")):
            extension = os.path.splitext(mapfile)[1]
            if extension == ".sam":
                fname = os.path.splitext(mapfile)[0]
                strippedmapfile = fname[len('align'):]
                if options.verbose:
                    subprocess.Popen("./htseq-count --mode={mode} --stranded={stranded} --order={order} --type={type} -a {minqual} --idattr={idattr} ".format(mode=options.mode, stranded=options.stranded, order=options.order, type=options.type, minqual=options.minqual, idattr=options.idattr) + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()
                else:
                    subprocess.Popen("./htseq-count --quiet --mode={mode} --stranded={stranded} --order={order} --type={type} -a {minqual} --idattr={idattr} ".format(mode=options.mode, stranded=options.stranded, order=options.order, type=options.type, minqual=options.minqual, idattr=options.idattr) + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()

        if options.cleanup:
            for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
                subprocess.Popen("rm " + os.path.join(analysislocation, "Bowtie", "{file}".format(file=file)), shell=True).wait()

        return
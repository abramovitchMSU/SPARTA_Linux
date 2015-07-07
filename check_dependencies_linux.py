__author__ = 'benkjohnson'

import imp
import os
import glob
import subprocess

class CheckDependencies(object):
    def __init__(self):
        self._answerstate = False


    def installdependencies(self):
        """Ask if the user would like to install dependencies."""

        self._answer = str(raw_input("Would you like SPARTA to try and download and install the missing dependencies? (Y or N):"))
        if self._answer.upper() == "Y" or self._answer.upper() == "YES":
            self._answerstate = True
        else:
            print "You need to install Java, NumPy, matplotlib, R and an appropriate compiler (e.g. gcc) before proceeding."
            print "Now quitting"
            quit()

        return

    def getanswerstate(self):
        """Return self._answerstate"""

        return self._answerstate


    def checkjava(self):
        """Check to see if Java is installed properly and included in the PATH"""

        try:
            subprocess.call(["java", "-version"], stderr=open(os.devnull, 'wb'))

        except:
            print "Couldn't find Java. It might not be installed or not included in the PATH"
            print "If it is not installed, please download and install the latest version of Java"
            print "It can be downloaded from http://www.java.com/en/"
            print "Trimming cannot be performed if Java is not installed. Quitting now."
            quit()


    def checkR(self):
        """Check to see if R is installed properly and included in the PATH"""

        try:
            subprocess.call(["R", "--version"], stdout=open(os.devnull, 'wb'))

        except:
            print "Couldn't find R. It might not be installed or not included in the PATH"
            print "If it is not installed, please install R from within the 'DE_analysis' folder within SPARTA"
            print "Differential gene expression cannot be performed if R is not installed. Quitting now."
            quit()
        return

    def checknumpy(self, options):
        """Check to see if NumPy module exists"""

        try:
            imp.find_module('numpy')
            self._foundnumpy = True

        except ImportError:
            self._foundnumpy = False

        if self._foundnumpy == False:
            print "You need to install 'NumPy' before proceeding, otherwise HTSeq will not work properly."
            if options.noninteractive:
                quit()
        return self._foundnumpy

    def checkhtseq(self):
        """Check to see if HTSeq module exists"""

        try:
            imp.find_module('HTSeq')
            self._foundhtseq = True

        except ImportError:
            self._foundhtseq = False

        if self._foundhtseq == False:
            print "HTSeq doesn't appear to be installed. An attempt will be made to install it for the local user."

        return self._foundhtseq

    # def checkmatplotlib(self):
    #     """Check to see if matplotlib module exists"""
    #
    #     try:
    #         imp.find_module('matplotlib')
    #         foundmatplotlib = True
    #
    #     except ImportError:
    #         foundmatplotlib = False
    #
    #     if foundmatplotlib == False:
    #         print "You need to install 'matplotlib' before proceeding, otherwise HTSeq will not work properly."
    #
    #     return foundmatplotlib

    def getNumPy(self, options):
        """Get latest NumPy iteration from sourceforge"""

        cd = CheckDependencies()
        subprocess.call(["curl", "-L", "-O", "http://sourceforge.net/projects/numpy/files/latest/download?source=files"])
        spartadir = cd.getSPARTAdir(options)
        #Make sure the user is in the SpartaDir
        subprocess.call(["mv", os.path.join(spartadir, "download?source=files"), os.path.join(spartadir, "numpy-latest.tar.gz")])
        subprocess.call(["tar", "-zxf", os.path.join(spartadir, "numpy-latest.tar.gz")])
        subprocess.call(["rm", os.path.join("numpy-latest.tar.gz")])
        return

    def installNumPy(self, options):
        """Install latest NumPy from source"""

        cd = CheckDependencies()
        current_numpy = glob.glob("numpy-*")[0]
        spartadir = cd.getSPARTAdir(options)
        os.chdir(spartadir + "/" + current_numpy)
        subprocess.call(["sudo python setup.py build install"], shell=True)
        return

    # def getmatplotlib(self):
    #     """Get latest NumPy iteration from sourceforge"""
    #     #IMPORTANT: for htseq-count, matplotlib is NOT required!
    #
    #     subprocess.call(["curl", "-L", "-O", "http://sourceforge.net/projects/matplotlib/files/latest/download?source=files"])
    #     subprocess.call(["mv download\?source\=files matplotlib-latest.tar.gz"], shell=True)
    #     subprocess.call(["tar", "-zxf", "matplotlib-latest.tar.gz"])
    #     subprocess.call(["rm", "matplotlib-latest.tar.gz"])
    #     return
    #
    # def installmatplotlib(self):
    #     """Install latest matplotlib from source"""
    #
    #     current_matplotlib = glob.glob("matplotlib-*")[0]
    #     spartadir = CheckDependencies.getSPARTAdir()
    #     os.chdir(spartadir + "/" + current_matplotlib)
    #     subprocess.call(["python setup.py build"], shell=True) #have not added 'install' yet because several dependencies
    #     #are required to install matplotlib...
    #     return


    def getpwd(self):
        """Get present working directory"""

        present_working_directory = subprocess.Popen("pwd", stdout=subprocess.PIPE).communicate()[0].strip("\n")
        return present_working_directory

    def getdesktoppath(self):
        """Get the path to the desktop"""

        desk_path = os.path.join(subprocess.Popen("echo $HOME", shell=True, stdout=subprocess.PIPE).stdout.readline().strip("\n"), "Desktop")
        return desk_path

    def getSPARTAdir(self, options):
        """Attempt to figure out where SPARTA is located. Default should be Desktop"""

        desk_path = os.path.join(subprocess.Popen("echo $HOME", shell=True, stdout=subprocess.PIPE).stdout.readline().strip("\n"), "Desktop")
        #This is explicitly coded to ensure that the rest of the functions are able to find the appropriate binaries
        try:
            if os.path.lexists(os.path.join(desk_path, "SPARTA_Linux")):
                sparta_dir = os.path.join(desk_path, "SPARTA_Linux")

            elif os.path.lexists(os.path.join(desk_path, "SPARTA_Linux-master")):
                sparta_dir = os.path.join(desk_path, "SPARTA_Linux-master")

        except:
            print "Couldn't find the SPARTA_Linux folder on the Desktop"

            if options.noninteractive:
                quit()

            while not os.path.lexists(sparta_dir):
                sparta_dir = str(raw_input("SPARTA_Linux folder is not on the Desktop. Please place the folder on the Desktop or enter the file path for the folder location or enter quit to exit the program:"))
                if sparta_dir.upper() == "Q" or sparta_dir.upper() == "QUIT":
                    quit()
                if not os.path.lexists(sparta_dir):
                    print("Invalid file path. The path you have selected does not exist or was not written correctly. \nAn example of path on Linux: /home/yourusername/Desktop/SPARTA_Linux")
        return sparta_dir

    def parseConfigFile(self, options):
        print "SPARTA is running in non-interactive mode."

        cd = CheckDependencies()
        spartadirloc = cd.getSPARTAdir(options)

        conditions_list = []
        #check and make sure that the ConfigFile.txt exists in the SPARTA directory
        if os.path.isfile(os.path.join(spartadirloc, "ConfigFile.txt")):
            with open(os.path.join(spartadirloc, "ConfigFile.txt"), "r") as configfile:
                with open(os.path.join(spartadirloc, 'conditions_input.txt'), "w") as conditions_input:
                    for line in configfile:
                        if line.startswith("Data") or line.startswith("Trimmomatic") or line.startswith("Bowtie") or line.startswith("HTSeq"):
                            parameters = line.split("->")[1]
                            paramlst = parameters.split(",")
                            if line.startswith("Data"):
                                dataloc = paramlst[0].strip()
                                inputname = paramlst[-1].strip()
                                data_path = os.path.join(subprocess.Popen("echo $HOME", shell=True, stdout=subprocess.PIPE).stdout.readline().strip("\n"), dataloc, inputname)
                            elif line.startswith("Trimmomatic"):
                                try:
                                    options.threads = paramlst[0].strip().split("=")[1]
                                    options.illuminaclip = paramlst[1].strip()[len("ILLUMINACLIP:"):]
                                    options.leading = paramlst[2].strip()[len("LEADING:"):]
                                    options.trailing = paramlst[3].strip()[len("TRAILING:"):]
                                    options.slidingwindow = paramlst[4].strip()[len("SLIDINGWINDOW:"):]
                                    options.minlentrim = paramlst[5].strip()[len("MINLEN:"):]
                                except Exception:
                                    print "There doesn't appear to be enough parameters associated with Trimmomatic."
                                    print "Make sure that there are 6: threads, ILLUMINACLIP, LEADING, TRAILING, SLIDINGWINDOW, MINLEN"
                                    print "Proceeding with defaults"
                            elif line.startswith("Bowtie"):
                                try:
                                    if paramlst[0].strip().split("=")[1] != '0' and paramlst[1].strip().split("=")[1] != "None":
                                        options.mismatch = paramlst[0].strip().split("=")[1]
                                        options.otherbowtieoptions = paramlst[1].strip().split("=")[1]
                                    elif paramlst[0].strip().split("=")[1] != '0':
                                        options.mismatch = paramlst[0].strip().split("=")[1]
                                    elif paramlst[1].strip().split("=")[1] != "None":
                                        options.otherbowtieoptions = paramlst[1].strip().split("=")[1]
                                    else:
                                        print "There doesn't appear to be any options specified for Bowtie."
                                        print "Proceeding with defaults"
                                except Exception:
                                    print "There seems to be an error with the options input for Bowtie."
                                    print "Proceeding with defaults"
                            elif line.startswith("HTSeq"):
                                try:
                                    options.stranded = paramlst[0].strip().split("=")[1]
                                    options.order = paramlst[1].strip().split("=")[1]
                                    options.minqual = paramlst[2].strip().split("=")[1]
                                    options.type = paramlst[3].strip().split("=")[1]
                                    options.idattr = paramlst[4].strip().split("=")[1]
                                    options.mode = paramlst[5].strip().split("=")[1]
                                except Exception:
                                    print "There appears to be an issue with the parameter input for HTSeq"
                                    print "Proceeding with defaults"

                        elif line.upper().startswith("REFERENCE") or line.upper().startswith("EXPERIMENTAL"):
                            #Generate the conditions_input.txt file
                            conditions_input.write(line)
                            #Strip off any newlines
                            conditions = line.strip('\n')
                            #Split on the colon and grab just the conditions
                            colsep = conditions.split(':')[1]
                            #Clean up the condition input by attempting to remove spaces and tabs if present
                            condition = colsep.split(',')
                            conditions_list.append(condition)


        return conditions_list, data_path
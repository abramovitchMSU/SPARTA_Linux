__author__ = 'benkjohnson'

import os
import subprocess
import numpy as np
import re

class DifferentialExpression(object):
    def __init__(self):

        return

    def removenoncountdata(self, analysislocation):
        '''Remove the last 5 lines of the HTSeq output
        This method works for the small text file output
        and is platform independent, but probably would be
        *very* slow for large files.'''

        countspath = os.path.join(analysislocation, 'HTSeq')
        depath = os.path.join(analysislocation, 'DEanalysis')
        countdata = os.listdir(countspath)
        for countfile in countdata:
            extension = os.path.splitext(countfile)[1]
            if extension == ".sam":
                with open(os.path.join(countspath, countfile), "r") as counts:
                    line = counts.readlines()
                with open(os.path.join(depath, countfile), "w") as newcounts:
                    for newcountsline in line:
                        if not newcountsline.startswith("__"):
                            newcounts.write(newcountsline)
        return

    def getuserinput(self, analysislocation):
        '''Take in user input for specifying the conditions to be tested'''

        depath = os.path.join(analysislocation, 'DEanalysis')
        print "SPARTA has these files:"
        number = 1
        for countsfile in os.listdir(depath):
            print "{num}) {file}".format(num=number, file=countsfile)
            number += 1
        moveforward = False
        while moveforward == False:
            conditionnumber = str(raw_input("How many conditions are there?:"))
            if conditionnumber in ['', '0', '1']:
                print "You can't compare a condition to itself."
            else:
                try:
                    int(conditionnumber)
                    usersure = str(raw_input("Are you sure that's how many conditions you would like to compare? (y/n):"))
                    if usersure.upper() == "Y":
                        moveforward = True
                    else:
                        continue
                except Exception:
                    print "Please enter a number (i.e. 2)"
        conditionnumber = int(conditionnumber)
        conditioncounter = 1
        with open(os.path.join(analysislocation, 'DEanalysis', 'conditions_input.txt'), "w") as conditions_input:
            while conditioncounter != conditionnumber+1:
                if conditioncounter == 1:
                    conditions_input.write("Reference_Condition_Files:\n")
                    conditioncounter += 1
                else:
                    outline = "Experimental_Condition_{val}_Files:\n".format(val=conditioncounter)
                    conditions_input.write(outline)
                    conditioncounter += 1

        exampleconditionnumber = 2
        exampleconditioncounter = 1
        with open(os.path.join(analysislocation, 'DEanalysis', 'conditions_input_example.txt'), "w") as exampleconditions_input:
            while exampleconditioncounter != exampleconditionnumber+1:
                if exampleconditioncounter == 1:
                    exampleconditions_input.write("Reference_Condition_Files: mapWT_A.sam, mapWT_B.sam, mapWT_C.sam\n")
                    exampleconditioncounter += 1
                else:
                    exampleoutline = "Experimental_Condition_{val}_Files: mapmutant2466_A.sam, mapmutant2466_B.sam, mapmutant2466_C.sam\n".format(val=exampleconditioncounter)
                    exampleconditions_input.write(exampleoutline)
                    exampleconditioncounter += 1

        print "Now we need to edit a text file to specify which files belong to a given condition and which files are replicates for each condition."
        print "The file names that you need to use are listed above under 'SPARTA has these files'."
        print "The text file you need to edit (NOT WITH MICROSOFT WORD) using a text editor like gedit, is on your Desktop in RNAseq_Data -> date of the current run -> DEanalysis -> conditions_input.txt"
        print "Enter the relevant file names, with replicates separated by a comma."
        print "As an example, please see the 'conditions_input_example.txt' in the DEanalysis folder."

        moveforward = False
        conditions_list = []

        #I want to make sure that this next section catches as many potential errors as possible without killing the program,
        #because at this point, it may have taken a rather long time to get here. Also, asking for user input can be trouble.
        #Attempt to catch and account for as many user errors as possible.
        while moveforward == False:
            proceedanswer = str(raw_input("Once you have entered the file names, hit Enter/Return:"))
            if proceedanswer == '':
                with open(os.path.join(analysislocation, 'DEanalysis', 'conditions_input.txt'), "r") as user_condition_input:
                    line = user_condition_input.readlines()
                    for conditions in line:
                        if conditions != '\n':
                            try:
                                conditions = conditions.strip('\n')
                                colsep = conditions.split(':')[1]
                                linesep = colsep.split(',')
                                reptry = linesep[1]
                                for userfile in linesep:
                                    if userfile != '':
                                        userfile = userfile.strip()
                                        os.path.lexists(os.path.join(analysislocation, 'DEanalysis', userfile)) == True
                                moveforward = True
                            except Exception:
                                print "Make sure that the replicates are separated with commas, start with 'map', and end with '.sam'.\n If you would like an example, please refer to conditions_input_example.txt in the DEanalysis folder.\n Please edit the conditions_input.txt and try again.\n"
                                print "Make sure you spelled the file name correctly as well! The names of the files must match exactly as the ones in the DEanalysis folder."
                                print "If you don't have a replicate for a given condition, add a comma after the name of the file."
                                print "As an example: Experimental_Condition_2_Files: mapmutant2466_A.sam,"
                                print "You want to proceed with caution though as the condition is not replicated."
                                break
                        #Strip off any newlines
                        conditions = conditions.strip('\n')
                        #Split on the colon and grab just the conditions
                        colsep = conditions.split(':')[1]
                        #Clean up the condition input by attempting to remove spaces and tabs if present
                        condition = colsep.split(',')
                        conditions_list.append(condition)

        return conditions_list

    def generatecontrasts(self, contrastrow):

        if contrastrow > 2:
            contrastcounter = 2
            contrastcol = 0
            while contrastcounter != contrastrow:
                contrastcol += contrastrow - contrastcounter
                contrastcounter += 1
        if contrastcounter != 0:
            contrast = np.zeros((contrastcol,contrastrow)).astype(int)
            contrastcounter = 2
            rowindex = 0
            while contrastcounter != contrastrow:
                iterationnumber = contrastrow - contrastcounter
                itercounter = 0
                while itercounter != iterationnumber:
                    contrast[rowindex][contrastcounter - 1] = 1
                    contrast[rowindex][(contrastcounter - 1)+(itercounter + 1)] = -1
                    itercounter += 1
                    rowindex += 1
                contrastcounter += 1

        return contrast

    def writeRscript(self, analysislocation, conditions_list):
        '''Write the R script, incorporating the user input'''

        #Generate the contrasts for each comparison beyond experimental to reference
        #This is a rather kludgy way to make this work, but it works, and scales nicely.
        contrastrow = len(conditions_list)
        if contrastrow > 2:
            de = DifferentialExpression()
            contrast = de.generatecontrasts(contrastrow)
        else:
            contrast = None

        #Writing this all is going to be rather difficult...

        with open(os.path.join(analysislocation, 'DEanalysis', 'DEexpression.r'), "w") as de_expression:
            #de_expression.write("source('http://bioconductor.org/biocLite.R')\n")
            #de_expression.write("biocLite('edgeR')\n")
            de_expression.write("library('edgeR')\n")

            Rvarlst = []
            for treatments in conditions_list:
                repnumcounter = 1
                for varnames in treatments:
                    if varnames != '':
                        varnames = varnames.strip().strip('\t')
                        varname = os.path.splitext(varnames)[0]
                        varname = varname[len('map'):]
                        filepathname = os.path.join(analysislocation, 'DEanalysis', varnames)
                        de_expression.write("{varname}_{repnum} <- read.table('{filepathname}', row.names=1)\n".format(varname=varname, repnum=repnumcounter, filepathname=filepathname))
                        de_expression.write("colnames({varname}_{repnum}) <- '{varname}_{repnum}'\n".format(varname=varname, repnum=repnumcounter))
                        Rvarlst.append("{varname}_{repnum}".format(varname=varname, repnum=repnumcounter))
                        repnumcounter += 1
            Rvarlst = map(str, Rvarlst)
            Rvarlstjoin = ",".join(Rvarlst)
            de_expression.write("alldata <- cbind(" + Rvarlstjoin + ")\n")

            Rgrouplengths = []
            Rgrouplst = []
            repnum = 1
            for treatments in range(len(conditions_list)):
                for group in range(len(conditions_list[treatments])):
                    if treatments != len(conditions_list)+1 and conditions_list[treatments][group] != '':
                        Rgrouplst.append(repnum)
                repnum += 1
                Rgrouplengths.append(len(Rgrouplst)) #eventually use this to try and test for batch effects in n-conditions
            Rgrouplst = map(str, Rgrouplst)
            Rgrouplstjoin = ",".join(Rgrouplst)
            de_expression.write("group <- factor(c(" + Rgrouplstjoin + "))\n")

            de_expression.write("y <- DGEList(counts=alldata, group=group)\n")

            halfofsamples = len(Rvarlst)/2
            de_expression.write("keep <- rowSums(cpm(y)>2)>={halfofsamples}\n".format(halfofsamples=halfofsamples))
            de_expression.write("y <- y[keep,]\n")

            de_expression.write("y$samples$lib.size <- colSums(y$counts)\n")

            de_expression.write("y <- calcNormFactors(y)\n")

            # if MDSmethod == False:
            #     de_expression.write("y <- calcNormFactors(y)\n")
            # else:
            #     de_expression.write("y <- calcNormFactors(y, method=BCV)\n")
            #     print "Using BCV method"
            # print "If there is large variation between library sizes, you may want to re-run the DE analysis with 'MDSmethod=BCV'"

            de_expression.write("y$samples\n")

            de_expression.write("png('" + os.path.join(analysislocation, 'DEanalysis', 'MDSplot.png') + "')\n")
            de_expression.write("mymdsobj <- plotMDS(y)\n")
            de_expression.write("dev.off()\n")

            de_expression.write("design <- model.matrix(~group)\n")

            de_expression.write("y <- estimateGLMCommonDisp(y, design, verbose=TRUE)\n")
            de_expression.write("y <- estimateGLMTrendedDisp(y, design)\n")
            de_expression.write("y <- estimateGLMTagwiseDisp(y, design)\n")

            de_expression.write("png('" + os.path.join(analysislocation, 'DEanalysis', 'BCVplot.png') + "')\n")
            de_expression.write("plotBCV(y)\n")
            de_expression.write("dev.off()\n")

            de_expression.write("fit <- glmFit(y, design)\n")
            de_expression.write("lrt <- glmLRT(fit, coef=2)\n")

            if contrast is None:
                numofcond = 2
            else:
                numofcond = len(contrast[0])
            condval = 2
            while condval != numofcond+1:
                de_expression.write("lrt <- glmLRT(fit, coef={condval})\n".format(condval=condval))
                de_expression.write("FDR <- p.adjust(lrt$table$PValue, method='BH')\n")
                de_expression.write("summary(dt <- decideTestsDGE(lrt))\n")
                de_expression.write("isDE <- as.logical(dt)\n")
                de_expression.write("DEnames <- rownames(y)[isDE]\n")
                de_expression.write("png('" + os.path.join(analysislocation, 'DEanalysis', 'ReferenceCondvsExpCond{condval}_scatter.png'.format(condval=condval)) + "')\n")
                de_expression.write("plotSmear(lrt, de.tags=DEnames, cex=0.5)\n")
                de_expression.write("abline(h=c(-1,1), col='blue')\n")
                de_expression.write("dev.off()\n")
                de_expression.write("outfile <- cbind(cpm(y), lrt$table, FDR)\n")
                de_expression.write("write.csv(outfile, file='" + os.path.join(analysislocation, 'DEanalysis', 'ReferenceCondvsExpCond{condval}_DE.csv'.format(condval=condval)) + "')\n")
                condval += 1

            if contrast is not None:
                for contr in contrast:
                    trts = np.nonzero(contr)
                    reftrt = trts[0][0] + 1
                    exptrt = trts[0][1] + 1
                    contr = tuple(contr)
                    de_expression.write("lrt <- glmLRT(fit, contrast=c{contr})\n".format(contr=contr))
                    de_expression.write("FDR <- p.adjust(lrt$table$PValue, method='BH')\n")
                    de_expression.write("summary(dt <- decideTestsDGE(lrt))\n")
                    de_expression.write("isDE <- as.logical(dt)\n")
                    de_expression.write("DEnames <- rownames(y)[isDE]\n")
                    de_expression.write("png('" + os.path.join(analysislocation, 'DEanalysis', 'ExpCond{reftrt}vsExpCond{exptrt}_scatter.png'.format(reftrt=reftrt, exptrt=exptrt)) + "')\n")
                    de_expression.write("plotSmear(lrt, de.tags=DEnames, cex=0.5)\n")
                    de_expression.write("abline(h=c(-1,1), col='blue')\n")
                    de_expression.write("dev.off()\n")
                    de_expression.write("outfile <- cbind(cpm(y), lrt$table, FDR)\n")
                    de_expression.write("write.csv(outfile, file='" + os.path.join(analysislocation, 'DEanalysis', 'ExpCond{reftrt}vsExpCond{exptrt}_DE.csv'.format(reftrt=reftrt, exptrt=exptrt)) + "')\n")
            #Batch effect test: compare the first two samples of the conditions tested based on MDS plot
            #Reasoning: if the first two samples
            grouplen = len(Rgrouplengths)
            batchcounter = 1
            while batchcounter != grouplen:
                batcheffectTestIdxA, batcheffectTestIdxB = (Rgrouplengths[0] - Rgrouplengths[0]) + 1, (Rgrouplengths[batchcounter] - Rgrouplengths[0]) + 1
                de_expression.write("ifelse((mymdsobj$cmdscale.out[{batcheffectTestIdxA}, {batcheffectTestIdxA}] > 0) == (mymdsobj$cmdscale.out[{batcheffectTestIdxB}, {batcheffectTestIdxA}] > 0), ifelse((mymdsobj$cmdscale.out[{batcheffectTestIdxA}, 1] > 0) != (mymdsobj$cmdscale.out[{batcheffectTestIdxB}, 1] > 0), 'Reference vs. Experimental {batchcounter}: Potential batch effect', 'Reference vs. Experimental {batchcounter}: All appears to be well'), 'Reference vs. Experimental {batchcounter}: All appears to be well')\n".format(batcheffectTestIdxA=batcheffectTestIdxA, batcheffectTestIdxB=batcheffectTestIdxB, batchcounter=batchcounter + 1))
                batchcounter += 1
            #End of batch effect test

        return


    def runRscript(self, analysislocation):
        '''Execute the pre-written R script'''

        #Read in the pre-written R script from the DEanalysis folder
        Rscriptloc = os.path.join(analysislocation, 'DEanalysis', 'DEexpression.r')

        #Run the script
        #You have to run as super user to install edgeR
        subprocess.Popen("R --vanilla --slave < " + Rscriptloc, shell=True).wait()

        #Send a final notice to the user
        print "Analysis complete. Thank you for using SPARTA."

        return

    def de_analysis(self, analysislocation):
        '''Call all of the DE functions and run them'''

        de = DifferentialExpression()
        de.removenoncountdata(analysislocation)
        conditions_list = de.getuserinput(analysislocation)
        de.writeRscript(analysislocation, conditions_list)
        de.runRscript(analysislocation)

        return

    def de_analysis_noninteractive(self, analysislocation, conditions_list):
        de = DifferentialExpression()
        de.removenoncountdata(analysislocation)
        de.writeRscript(analysislocation, conditions_list)
        de.runRscript(analysislocation)

        return
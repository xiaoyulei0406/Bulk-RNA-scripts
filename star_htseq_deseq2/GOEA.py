#!/usr/bin/env python
from __future__ import print_function
import sys
import gzip
import os
import re
import csv
from operator import itemgetter
from openpyxl import load_workbook
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.colors as mcolors
import pandas as pd
import time
import datetime


from collections import namedtuple
setting = namedtuple('setting', ['geneSCF', 'goobo', 'ensembl'])

def KEGGenrichment(inputFile,ensemToGene,outputPrefix,inputOrg,totalGeneN,barN,excels):
    inputHandle=open(inputFile,"rU")
    #basenameInput=os.path.basename(inputFile)
    outputFileGene=outputPrefix+".geneID.list"
    output=open(outputFileGene,"w")
    for inline in inputHandle:
        infor=inline.strip()
        if ensemToGene.has_key(infor):
            print (ensemToGene[infor],file=output)
    inputHandle.close()
    output.close()
    outputFolder="/".join(outputFileGene.split("/")[:-1])
    if inputOrg:
        cmd=setting.geneSCF+" -m=normal -t=gid -db=KEGG -p=no -i="+outputFileGene+" -o="+outputFolder+" -bg="+str(totalGeneN)+" -org="+inputOrg
        print (cmd)
        import subprocess
        subprocess.call([cmd],shell=True)
        #print (outputFileGene)
        os.rename(outputFileGene+"_KEGG_"+inputOrg+"_functional_classification.tsv",outputPrefix+"-kegg")
        (geneList,geneCount,pvalue,excels)=processGENEscf(outputPrefix+"-kegg",excels)
        if geneList:
            Boxplot(geneList[:barN],geneCount[:barN],pvalue[:barN],outputPrefix+"-keggEnrichment.png","Top "+str(barN)+" KEGG enrichment")
        return(excels)

def processGENEscf(inputTSV,excels):
    inputF=csv.reader(open(inputTSV), delimiter="\t")
    headers=next(inputF,None)
    sortList=sorted(inputF, key = lambda x: float(x[5]))
    outputFile=open(inputTSV+"significant","w")
    sortFile=open(inputTSV+".sort","w")
    geneList=[]
    geneCount=[]
    pvalue=[]
    print ("\t".join(headers),file=sortFile)
    print ("\t".join(headers),file=outputFile)
    for infor in sortList:
        print ("\t".join(infor),file=sortFile)
        if float(infor[5])<0.05:
            print ("\t".join(infor),file=outputFile)
            geneList.append(infor[1])
            geneCount.append(int(infor[2]))
            pvalue.append(float(infor[5]))
    outputFile.close()
    sortFile.close()
    os.rename(inputTSV+".sort",inputTSV)
    excels.append(inputTSV+".xlsx")
    excels.append(inputTSV+"significant.xlsx")
    df_new = pd.read_csv(inputTSV, sep=r'\t', engine='python')
    writer = pd.ExcelWriter(inputTSV+".xlsx", engine='xlsxwriter')
    df_new.to_excel(writer, index = True)
    writer.save()
    df_sig =pd.read_csv(inputTSV+"significant", sep=r'\t', engine='python')
    writerSig = pd.ExcelWriter(inputTSV+"significant.xlsx", engine='xlsxwriter')
    df_sig.to_excel(writerSig, index = False)
    writerSig.save()
    os.remove(inputTSV)
    os.remove(inputTSV+"significant")
    return(geneList,geneCount,pvalue,excels)

def keggOrgMapper():
    keggOrg=dict()
    keggOrgList=set()
    orgFile="/".join(setting.geneSCF.strip().split("/")[:-1])+"/org_codes_help/KEGG_organism_codes.txt"
    orgFileHandle=open(orgFile,"rU")
    for inline in orgFileHandle:
        if not inline.startswith("#"):
            infor=inline.strip().split("\t")
            g=re.search("\((.+)\)",inline)
            if g:
                keggOrg[g.group(1).lower()]=infor[1]
                keggOrgList.add(infor[1])
                #print (g.group(1).lower()+" "+infor[1])
    orgFileHandle.close()
    return (keggOrg,keggOrgList)


def Boxplot(funTerm,funCount,funPvalue,outputFile,PlotTitle):
    plt.switch_backend('agg')
    plt.rcdefaults()
    fig, ax = plt.subplots()
    y_pos= np.arange(len(funTerm))
    normalize = mcolors.Normalize(vmin=min(funPvalue), vmax=max(funPvalue))
    colormap = cm.jet
    newcolor=[]
    ylabMaxLen=len(max(funTerm, key=len))
    for n in funPvalue:
        newcolor.append(colormap(normalize(n)))
    ax.barh(y_pos, funCount,align='center',color=newcolor)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(funTerm)
    ax.invert_yaxis()
    ax.set_title(PlotTitle)
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(funPvalue)
    cbar = plt.colorbar(scalarmappaple)
    cbar.set_label("p-value")
    fig.savefig(outputFile,bbox_inches = "tight")
    plt.close(fig)

    
def loadGOdb():
    obo_fname = download_go_basic_obo(obo=setting.goobo)
    obodag = GODag(setting.goobo)
    return (obodag)

def loadDB(inputOrg):
    db=inputOrg
    """
    if inputOrg=="human":
        db=setting.humanEnsemble
    elif inputOrg=="mouse":
        db=setting.mouseEnsemble
    elif inputOrg=="pig":
        db=setting.pigEnsemble
    elif inputOrg=="cow":
        db=setting.cowEnsemble
    """
    db=setting.ensembl
    assoc = {}
    allGene=set()
    ensTogene={}
    with gzip.open(db, 'rt') as f:
        for inline in f:
            if (not inline.startswith("Gene")) and (not inline.startswith("#")):
                infor=inline.strip().split()
                #if len(infor)<=3:
                allGene.add(infor[0])
                #if len(infor)==3:
                for index,eachItem in enumerate(infor):
                    if index>0:
                        if eachItem.startswith("GO:"):
                            if not assoc.has_key(infor[0]):
                                assoc[infor[0]]=set()
                            assoc[infor[0]].add(eachItem)
                        else:
                            if not ensTogene.has_key(infor[0]):
                                ensTogene[infor[0]]=eachItem
    return (allGene, assoc,ensTogene)
    


def goea(inputFile,methods,bkGene,assocDict,obodag,outputFile,excels,depth,barN):
    methodsList=["bonferroni", "sidak", "holm", "fdr"]
    if methods=="all":
        Inputmethods =methodsList
    elif "," in methods:
        Inputmethods=methods.split(",")
        if not (all([z in methodsList for z in Inputmethods])):
            print ("The chosen method is not N/A")
            sys.exit()
    elif methods in methodsList:
        Inputmethods=[methods]
    else:
        print ("The chosen method is not N/A")
        sys.exit()
    goeaobj = GOEnrichmentStudy(
        bkGene, # List of mouse protein-coding genes
        assocDict, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = True,
        alpha = 0.05, # default significance cut-off
        methods = Inputmethods)
    study=[]
    ana=inputFile
    if os.path.isfile(ana):
        inputF=open(ana,"rU")
        for inline in inputF:
            study.append(inline.strip())
        inputF.close()
    if study and assocDict and bkGene:
        goea_results_all = goeaobj.run_study(study)
        goeaobj.wr_xlsx(outputFile+"-allGO.xlsx", goea_results_all)
        excels.append(outputFile+"-allGO.xlsx")

        sig=[r for r in goea_results_all if r.enrichment=="e" and r.p_uncorrected<0.05]
        goeaobj.wr_xlsx(outputFile+"-GOsignificant.xlsx",sig)
        excels.append(outputFile+"-GOsignificant.xlsx")
        
        filterBP=[r for r in goea_results_all if r.enrichment=="e" and r.NS=="BP" and r.depth>=depth and r.p_uncorrected<0.05]
        goeaobj.wr_xlsx(outputFile+"-GOBiologicalProcess.xlsx",filterBP)
        fileName=outputFile+"-GOBiologicalProcess.png"
        plotTitle="GO Biological Process"
        plotGO(filterBP[:barN],fileName,plotTitle)
        excels.append(outputFile+"-GOBiologicalProcess.xlsx")

        filterCC=[r for r in goea_results_all if r.enrichment=="e" and r.NS=="CC" and r.depth>=depth and r.p_uncorrected<0.05]
        goeaobj.wr_xlsx(outputFile+"-GOCellularComponent.xlsx",filterCC)
        fileName=outputFile+"-GOCellularComponent.png"
        plotTitle="GO Cellular Component"
        plotGO(filterCC[:barN],fileName,plotTitle)
        excels.append(outputFile+"-GOCellularComponent.xlsx")
        filterMF=[r for r in goea_results_all if r.enrichment=="e" and r.NS=="MF" and r.depth>=depth and r.p_uncorrected<0.05]
        fileName=outputFile+"-GOMolecularFunction.png"
        plotTitle="GO Molecular Function"
        goeaobj.wr_xlsx(outputFile+"-GOMolecularFunction.xlsx",filterMF)
        plotGO(filterMF[:barN],fileName,plotTitle)
        excels.append(outputFile+"-GOMolecularFunction.xlsx")
    return(excels)

def combineExc(excelList,outputFile):
    writer=pd.ExcelWriter(outputFile,engine='xlsxwriter')
    for i in excelList:
        df1=pd.read_excel(i,"Sheet1")
        sheetName1=os.path.basename(i).replace(".xlsx","").split("-")[-1]
        os.remove(i)
        df1.to_excel(writer,sheet_name=sheetName1[:25])
    writer.save()

    
def plotGO(filterGO,fileName,plotTitle):
    if filterGO:
        FunTerm=[]
        FunCount=[]
        FunPval=[]
        for each in filterGO:
            FunTerm.append(each.name)
            FunCount.append(each.study_count)
            FunPval.append(each.p_uncorrected)
        if FunTerm:
            Boxplot(FunTerm,FunCount,FunPval,fileName,plotTitle)

def main():
    inputFile=None
    methods = "bonferroni"
    org=None
    bk=None
    outputFile=None
    depth=2
    prefix="result"
    barN=15
    excels=[]
    Keggorg=None
    kegg="on"
    go="on"
    logFile=None
    geneSCF=None
    goobo=None
    ensembl=None
    #methods = ["bonferroni", "sidak", "holm", "fdr"]
    usage="python GOEA.py -in <ensemble gene list> -org <organism: human;mouse;pig;cow.> -out <output folder>\n\n"
    usage+="parameters:\n"
    usage+='%-10s' % ("-in:")+"input file contains one column of ensemble gene list\n"
    usage+='%-10s' % ("-org:")+"the pre-defined database that include human;mouse;pig;cow.\n"
    usage+='%-10s' % ("")+"Other org anism can be define in self-defined database \"-db\".\n"
    usage+='%-10s' % ("")+"If the pre-defined org is set, the self-defined db will not be considered.\n"
    usage+='%-10s' % ("-out:")+"output folder \n\n"
    usage+="parameters[optional]:\n"
    usage+='%-10s' % ("-prefix:")+"output file prefix [Default:result]\n"
    usage+='%-10s' % ("-db:")+"the self-defined database. Check the format in README (For GO)\n"
    usage+='%-10s' % ("-keggOrg:")+"the self-defined organism (For KEGG)\n"
    usage+='%-10s' % ("")+"The organism code must be in KEGG org list\n"
    usage+='%-10s' % ("-method:")+"the methods chosen from [bonferroni,sidak,holm,fdr] (For GO).\n"
    usage+='%-10s' % ("")+"For example: bonferroni,sidak or all. [Default: bonferroni]\n"
    usage+='%-10s' % ("-kegg:")+"run kegg enrichment pathway (\"on\") or not (\"off\").  (Default: on)\n"
    usage+='%-10s' % ("-go:")+"run GO enrichment pathway (\"on\") or not (\"off\").  (Default: on)\n\n"
    usage+="paramter[optional for plot]:\n"
    usage+='%-10s' % ("-barN:")+"show top identified GO categories/KEGG pathways [Default:15]\n"
    usage+='%-10s' % ("-depth:")+"report the specific depth of GO categories [Default:2]\n"
    usage+='%-10s' % ("-geneSCF:")+"report the specific depth of GO categories [Default:2]\n"
    usage+='%-10s' % ("-goobo:")+"report the specific depth of GO categories [Default:2]\n"
    usage+='%-10s' % ("-ensembl:")+"report the specific depth of GO categories [Default:2]\n"
    
    for idx in range(len(sys.argv)):
        if (sys.argv[idx] == "-in") and (len(sys.argv) > idx + 1):
            inputFile=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-h") or (sys.argv[idx] == "-help"):
            print (usage)
            sys.exit()
        elif (sys.argv[idx] == "-org") and (len(sys.argv) > idx + 1):
            org=(sys.argv[idx + 1]).lower()
        elif (sys.argv[idx] == "-method") and (len(sys.argv) > idx + 1):
            methods=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-db") and (len(sys.argv) > idx + 1):
            bk=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-out") and (len(sys.argv) > idx + 1):
            outputFile=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-depth") and (len(sys.argv) > idx + 1):
            depth=int(sys.argv[idx + 1])
        elif (sys.argv[idx] == "-prefix") and (len(sys.argv) > idx + 1):
            prefix=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-barN") and (len(sys.argv) > idx + 1):
            barN=int(sys.argv[idx + 1])
        elif (sys.argv[idx] == "-keggOrg") and (len(sys.argv) > idx + 1):
            Keggorg=(sys.argv[idx + 1]).lower()
        elif (sys.argv[idx] == "-kegg") and (len(sys.argv) > idx + 1):
            kegg=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-go") and (len(sys.argv) > idx + 1):
            go=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-geneSCF") and (len(sys.argv) > idx + 1):
            geneSCF=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-goobo") and (len(sys.argv) > idx + 1):
            goobo=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-ensembl") and (len(sys.argv) > idx + 1):
            ensembl=sys.argv[idx + 1]
    flag=0
    if (not inputFile) or (not outputFile):
        flag=1
    if org:
        if (not geneSCF or not goobo or not ensembl):
            flag = 1
        else:
            setting.geneSCF = geneSCF
            setting.goobo = goobo
            setting.ensembl = ensembl
        if (not org=="human") and (not org=="cow") and (not org=="mouse") and (not org=="pig"):
            print ("Error: the selected org is NA!")
            flag=1
        else:
            (allGene, assoc,ensToGene)=loadDB(org)
    else:
        if (not bk):
            flag=1
        else:
             (allGene, assoc,ensToGene)=loadDB(bk)
             
    if flag==1:
        print (usage)
        sys.exit()
    else:
        outputFolder=outputFile
        outputFile=outputFile+"/"+prefix
        logFile=outputFile+".log"
        loghandle=open(logFile,"w")
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print ("Starts at ",end=" ")
        print(st, file=loghandle)
        print ("#command line: "+" ".join(sys.argv),file=loghandle,end="\n")
        loghandle.flush()
        GOdb=loadGOdb()
        if go=="on":
            print("Run go enrichment analysis ...", file=loghandle)
            print("Run go enrichment analysis ...")
            loghandle.flush()
            excels=goea(inputFile,methods,allGene,assoc,GOdb,outputFile,excels,depth,barN)
        if kegg=="on":
            print("Run KEGG enrichment analysis ...", file=loghandle)
            print("Run KEGG enrichment analysis ...")
            loghandle.flush()
            keggOrgInput=None
            print (Keggorg)
            (keggOrgDict,keggOrgList)=keggOrgMapper()
            if keggOrgDict.has_key(org):
                keggOrgInput=keggOrgDict[org]
            elif Keggorg in keggOrgList:
                keggOrgInput=Keggorg
            else:
                print ("Error: the organism information you provided cannot be found in KEGG org list!",file=loghandle)
                #sys.exit()

            excels=KEGGenrichment(inputFile,ensToGene,outputFile,keggOrgInput,len(allGene),barN,excels)
        if excels:
            combineExc(excels,outputFile+"-AnalysisResult.xlsx")
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print ("Finish at ",file=loghandle,end=" ")
        print(st, file=loghandle)

        loghandle.close()
if __name__ == '__main__':
    main()

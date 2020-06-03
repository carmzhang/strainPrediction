import os.path
import xlwt
import os

from optparse import OptionParser
import subprocess
import sys
import difflib
import xlsxwriter
import re
import math
path="D:\\Duke\\Duke_201904\\CARD_results"
path3="D:\\Duke\\Duke_201904\\CARD"

def getwork():
    files=os.listdir(path)
    wb=xlsxwriter.Workbook((path3+"\\"+'Card_results20190425_640.xlsx'))
    ws1 = wb.add_worksheet('drug class')
    ws1.write (0,0,'ID')
    dic_tmp={}
    dic_drug={}
    dic_results={}
    idlist={}
    results=[]
    bt_results=[]
    druglist=[]
    btlist=[]
    number=1
    nb=1
    nbbt=1
    druglist=["cephalosporin","diaminopyrimidine antibiotic","fluoroquinolone antibiotic","acridine dye",'aminoglycoside antibiotic','penam','tetracycline antibiotic','phenicol antibiotic','nucleoside antibiotic','fusidic acid','macrolide antibiotic','lincosamide antibiotic','streptogramin antibiotic','pleuromutilin antibiotic','monobactam','carbapenem','cephamycin','glycylcycline','peptide antibiotic','ADC beta-lactamase','OXA beta-lactamase','Bc beta-lactamase']
    for drug in druglist:
        ws1.write(0,nb,drug)
        dic_drug[drug]=nb
        nb += 1
    for i in range(len(files)):
        workfile = path+"\\"+files[i]
        file=files[i]
        file_name=file.split('.')
        infile = open(workfile, "r").readlines()
        idlist[file_name[0]]=number
        ws1.write(number,0,file_name[0])
        number += 1    
        for m in druglist:
            dic_tmp[m]=[]
            for n in range(1,len(infile)):
                line=infile[n].strip().split("\t")
                drugl=line[14].split('; ')
                btdata=line[16].split(' ')
                gene=line[8]
                identity=line[9]
                coverage=line[20]
                prolength=len(line[19])
                for dr in drugl:
                    if dr == m:
                        result=[gene,str(identity),str(coverage),str(prolength),'|']
                        dic_tmp[dr].extend(result)
                        result=[]
                    else:
                        continue
                for mn in range(0,len(btdata)):
                    if btdata[mn] == 'beta-lactamase':
                        bt=str(btdata[mn-1])+' beta-lactamase'
                        if m==bt:
                            bt_result=[gene,str(identity),str(coverage),str(prolength),'|']
                            dic_tmp[bt].extend(bt_result)
                            bt_result=[]
                    else:
                        continue
        for ider in dic_tmp.keys():
            tmp=(' ').join(dic_tmp[ider])
            ws1.write(idlist[file_name[0]],dic_drug[ider],tmp)                        
    wb.close()
        
if __name__=='__main__':
    getwork()

       
        


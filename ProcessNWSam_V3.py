#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
import time
import gzip
import itertools
import sys
import json
import pysam
from Levenshtein import distance

samp=sys.argv[1]
RE=sys.argv[2]
prbsALL = json.load(open(sys.argv[3]))
RefGenome = sys.argv[4]
inputDir = sys.argv[5]
outputDirSam =sys.argv[6]
outputDirStats = sys.argv[7]

useUMI=False


#if len(sys.argv) > 5: 
 #  useUMI=sys.argv[5] == "Y"


#$reffa $inputDir $outputDirStats

print ("Parameters")
print (sys.argv)



#samp="NW009"
#RE="AAGCTT"
#useUMI=True

#samp="NW047.NW28HRPB009tech22_MID06"
#RE="AAGCTT"
#useUMI=False


def distancex(seq1, seq2):
    oneago = None
    thisrow = range(1, len(seq2) + 1) + [0]
    for x in xrange(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]



def parseCigar(cigar):
    chars=["M","N","D","I","S","H","P","X","="]    
    cts={}
    temp=""
    x=0
    for i in range(0,len(cigar)): 
        if not cigar[i] in chars: 
           temp +=  cigar[i]
        else:            
           cts[len(cts)] = [cigar[i],int(temp)]
           if cigar[i] in cts: 
              cts[cigar[i]] += int(temp)
           else: 
              cts[cigar[i]] = int(temp)  
           temp=""
    return cts


############################################################################################################################
## LOAD REQUIRED DATASETS
############################################################################################################################

refGenome={}
for chromo in SeqIO.parse(open(RefGenome,"r"), "fasta"):
    refGenome[chromo.id] = str(chromo.seq).upper()
    print(chromo.id)


lstBlackListed = {}
for x in open("ENCODE_Blacklist.regions.tsv"):
    [SRC,chrom,chromStart,chromEnd,name] = x[:-1].split("\t")
    if SRC != 'SRC': 
       a = int((int(chromStart)/1000) - 1)
       b = int((int(chromEnd)/1000) + 2 )
        
       for i in range(a,b): 
           lstBlackListed[chrom+":"+ str(i)] = True

   



probes={}
for xprb in prbsALL: 
    prb = prbsALL[xprb]
    szPrb = prb[4] - prb[3]
    if prb[2] == "-": 
       posRE = refGenome[prb[6]+":" +prb[1]].rfind(RE,0,prb[3])
       xprbStart=prb[4]
       xprbEnd=prb[3]
       xfragEnd=posRE 
    else: 
       posRE = refGenome[prb[6]+":"+prb[1]].find(RE,prb[4])
       xprbStart=prb[3]
       xprbEnd=prb[4]
       xfragEnd=posRE + 6
    entry = [prb[0], xprbStart, xprbEnd, xfragEnd, prb[2]]
    for i in range( int(prb[3]/1000) - 3, int(prb[3]/1000) + 4):   
        binid=prb[1] + ":" + str(i)    
        if not binid in probes: 
           probes[binid] = []
        probes[binid].append(entry)






############################################################################################################################
##Analysis Functions 
############################################################################################################################



def SamReader(fl):
    i = 0
    for rd in gzip.open(fl,"rt"): 
      if rd[0:3] != "@SQ": 
        #M00559:63:000000000-A5KK7:1:1101:10703:1212     16      HG19:chr16      49386963        179     106M45S *       0       0       SEQ QUAL AS:i:102        XS:i:0  XF:i:3  XE:i:2  NM:i:2
        cols = rd[:-1].split("\t")
        [segid,flag,chromo,pos,mapq,cigar] = cols[0:6]
        seq=cols[9]
        mapqbin="UNKNOWN"
        mapq = int(mapq)
        if flag == "4": 
           mapqbin = "UNMAPPED"
        elif mapq == 0: 
           mapqbin="AMB"
        elif mapq < 10:
           mapqbin = "LT10"
        elif mapq < 30:     
           mapqbin = "LT30"
        else:  
           mapqbin = "GT30"                
        if mapq > 0: 
           matched = parseCigar(cigar)['M']
           [genome,chromo] = chromo.split(":")
           if chromo+":"+str(int(pos)/1000) in lstBlackListed: 
              mapqbin = "BlackHole"
        else: 
           matched = 0
           [genome,chromo] = ["NA","NA"]
        i += 1      
        yield  [segid,flag,genome,chromo,int(pos), int(pos) + matched,mapq,mapqbin, seq, rd] 
  


def SamReader2(fl):
    i = 0
    #print(fl)
    samfile = pysam.AlignmentFile(fl, "rb")
    for read in samfile.fetch():
      rd = read.to_string()  
      if rd[0:3] != "@SQ": 
        #M00559:63:000000000-A5KK7:1:1101:10703:1212     16      HG19:chr16      49386963        179     106M45S *       0       0       SEQ QUAL AS:i:102        XS:i:0  XF:i:3  XE:i:2  NM:i:2
        cols = rd[:-1].split("\t")
        [segid,flag,chromo,pos,mapq,cigar] = cols[0:6]
        seq=cols[9]
        mapqbin="UNKNOWN"
        mapq = int(mapq)
        if flag == "4": 
           mapqbin = "UNMAPPED"
        elif mapq == 0: 
           mapqbin="AMB"
        elif mapq < 10:
           mapqbin = "LT10"
        elif mapq < 30:     
           mapqbin = "LT30"
        else:  
           mapqbin = "GT30"                
        if mapq > 0: 
           matched = parseCigar(cigar)['M']
           [genome,chromo] = chromo.split(":")
           if chromo+":"+str(int(pos)/1000) in lstBlackListed: 
              mapqbin = "BlackHole"
        else: 
           matched = 0
           [genome,chromo] = ["NA","NA"]
        i += 1      
        yield  [segid,flag,genome,chromo,int(pos), int(pos) + matched,mapq,mapqbin, seq, rd] 


def CheckProbe(chromo,xstart,xend, flag):
    pos1kbp = chromo+":"+str(int(xstart/1000))
    if pos1kbp not in probes:
       return ["UNKNOWN","UNKNOWN"]
    bestDst = 10000
    #Find the probe closest to the read start
    bestEntry = []
    for [prbid,prbStart,prbEnd,prbFragEnd,prbTp] in probes[pos1kbp]: 
        if prbTp == "+": 
           xdst = abs(xstart - prbStart) 
        else:
           xdst = abs(xend - prbStart) 
        if xdst < bestDst: 
           bestDst = xdst
           bestEntry = [prbid,prbStart,prbEnd,prbFragEnd,prbTp]
    [prbid,prbStart,prbEnd,prbFragEnd,prbTp] = bestEntry
    if prbTp == "+": 
       if flag != "0":
          return [prbid,"FLAG-Mismatch:"+flag]
       else:  
          if abs(xstart - prbStart) < 5 and abs(xend - prbFragEnd) < 5: 
             return [prbid,"VALID"]
          elif abs(xstart - prbStart) < 5 and abs(xend - prbEnd) < 5:  
             return [prbid,"MISSANEAL"] 
          elif abs(xstart - prbStart) < 5 and xend > prbFragEnd:  
             return [prbid,"UNDIG"]
          else: 
             return [prbid,"NO-MATCH"]
    else: 
       if flag != "16":
          return [prbid,"FLAG-Mismatch:"+flag]
       else: 
          if abs(xend - prbStart) < 5 and abs(xstart - prbFragEnd) < 5: 
             return [prbid,"VALID"]
          elif abs(xend - prbStart) < 5 and abs(xstart - prbEnd) < 5:  
             return [prbid,"MISSANEAL"] 
          elif abs(xend - prbStart) < 5 and xstart < prbFragEnd:  
             return [prbid,"UNDIG"]
          else: 
             return [prbid,"NO-MATCH"]    


def getREFrag(chromo,  pos, resEnz): 
   B = refGenome[chromo].find(resEnz, pos -  len(resEnz) + 1)
   if B > 0: 
      B = B + len(resEnz) 
   else: 
      B = len(refGenome[chromo]) 
   A = -1
   offset = 0
   while A < 0 and pos > offset: 
      A = refGenome[chromo].rfind(resEnz, pos - offset - 1010, pos - offset )
      offset += 1000
   A=min(max(1,A), len(refGenome[chromo]))
   B=min(max(1,B), len(refGenome[chromo]))
   return [chromo + ":" + str(A) + "+" + str(B-A), A,B]



############################################################################################################################
## LOAD Sam files
############################################################################################################################

valid={}
statsLib={}
statsFrag={}
statsPrbs={}
umis={}

#### Count UMI's
iprog=0
if useUMI: 
   #for u in SeqIO.parse(gzip.open("../UPLOADS/NW_ST1/"+samp+"_I2.fastq.gz", "rt"),"fastq"): 	
   for u in SeqIO.parse(gzip.open(inputDir+"/"+samp+"_I2.fastq.gz", "rt"),"fastq"): 
       umis[u.id] = str(u.seq)
       if iprog % 100000 == 0: 
          print (["Loading UMI", samp, iprog])
       iprog += 1

#### Count input number of reads from FASTQ
rd2ct=0
iprog=0
#for x in gzip.open("../UPLOADS/NW_ST1/"+samp+"_R2.fastq.gz", "rt"): 
for x in gzip.open(inputDir+"/"+samp+"_R2.fastq.gz", "rt"): 
    #if iprog > 1000000: 
    #   break; 
    rd2ct += 0.25
    if iprog % 1000000 == 0: 
       print (["Counting R2", samp, iprog])
    iprog += 0.25


statsPrbs["TotReads\tFastq"]=rd2ct
print("the printng value of TotalReads tab delimated Fastq")
print(rd2ct)


#### Process R2 (Determine if probe was present and if it was digested and completed properly ) 
iprog=0
#fl=outputDirSam+"/"+samp+"_R2.sam.gz"
fl=outputDirSam+"/"+samp+"_R2.sam.gz"
for [rdid,flag,genome, chromo,xstart, xend,mapq,mapqbin, seq, sam] in SamReader(fl): 
    #if iprog > 1000000: 
    #   break; 
    if iprog % 1000000 == 0: 
       print( "\t".join(map(str, ["Loading R2", samp, iprog, str(int(100*iprog/rd2ct))+"%", ]))) 
    iprog += 1    
    [a,b] =  CheckProbe(chromo,xstart,xend, flag)
    if b == "VALID":
       valid[rdid] = a
    if a != "UNKNOWN": 
       if a+"\t"+b in statsPrbs: 
          statsPrbs[a+"\t"+b] += 1
       else: 
          statsPrbs[a+"\t"+b] = 1    

#### Process R1 (Look for the complementary sequence) 
#print(statsPrbs)
#fl="../UPLOADS/NW_ST2/"+samp+"_R1.sam.gz"
fl=outputDirSam+"/"+samp+"_R1.sam.gz"
statsFrags={}
iprog=1
ValidCov={}
for [rdid,flag,genome, chromo,xstart, xend,mapq, mapqbin, seq, sam] in SamReader(fl): 
    #if iprog > 1000000: 
    #   break; 
    if iprog % 1000000 == 0: 
       print( "\t".join(map(str, ["Loading R1", samp, iprog, str(int(100*iprog/rd2ct))+"%", ]))) 
    iprog += 1
    isValid = rdid in valid
    statid="\t".join([samp,"R1",genome, chromo, str(1000000*((xstart + xend)/2000000)), mapqbin, str(isValid)])
    if statid in statsLib: 
       statsLib[statid] += 1
    else: 
       statsLib[statid] = 1
    if isValid and mapqbin not in ["BlackHole", "UNKNOWN","UNMAPPED", "AMB", "LT10"]:
       frg = getREFrag(genome+":"+chromo,int((xstart+xend)/2),RE)
       prb = valid[rdid]       
       vstat = "\t".join( map(str, [prb,flag,genome, chromo,xstart, xend,frg[0]] )  )
       if not vstat in ValidCov: 
          ValidCov[vstat] = 1
       else: 
          ValidCov[vstat] += 1 
       entryid = frg[0] + "\t" + prb
       if not entryid in statsFrags: 
          statsFrags[entryid] = {"LT10":0, "LT30":0, "GT30":0, "RE_Start":frg[1], "RE_End":frg[2],"cov":{},"umi":{}}      
       if useUMI and mapq > 10: 
          u = umis[rdid]
          if u in statsFrags[entryid]["umi"]:
             statsFrags[entryid]["umi"][u] += 1 
          else: 
             statsFrags[entryid]["umi"][u] = 1 
       statsFrags[entryid][mapqbin] += 1
       iPos="XX"
       if flag == "0":           
          iPos = "F:"+ str(frg[2] - xstart)
       if flag == "16": 
          iPos = "R:"+ str(xend - frg[1])
       if iPos in statsFrags[entryid]["cov"]:
          statsFrags[entryid]["cov"][iPos] += 1
       else: 
          statsFrags[entryid]["cov"][iPos] = 1


############################################################################################################################
## Write Results
############################################################################################################################
#print(statsFrags)
#print ("Writing Stats: " +samp)


out = open(outputDirStats+"/"+samp+"_ProbeStats.tab","w")
for x in statsPrbs: 
    out.write(samp + "\t" + x + "\t" + str(statsPrbs[x]) + "\n")

out.close()

#out = open("../UPLOADS/STATS/"+samp+"_ReadStats.tab","w")
out = open(outputDirStats+"/"+samp+"_ReadStats.tab","w")
for x in statsLib: 
    out.write(x + "\t" + str(statsLib[x]) + "\n")

out.close()

#out = open("../UPLOADS/STATS/"+samp+"_FragStats.tab","w")
out = open(outputDirStats+"/"+samp+"_FragStats.tab","w")
outUMI = open(outputDirStats+"/"+samp+"_UmiStats.tab","w")
out.write( "\t".join(["SAMP","Frag",  "GENOME","CHROMO", "FragStart","FragEnd","PROBE","ProbeChromo","ProbePosition","ProbeGenome" ,"LT10","LT30","GT30", "ctpos","ctPos5","ctPos10","ctumi","ctUmi5","ctUmi10","maxPrbMatches", "maxPrbOverPos", "maxPrbTp","Positions","Umis"]) + "\n")

for x in statsFrags:
    [frag,prb ]=x.split("\t")
    [genome, chromo, RE_frag] = frag.split(":")
    frgStats = statsFrags[x]
    LT10=frgStats["LT10"]
    LT30=frgStats["LT30"]
    GT30=frgStats["GT30"]
    FragStart=frgStats["RE_Start"]
    FragEnd=frgStats["RE_End"]    
    ctpos=len(frgStats["cov"])
    if useUMI: 
       ctumi=len(frgStats["umi"])
    else: 
       ctumi=0
    maxPrbMatches=0 
    maxPrbOverPos=0 
    maxPrbTp="UNK"
    seq = refGenome[genome+":"+chromo][FragStart:FragEnd] 
    maxPrbMatches=0
    maxPrbOverPos=-1
    maxPrbTp="NONE"
    probe = prbsALL[prb]
    prb_seq = refGenome[probe[6] + ":"+probe[1] ][probe[3]:probe[4]]
    for i in range(1, len(seq) - len(prb_seq) - 1):
           subseq= seq[i:(len(prb_seq)+i)  ]
           dst = 1.0*distance(prb_seq, subseq) 
           score = 1 - dst/len(prb_seq)
           if score > maxPrbMatches:
              maxPrbMatches = score
              maxPrbOverPos = i + FragStart
              maxPrbTp="FOR"
    if len(seq) > 500: 
       seq="Too Long"
       covct="Too Long"
    ctPos5 = 0
    ctPos10 = 0
    xpos = ""
    for x in frgStats["cov"]:  
      if frgStats["cov"][x] > 5: 
         xpos += x+":"+str(frgStats["cov"][x])+","
         ctPos5 += 1
      if frgStats["cov"][x] > 10: 
         ctPos10 += 1
    xpos = xpos[0:2000]
    xumi = ""
    ctUmi5 = 0
    ctUmi10 = 0
    for u in frgStats["umi"]:  
      xx=[samp,frag,prb,u,str(frgStats["umi"][u]), str(ctpos) ]
      outUMI.write("\t".join(xx)+"\n") 
      if frgStats["umi"][u] > 5: 
         xumi += u+":"+str(frgStats["umi"][u])+","
         ctUmi5 += 1
      if frgStats["umi"][u] > 10: 
         ctUmi10 += 1
    out.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,prb,probe[1],probe[3],probe[6],LT10,LT30,GT30, ctpos,ctPos5,ctPos10, ctumi,ctUmi5,ctUmi10, maxPrbMatches, maxPrbOverPos, maxPrbTp, xpos, xumi])) + "\n")


out.close()
outUMI.close()

out = open(outputDirStats+"/"+samp+"_ValidCovStats.tab","w")
for x in ValidCov: 
  ln= x+"\t"+str(ValidCov[x])+"\n"
  out.write(ln)


out.close()


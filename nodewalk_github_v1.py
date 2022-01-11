#!/usr/bin/env python
# coding: utf-8

# In[14]:


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
import os


# In[13]:


baitseq='GTCTATGTACTTGTGAATTATTTCACGT'

samp=sys.argv[1]
RE=sys.argv[2]
prbsALL = json.load(open(sys.argv[3]))
RefGenome = sys.argv[4]

inputDir = sys.argv[5]
outputDirSam =sys.argv[6]
outputDirStats = sys.argv[7]

#inputDir = "../UPLOADS/NW_ST1/"
#outputDirSam ="../UPLOADS/NW_ST2/"
#outputDirStats = "../UPLOADS/STATS/"




#print(prbsALL)

print('argument',samp,'argument2',RE,'argument3','not_include','arg4',RefGenome,'inputdir',inputDir,'samdir',outputDirSam,'outstats',outputDirStats)
#exit()

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


# In[3]:


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

   


#------------------------------------Brought changings here----------------------------------
probes={}
for xprb in prbsALL: 
    prb = prbsALL[xprb]
    szPrb = prb[4] - prb[3]
    if prb[2] == "-": 
       posRE = refGenome[prb[6]+":" +prb[1]].rfind(RE,0,prb[3])
      # print("number of hind3 sites"+str(refGenome[prb[6]+":" +prb[1]].count(RE,prb[3],posRE)))
      # print(posRE)
       xprbStart=prb[4]
       xprbEnd=prb[3]
       xfragEnd=posRE 
    else: 
       posRE = refGenome[prb[6]+":"+prb[1]].find(RE,prb[4])
      # print("number of hind3 sites-positive"+str(refGenome[prb[6]+":"+prb[1]].count(RE,prb[4],posRE)))
       xprbStart=prb[3]
       xprbEnd=prb[4]
       xfragEnd=posRE + 6
    entry = [prb[0], xprbStart, xprbEnd, xfragEnd, prb[2]]
    entry2 = prb[1]+"\t"+str(xprbStart)+"\t"+str(xfragEnd)+"\t"+str(prb[0])+"\t"+str(prb[2])+"\n"
    
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
           #matched = len(seq)
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


def CheckProbe(chromo,xstart,xend, flag,sam=""):
    pos1kbp = chromo+":"+str(int(xstart/1000))
    if pos1kbp not in probes:
       return ["UNKNOWN","UNKNOWN",0,0,0]
    bestDst = 10000
    #Find the probe closest to the read start
    bestEntry = []
    
    #[['RPB002_ARNTL', 13317801, 13317823, 13317849, '+']] probe structure is here
    #print("the pos1kbp value is ="+pos1kbp)
    
    for [prbid,prbStart,prbEnd,prbFragEnd,prbTp] in probes[pos1kbp]: 
        if prbTp == "+": 
           xdst = abs(xstart - prbStart) 
        else:
           xdst = abs(xend - prbStart) 
        if xdst < bestDst: 
           bestDst = xdst
           bestEntry = [prbid,prbStart,prbEnd,prbFragEnd,prbTp]
          # print('the best entry is',bestEntry)
    [prbid,prbStart,prbEnd,prbFragEnd,prbTp] = bestEntry
    if prbTp == "+": 
       if flag != "0":
          return [prbid,"FLAG-Mismatch:"+flag,prbStart,prbFragEnd,prbTp]
       else:  
          if abs(xstart - prbStart) < 5 and abs(xend - prbFragEnd) < 5: 
             return [prbid,"VALID",prbStart,prbFragEnd,prbTp]
          elif abs(xstart - prbStart) < 5 and abs(xend - prbEnd) < 5:
             #print("opening the file")   
            # with open('Read2p.sam', 'a+') as the_file:
             #   the_file.write(sam)   
             return [prbid,"MISSANEAL",prbStart,prbFragEnd,prbTp]
          elif abs(xstart - prbStart) < 5 and xend > prbFragEnd:
            # with open('Read3p.sam', 'a+') as the_file2:
             #   the_file2.write(sam)    
             return [prbid,"UNDIG",prbStart,prbFragEnd,prbTp]
          else:
            # print("no match found")
             return [prbid,"NO-MATCH",prbStart,prbFragEnd,prbTp]
        
    else: 
       if flag != "16":
          return [prbid,"FLAG-Mismatch:"+flag,prbStart,prbFragEnd,prbTp]
       else: 
          if abs(xend - prbStart) < 5 and abs(xstart - prbFragEnd) < 5: 
             return [prbid,"VALID",prbStart,prbFragEnd,prbTp]
          elif abs(xend - prbStart) < 5 and abs(xstart - prbEnd) < 5:
            # print("opening the file")    
            # with open('Read2m.sam', 'a+') as the_file:
             #   the_file.write(sam)   
             return [prbid,"MISSANEAL",prbStart,prbFragEnd,prbTp]
          elif abs(xend - prbStart) < 5 and xstart < prbFragEnd:
             #with open('Read3m.sam', 'a+') as the_file2:
              #  the_file2.write(sam) 
             return [prbid,"UNDIG",prbStart,prbFragEnd,prbTp]
          else: 
             #print("no match found")
             return [prbid,"NO-MATCH",prbStart,prbFragEnd,prbTp]  



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



# In[4]:


valid={}
statsLib={}
statsFrag={}
statsPrbs={}
umis={}
validr1={}
statsLibr1={}
statsFragr1={}
statsPrbsr1={}
umisr1={}
#### Count UMI's
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


# In[5]:


#baitseq='GTCTATGTACTTGTGAATTATTTCACGT'
valid_R2 ={}
valid_R2_seq={}
bait_Valid_R2= {}
valid ={}
valid_R2_getstate ={}
#### Process R2 (Determine if probe was present and if it was digested and completed properly ) 
iprog=0
#fl=outputDirSam+"/"+samp+"_R2.sam.gz"
fl=outputDirSam+"/"+samp+"_R2.sam.gz"
for [rdid,flag,genome, chromo,xstart, xend,mapq,mapqbin, seq, sam] in SamReader(fl): 
    #if iprog > 1000000: 
    #break; 
    if iprog % 1000000 == 0: 
       print( "\t".join(map(str, ["Loading R2", samp, iprog, str(int(100*iprog/rd2ct))+"%", ]))) 
    iprog += 1 
    #[a,b,frgstart,fragend,tp]
    [a,b,frgstart,fragend,tp] =  CheckProbe(chromo,xstart,xend, flag,seq)
    if b == "VALID":
       if(sam.count(baitseq)>0):
        if(sam.count(RE)>0):
           bait_Valid_R2[rdid]= seq.split(RE)[0]+'\t'+str(seq.find(RE))
           #with open('output/'+samp+'.sam', 'a+') as the_file:
             #   the_file.write(rdid+'\t'+seq.split(RE)[0]+'\t'+str(seq.find(RE))+'\n')
        
        
        #sam= sam.replace("GTCTATGTACTTGTGAATTATTTCACGT",'')
        #print("there is sequence...............");
      # valid[rdid] = a
       valid_R2_getstate[rdid] = a
       #if (seq.count(RE)>0):     
       # valid_R2[rdid]= chromo+'\t'+str(xstart)+'\t'+str(xend)+'\t'+str(mapq)+'\t'+str(mapqbin)+"\t"+str(chromo)+"\t"+genome+"\t"+flag+"\t"+"1"
       # valid_R2_seq[rdid] = sam
     #  else:
      #  valid_R2[rdid]= chromo+'\t'+str(xstart)+'\t'+str(xend)+'\t'+str(mapq)+'\t'+str(mapqbin)+"\t"+str(chromo)+"\t"+genome+"\t"+flag+"\t"+"0"
      #  valid_R2_seq[rdid] = sam
       #print(chromo)
       #if(mapqbin!="GT30" and mapqbin!="LT30" ):
          #  print(mapqbin)
    if a != "UNKNOWN": 
       if a+"\t"+b in statsPrbs: 
          statsPrbs[a+"\t"+b] += 1
       else: 
          statsPrbs[a+"\t"+b] = 1    


# In[6]:



print("Extracting valid reads and remapping to "+RefGenome)
t=0
ft=0
ncount =0
bait_flag=False;
bait_tempID=''


file_fq = gzip.open(inputDir+'/'+samp+'_Chimeric.fastq.gz',"wt") 
for x in gzip.open(inputDir+"/"+samp+"_R2.fastq.gz", "rt"):
 
   # if(ft==4):break
        
    if(t==0):
        tid = x
        x = x.split('@')[1].split(' ')[0] 
        #print("start of new read"+"the value of t is "+str(t)+"\t"+x)
        if(x in bait_Valid_R2):
            bait_flag=True
            bait_tempID=x
            #print('valid read is there')
            ncount= ncount+1         
            sposition = int(bait_Valid_R2[x].split('\t')[1])
            #print(x)
            #print(bait_Valid_R2[x].split('\t')[0])
            #x="@"+x
            file_fq.write(tid)
            
        t=t+1
    elif(t==3):
        if(bait_flag==True):
            #pqbit = x[len(x)-sposition-1:len(x)]+'\t'+str(len(x[len(x)-sposition-1:len(x)]))
            ###pqbit = x[len(x)-sposition-1:len(x)]
            ###print(pqbit)
            x = x[len(x)-sposition-1:len(x)]
            bait_flag=False
            file_fq.write(x)
        t=0
    elif(t==1):
        if(bait_flag==True):
            #pseq = x[len(x)-sposition-1:len(x)]+'\t'+str(len(x[len(x)-sposition-1:len(x)]))
           ### pseq = x[len(x)-sposition-1:len(x)]
           ### print(pseq)
            x= x[len(x)-sposition-1:len(x)]
            file_fq.write(x)
            #print(pseq+'\t'+str(len(pseq)))
           # print(bait_Valid_R2[bait_tempID]+'\t'+str(len(bait_Valid_R2[bait_tempID].split('\t')[0]))+'end here \n\n\n')
           
            
        t=t+1
    else:
        if(bait_flag==True):
            file_fq.write(x)
           ### print(x)
        
            xxx=0
        #print("the value of t is \t"+str(t)+"\t"+x)
        #if(t==2):
        #print(x)
        t=t+1
        
    ft=ft+1;
file_fq.close()


# In[7]:


print("Remapping the chimeric reads to "+RefGenome)
cfile = inputDir+'/'+samp+'_Chimeric.fastq.gz '
cmd = "./bwa bwasw -t 32 "+RefGenome + " "+ cfile + " | gzip - >  "+inputDir+'/'+samp+'_Chimeric.sam.gz'
t = os.system(cmd)
os.remove(inputDir+'/'+samp+'_Chimeric.fastq.gz')


# In[8]:


fl= "R2_map_again.sam.gz"
fl=inputDir+'/'+samp+'_Chimeric.sam.gz'

valid={}
valid_R2 ={}
valid_R2_seq={}
bait_Valid_R2= {}
#### Process R2 (Determine if probe was present and if it was digested and completed properly ) 
iprog=0
#fl=outputDirSam+"/"+samp+"_R2.sam.gz"
#fl=outputDirSam+"/"+samp+"_R2.sam.gz"
for [rdid,flag,genome, chromo,xstart, xend,mapq,mapqbin, seq, sam] in SamReader(fl): 
    #if iprog > 1000000: 
    #break; 
    if iprog % 1000000 == 0: 
       print( "\t".join(map(str, ["Loading R2", samp, iprog, str(int(100*iprog/rd2ct))+"%", ]))) 
    iprog += 1 
    
    #[a,b,frgstart,fragend,tp]
    [a,b,frgstart,fragend,tp] =  CheckProbe(chromo,xstart,xend, flag,seq)
    
    if chromo != "NA":
       valid[rdid] = valid_R2_getstate[rdid]
       a = valid_R2_getstate[rdid]
       valid_R2[rdid]= chromo+'\t'+str(xstart)+'\t'+str(xend)+'\t'+str(mapq)+'\t'+str(mapqbin)+"\t"+str(chromo)+"\t"+genome+"\t"+flag+"\t"+"0"
       valid_R2_seq[rdid] = sam
       #print(chromo)
       #if(mapqbin!="GT30" and mapqbin!="LT30" ):
          #  print(mapqbin)
    if a != "UNKNOWN": 
       if a+"\t"+b in statsPrbs: 
          statsPrbs[a+"\t"+b] += 1
       else: 
          statsPrbs[a+"\t"+b] = 1 


# In[9]:


#### Process R1 (Look for the complementary sequence) 
#print(statsPrbs)
#fl="../UPLOADS/NW_ST2/"+samp+"_R1.sam.gz"


R2_duplicated = {}
rfilename = "_R1.sam.gz";
filewritename = "R1"

fl=outputDirSam+"/"+samp+"_R1.sam.gz"
statsFrags={}
iprog=1
ValidCov={}
flag_c=0
flag_c2=0
for [rdid,flag,genome, chromo,xstart, xend,mapq, mapqbin, seq, sam] in SamReader(fl):
    
    flag_R1_n =0
    flag_R1_ex =0;
    flag_R2_n=0;
    flag_R2_ext = 0;
    
    
    cflag_n=0
    cflag_ex =0
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
        
        
    
    #---------------Upgrading according to the last meeting final ------------------------------------------------------------------------------------------
       if rdid not in R2_duplicated:
            R2_duplicated[rdid]=1
       else:
            R2_duplicated[rdid]= R2_duplicated[rdid]+1
       frg = getREFrag(genome+":"+chromo,int((xstart+xend)/2),RE)
    
    
       frg2 = getREFrag(genome+":"+chromo,int(frg[1]-4),RE)
       frg3 = getREFrag(genome+":"+chromo,int(frg[2]+4),RE)    
    
       frg_ext = [frg[0].split(':')[0]+":"+frg[0].split(':')[1]+":"+str(frg2[1])+"+"+str(frg3[2]-frg2[1]), frg2[1],frg3[2]]
        
    
       R2_start = int(valid_R2[rdid].split('\t')[1])
       R2_end = int(valid_R2[rdid].split('\t')[2])
       R2_chromo = valid_R2[rdid].split('\t')[5]
       R2_genome =  valid_R2[rdid].split('\t')[6]
       R2_mapbin =  valid_R2[rdid].split('\t')[4]
       R2_flag =   valid_R2[rdid].split('\t')[7]
       R2_inside_hint =   valid_R2[rdid].split('\t')[8]
       R2_seq =  valid_R2_seq[rdid]
      # R2_seq = R2_seq[9:30]
    
       lt30=0
       gt30=0
        
       if(R2_mapbin=='LT30'):
        lt30=1
       if(R2_mapbin=='GT30'):
        gt30=1
       #if R2_mapbin != 'LT30' and R2_mapbin != 'GT30':
       # print(R2_mapbin) 
    
       if(xstart >=frg[1] and xend <= frg[2]):
            flag_R1_n=1
       elif(xstart >=frg_ext[1] and xend <= frg_ext[2]):
            flag_R1_ex=1
       if(R2_start >=frg[1] and R2_end <= frg[2]):
            flag_R2_n =1
       elif(R2_start >=frg_ext[1] and R2_end <= frg_ext[2]):
            flag_R2_ext =1   
    
    
    # ---------------------------------------------------------End of Last meeting upgradiation-------------------------------------------------------------
    
    #Writing Majoe code here--------------------------------------------------------------------------

       if (xstart >=frg[1] and xend <= frg[2] and R2_start >=frg[1] and R2_end <= frg[2]):
            flag_c = flag_c+1
            cflag_n =1
            
            
    # now looking at the neigbhoring fragment ...........................................................
       frg = getREFrag(genome+":"+chromo,int((xstart+xend)/2),RE)  
        
       frg2 = getREFrag(genome+":"+chromo,int(frg[1]-4),RE)
    
       frg3 = getREFrag(genome+":"+chromo,int(frg[2]+4),RE)
        
       frg = [frg[0].split(':')[0]+":"+frg[0].split(':')[1]+":"+str(frg2[1])+"+"+str(frg3[2]-frg2[1]), frg2[1],frg3[2]]
        
       if (xstart >=frg[1] and xend <= frg[2] and R2_start >=frg[1] and R2_end <= frg[2]):
            flag_c2 = flag_c2+1
            cflag_ex =1
        
       frg = getREFrag(genome+":"+chromo,int((xstart+xend)/2),RE)
    
       #frg = [frg[0].split(':')[0]+":"+frg[0].split(':')[1]+":"+str(frg2[1])+"+"+str(frg3[2]-frg2[1]), frg2[1],frg3[2]]
        
       prb = valid[rdid]       
       vstat = "\t".join( map(str, [prb,flag,genome, chromo,xstart, xend,frg[0]] )  )
       if not vstat in ValidCov: 
          ValidCov[vstat] = 1
       else: 
          ValidCov[vstat] += 1 
       entryid = frg[0] + "\t" + prb
       if not entryid in statsFrags: 
          statsFrags[entryid] = {"LT10":0, "LT30":0, "GT30":0, "RE_Start":frg[1], "RE_End":frg[2], "eRE_Start":frg2[1], "eRE_End":frg3[2],"R1_nFrg":flag_R1_n,"R1_extFrg":flag_R1_ex,"R2_nFrg":flag_R2_n,"R2_extFrg":flag_R2_ext,"cov":{},"cov2":{},"umi":{},"nfrag":cflag_n,"efrag":cflag_ex,"R2_fragment":{}}
       else:
          statsFrags[entryid]["nfrag"] = statsFrags[entryid]["nfrag"]+cflag_n
          statsFrags[entryid]["efrag"] = statsFrags[entryid]["efrag"]+cflag_ex        
          statsFrags[entryid]["R1_nFrg"] = statsFrags[entryid]["R1_nFrg"]+flag_R1_n
          statsFrags[entryid]["R1_extFrg"] = statsFrags[entryid]["R1_extFrg"]+flag_R1_ex
          statsFrags[entryid]["R2_nFrg"] = statsFrags[entryid]["R2_nFrg"]+flag_R2_n
          statsFrags[entryid]["R2_extFrg"] = statsFrags[entryid]["R2_extFrg"]+flag_R2_ext
        
        
       #if (cflag_n==0):
       #   with open('output/'+samp+'.sam', 'a+') as the_file:
        #     the_file.write((R2_seq))
       #-----------------------------Updating the files according to anita told Us
       #if(R2_chromo!=chromo):
            #print("there are difference in chromosoms"+R2_chromo+"::"+chromo)
            #break      
            
       frg_r2 = getREFrag(R2_genome+":"+R2_chromo,int((R2_start+R2_end)/2),RE)
        
       entryid2 = frg_r2[0] + "\t" + prb
    #----Codefor cttoc values extractions from the data-----------------------------------------------------------------------------Important last change ---------------------
        
       iPos_r2="XX"
       if R2_flag == "0":           
          iPos_r2 = "F:"+ str(frg_r2[2] - R2_start)
       if R2_flag == "16": 
          iPos_r2 = "R:"+ str(R2_end - frg_r2[1])
       iPos_r2= entryid2+":"+iPos_r2
       #print(iPos_r2)
       if iPos_r2 in statsFrags[entryid]["cov2"]:
          statsFrags[entryid]["cov2"][iPos_r2] += 1
       else: 
          statsFrags[entryid]["cov2"][iPos_r2] = 1

    #----------------------------------------------------------------------------------End of the important cttoc calculation---------------------------------------------------
    
       if entryid2 not in statsFrags[entryid]["R2_fragment"]:
            entrval = str(1)+":"+str(lt30)+":"+str(gt30)+":"+str(R2_inside_hint)+':'+str(cflag_n)+":"+str(cflag_n)+":"+str(cflag_ex)
            statsFrags[entryid]["R2_fragment"][entryid2] = entrval
       else:
            # 1:LT30:GT30
            nval = int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[0])
            lt30x = int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[1])
            gt30x = int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[2])
            r2hind3 = int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[3])
            chimericflag = int(int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[4]))
            chimericcount = int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[5])
            chimericcount_ext = int(statsFrags[entryid]["R2_fragment"][entryid2].split(':')[6])
           # print(str(chimericcount+cflag_n))
            
           # if cflag_n != chimericflag :
           #    print("Chimeric conflict Error...... in code .......")            
            entrval = str(nval+1)+":"+str(lt30x+lt30)+":"+str(gt30x+gt30)+":"+str(r2hind3+int(R2_inside_hint))+":"+str(cflag_n)+":"+str(chimericcount+cflag_n)+":"+str(chimericcount_ext+cflag_ex)
            
            statsFrags[entryid]["R2_fragment"][entryid2] = entrval
                        #with open('output/R2seq.sam', 'a+') as the_file:
             #   the_file.write(validprob[rdid]+"\n")
            #statsFrags[entryid]["R2_fragmentbin"][R2_mapbin]+=1
            
            #statsFrags[entryid]["R2_fragmentbin"][R2_mapbin] += 1
       
            
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

print(flag_c)


# In[11]:


############################################################################################################################
## Write Results
############################################################################################################################
#print(statsFrags)
#print ("Writing Stats: " +samp)


out = open(outputDirStats+"/"+samp+filewritename+"_ProbeStats.tab","w")
for x in statsPrbs: 
    out.write(samp + "\t" + x + "\t" + str(statsPrbs[x]) + "\n")

out.close()

out_chimeric = open(outputDirStats+"/"+samp+filewritename+"_Chimaric_fragstats.tab","w")


out_chimeric.write( "\t".join(["SAMP","R1_Frag",  "R1_GENOME","R1_CHROMO", "R1_FragStart","R1_FragEnd","R2_Frag","R2_FragStart","R2_FragEnd","R2_GT30",'R2_LT30','R2_ctTot',"R2_ctPos","ctPos_c","R2_hind3_reads","Is chimeric","\n"]))



#out = open("../UPLOADS/STATS/"+samp+"_ReadStats.tab","w")
out = open(outputDirStats+"/"+samp+filewritename+"_ReadStats.tab","w")
for x in statsLib: 
    out.write(x + "\t" + str(statsLib[x]) + "\n")

out.close()

#out = open("../UPLOADS/STATS/"+samp+"_FragStats.tab","w")
out = open(outputDirStats+"/"+samp+filewritename+"_FragStats.tab","w")
outUMI = open(outputDirStats+"/"+samp+filewritename+"_UmiStats.tab","w")
#out.write( "\t".join(["SAMP","Frag",  "GENOME","CHROMO", "FragStart","FragEnd","PROBE","ProbeChromo","ProbePosition","ProbeGenome" ,"LT10","LT30","GT30", "ctpos","ctPos5","ctPos10","ctumi","ctUmi5","ctUmi10","maxPrbMatches", "maxPrbOverPos", "maxPrbTp","Positions","Umis"]) + "\n")

# comments of the before meeting with anita 29 dec, 2021
#out.write( "\t".join(["SAMP","Frag",  "GENOME","CHROMO", "FragStart","FragEnd","PROBE","ProbeChromo","ProbePosition","ProbeGenome","ctTot", "ctpos","R1 Nfrag","R1 ExtFrag","R2 Nfrag","R2_ExtFrag","Is Chimeric","\n"]))
out.write( "\t".join(["SAMP","Frag",  "GENOME","CHROMO", "R1_FragStart","R1_FragEnd","PROBE","ProbeChromo","ProbePosition","ProbeGenome","R1 ctTot", "R1 ctPos","Is Chimeric","\n"]))


     #     statsFrags[entryid]["nfrag"] = statsFrags[entryid]["R1_nFrg"]+flag_R1_n
       #   statsFrags[entryid]["efrag"] = statsFrags[entryid]["R1_extFrg"]+flag_R1_ex
        #  statsFrags[entryid]["nfrag"] = statsFrags[entryid]["R2_nFrg"]+flag_R2_n
         # statsFrags[entryid]["efrag"] = statsFrags[entryid]["R2_extFrg"]+flag_R2_ext
        
        
for x in statsFrags:
    [frag,prb ]=x.split("\t")
    [genome, chromo, RE_frag] = frag.split(":")
    frgStats = statsFrags[x]
    LT10=frgStats["LT10"]
    LT30=frgStats["LT30"]
    GT30=frgStats["GT30"]
    FragStart=frgStats["RE_Start"]
    FragEnd=frgStats["RE_End"]    
    
    FragStarte=frgStats["eRE_Start"]
    FragEnde=frgStats["eRE_End"] 
    nfrag_val = frgStats["nfrag"]
    efrag_val = frgStats["efrag"]
    R2_ctPos = len( frgStats["cov2"])
    R2_frg_start = 0
    R2_frg_end = 0;
    
    R1_nfrg_val = frgStats["R1_nFrg"]
    R1_extfrg_val = frgStats["R1_extFrg"]    
    
    R2_nfrg_val = frgStats["R2_nFrg"]
    R2_extfrg_val = frgStats["R2_extFrg"]  

    if (int(nfrag_val)>0 or int(efrag_val)>0):
        nfrag=False
    else:
        nfrag=True
 

    ccflag=0
    
    R1_ctpos_cvalue = len(frgStats["cov"].values())
    for item in frgStats["R2_fragment"]:
        
        #print(item+":"+str(frgStats["R2_fragment"][item])+"length of key is "+str(len(frgStats["R2_fragment"])))
        # print(item.split(':')[2].split('\t')[0].split('+'))
        R2_frg_start = item.split(':')[2].split('\t')[0].split('+')[0]
        R2_frg_end = int(item.split(':')[2].split('\t')[0].split('+')[0])  + int(item.split(':')[2].split('\t')[0].split('+')[1])
        R2_frag_info  = item.split('\t')[0];
        R2_GT30 = frgStats["R2_fragment"][item].split(':')[2] 
        R2_LT30 = frgStats["R2_fragment"][item].split(':')[1]
        R2_hint3_reads = frgStats["R2_fragment"][item].split(':')[3]
        Ischimeric_val = int(frgStats["R2_fragment"][item].split(':')[4])
        
        Ischimeric_count = int(frgStats["R2_fragment"][item].split(':')[5])
        Ischimeric_count_ext = int(frgStats["R2_fragment"][item].split(':')[6])
        #print(Ischimeric_count)
        
        cov2value = len(list(v for k,v in frgStats["cov2"].items() if item in k))
        
        
        if Ischimeric_count <1 and Ischimeric_count_ext <1:
            Ischimeric_flag=True
            both_ctpos = (R1_ctpos_cvalue+cov2value)
        else:
            Ischimeric_flag=False
            both_ctpos = (cov2value)
        
        
        if(ccflag==0):
            if(int(R2_GT30)>0 or int(R2_LT30)>0):
                out_chimeric.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,R2_frag_info,R2_frg_start,R2_frg_end,R2_GT30,R2_LT30,(int(R2_LT30)+int(R2_GT30)), cov2value,both_ctpos,R2_hint3_reads,Ischimeric_flag])) + "\n")
            ccflag=1
        else:
           # cov2value = len(list(v for k,v in frgStats["cov2"].items() if item in k))
            if(int(R2_GT30) >0 or int(R2_LT30)>0):
                out_chimeric.write( "\t".join(map(str,['...', '...',  '...', '...','...','...',R2_frag_info,R2_frg_start,R2_frg_end,R2_GT30,R2_LT30,(int(R2_LT30)+int(R2_GT30)), cov2value,both_ctpos,R2_hint3_reads,Ischimeric_flag])) + "\n")
            

    
    

    
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
            
     #out_chimeric.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,R2_frag_info,R2_frg_start,R2_frg_end,(int(R2_LT30)+int(R2_GT30)), R2_ctPos,(R2_ctPos+ctpos),R2_hint3_reads,nfrag])) + "\n")

   # out_chimeric.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,FragStarte,FragEnde,(LT30+GT30), ctpos,])) + "\n")
    
            
    #out.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,prb,probe[1],probe[3],probe[6],LT10,LT30,GT30, ctpos,ctPos5,ctPos10, ctumi,ctUmi5,ctUmi10, maxPrbMatches, maxPrbOverPos, maxPrbTp, xpos, xumi,nfrag,str(int(efrag)+int(nfrag))])) + "\n")
   
    #out.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,prb,probe[1],probe[3],probe[6],(LT30+GT30), ctpos,nfrag,str(int(efrag)+int(nfrag))])) + "\n")
   # out.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,prb,probe[1],probe[3],probe[6],(LT30+GT30), ctpos,R1_nfrg_val,R1_extfrg_val,R2_nfrg_val,R2_extfrg_val,nfrag])) + "\n")
    # lastly updated after anita meeting 27 dec
    out.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,prb,probe[1],probe[3],probe[6],(LT30+GT30), ctpos,nfrag])) + "\n")

   # out.write( "\t".join(map(str,[samp, frag,  genome, chromo, FragStart,FragEnd,prb,probe[1],probe[3],probe[6],LT10,LT30,GT30, ctpos,ctPos5,ctPos10, ctumi,ctUmi5,ctUmi10, maxPrbMatches, maxPrbOverPos, maxPrbTp, xpos, xumi])) + "\n")


out.close()
outUMI.close()

out = open(outputDirStats+"/"+samp+"_ValidCovStats.tab","w")
for x in ValidCov: 
  ln= x+"\t"+str(ValidCov[x])+"\n"
  out.write(ln)

out.close()
out_chimeric.close()
os.remove(inputDir+'/'+samp+'_Chimeric.sam.gz')
print("Done")


#!/usr/bin/env python3
#syntax: python run_XMMspectra.py
#written by Ben,modified by Tong
# written on 2018-6-13 to run XMM products on SgrA PWN for 5 XMM obsID
# most recently updated in May 2019 for Sgr B2 XMM data sets
# point source or extended source: change flag in afrgen: extendedsource=yes/no; currently set to no for point source

import sys
import os
import string
from pathlib import Path

# ------obsID List---------------
path="/Users/baotong/xmm"
obsID1 = "0111310101"
obsID2 = "0111970701"
obsID3 = "0694641401"
# -------------------------------

# ------choose obsID-------------
obsList = [obsID2]
# -------------------------------

# ------detName List-------------
det1 = "mos1"
det2 = "mos2"
det3 = "pn"
# -------------------------------

# ------choose det --------------
detList = [det1,det2,det3]
# -------------------------------
# export SAS_ODF=/Users/baotong/xmm/0111310101/ODF
# cifbuild
# odfingest
# setenv SAS_CCF /Users/baotong/xmm/0111310101/cal/ccf.cif
# setenv SAS_ODF /Users/baotong/xmm/0111310101/cal/0494_0111310101_SCX00000SUM.SAS
# emchain
# epchain
# barycen withtable=yes table=pn_bary.fits: withsrccoordinates=yes srcra=32.46125 srcdec=-63.311111
process=1
spectra=1
ra=174.6117668;dec=3.3686417
for obsID in obsList:
   os.chdir(path+"/"+obsID)
   mypath=Path("./cal")

   if process:
      if mypath.is_dir():
         print("continued")
      else:
         os.system("mkdir cal")
         #------------ set environment -----------------
         os.chdir("./cal")
         #cmd = "alias SAS_ODF=" + path +"/"+ obsID + "/ODF"
         os.environ['SAS_ODF'] = path +"/"+ obsID + "/ODF"
         os.system("cifbuild")
         os.environ['SAS_CCF'] = path+"/"+obsID+"/cal/ccf.cif"
         os.system("odfingest")
         os.system("rm SAS.txt")
         os.system("ls *.SAS >SAS.txt")
         with open("SAS.txt",'r') as f:
            sasname=f.readline()
         print(sasname)
         os.environ['SAS_ODF'] = path + "/" + obsID + "/cal/"+sasname[0:-1]
         # # ---------------------------------------------
         # -----------processing for evt2 file----------
         os.system("emchain")
         os.system("epchain")
         # # ---------------------------------------------
         #------------rename file-----------------------
         os.system("rm *FIT.txt")
         v1=os.popen("ls *M1*EVLI*.FIT")
         mos1name=v1.read()[0:-1]
         v2=os.popen("ls *M2*EVLI*.FIT")
         mos2name=v2.read()[0:-1]
         v3=os.popen("ls *PI*EVLI*.FIT")
         pnname=v3.read()[0:-1]

         cmd="cp "+mos1name+" mos1.fits"
         os.system(cmd)
         cmd="cp "+mos2name+" mos2.fits"
         os.system(cmd)
         cmd="cp "+pnname+" pn.fits"
         os.system(cmd)
         # # ---------------------------------------------
         # #-----------barycen----------------------------
         for det in detList:
            cmd="cp "+det+".fits "+det+"_bary.fits"
            print(cmd)
            os.system(cmd)
            cmd="barycen withtable=yes table="+det+"_bary.fits"+": withsrccoordinates=yes srcra="+str(ra)+" srcdec="+str(dec)
            print(cmd)
            os.system(cmd)
         # # ---------------------------------------------
   #---------choose region-------------------
   #You should open pn.fits and choose the region, saved in physical coordinate.
   # # ---------------------------------------------
   #---------reg def------------------------------------
   # # ---------------------------------------------
   srcName = "QZ_VIR"
   srcReg = "circle(26354.489,27891.2,600.000)"
   bkgReg = "annulus(26354.489,27891.2,1200.000,2400.000)"

   if spectra:
      for det in detList:
         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)

         datapath = path+"/"+obsID+"/cal/"
         print(datapath)
         if det == "pn":
             cmd = "evselect table="+datapath+det+".fits withfilteredset=yes expression='(PATTERN <= 12)&&(PI in [200:15000])&&#XMMEA_EP' filteredset="+datapath+det+"_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
             cmd = "evselect table="+datapath+det+".fits withfilteredset=yes expression='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM' filteredset="+datapath+det+"_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         print(" ")
         print("1 filter by energy")
         print(cmd)
         os.system(cmd)
         cmd = "evselect table="+datapath+det+"_filt.fits withimageset=yes imageset="+datapath+det+"_filt.img xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600"
         print(" ")
         print("2 create det image")
         print(cmd)
         os.system(cmd)
         cmd = "evselect table='"+datapath+det+"_filt.fits:EVENTS' destruct=false withfilteredset=true withimageset=true imageset="+datapath+det+"_detmap.ds xcolumn=DETX ycolumn=DETY #withxranges=true ximagemin=-1500 ximagemax=1500 withyranges=true #yimagemin=-1500 yimagemax=1500 imagebinning='imageSize' #ximagesize=20 yimagesize=20 #writedss=true updateexposure=true"
         print(" ")
         print("3 create det map")
         print(cmd)
         os.system(cmd)
         if det == "pn":
             cmd = "evselect table='"+datapath+det+"_filt.fits:EVENTS' withspectrumset=yes spectrumset="+datapath+det+"_"+srcName+".pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='#XMMEA_EP && (PATTERN<=12) && ((X,Y) IN "+srcReg+")'"
         else:
             cmd = "evselect table='"+datapath+det+"_filt.fits:EVENTS' withspectrumset=yes spectrumset="+datapath+det+"_"+srcName+".pha energycolumn=PI spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN "+srcReg+")'"
         print(" ")
         print("4 extract source spectrum")
         print(cmd)
         os.system(cmd)
         if det == "pn":
             cmd = "evselect table='"+datapath+det+"_filt.fits:EVENTS' withspectrumset=yes spectrumset="+datapath+det+"_BKG_for"+srcName+".pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='#XMMEA_EP && (PATTERN<=12) && ((X,Y) IN "+bkgReg+")'"
         else:
             cmd = "evselect table='"+datapath+det+"_filt.fits:EVENTS' withspectrumset=yes spectrumset="+datapath+det+"_BKG_for"+srcName+".pha energycolumn=PI spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN "+bkgReg+")'"
         print(" ")
         print("5 extract background spectrum")
         print(cmd)
         os.system(cmd)

         cmd = "backscale spectrumset="+datapath+det+"_"+srcName+".pha badpixlocation="+datapath+det+"_filt.fits"
         print(" ")
         print("6 create source backscale keyword")
         print(cmd)
         os.system(cmd)
         cmd = "backscale spectrumset="+datapath+det+"_BKG_for"+srcName+".pha badpixlocation="+datapath+det+"_filt.fits"
         print(" ")
         print("7 create background backscale keyword")
         print(cmd)
         os.system(cmd)

         cmd = "rmfgen spectrumset="+datapath+det+"_"+srcName+".pha rmfset="+datapath+det+"_"+srcName+".rmf"
         print(" ")
         print("8 create rmf")
         print(cmd)
         os.system(cmd)
         cmd = "arfgen spectrumset="+datapath+det+"_"+srcName+".pha arfset="+datapath+det+"_"+srcName+".arf withrmfset=yes rmfset="+datapath+det+"_"+srcName+".rmf badpixlocation="+datapath+det+"_filt.fits extendedsource=yes detmaptype=dataset detmaparray="+datapath+det+"_detmap.ds"
         print(" ")
         print("9 create arf")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " + datapath+det+"_BKG_for"+srcName+".pha " +datapath+det+"_"+srcName+".pha " +"BACKFILE add=yes"
         print("10 add key")
         print(" ")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " + datapath+det+"_"+srcName+".rmf " +datapath+det+"_"+srcName+".pha " +"RESPFILE add=yes"
         print(" ")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " + datapath+det+"_"+srcName+".arf " + datapath + det + "_" + srcName + ".pha " + "ANCRFILE add=yes"
         print(" ")
         print(cmd)
         os.system(cmd)
         print(" ")
         print(" ")


import argparse
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

def calculdistance(line1,line2):
       dist=((((float(line2[8]))-(float(line1[8])))**2)+((float(line2[9]))-(float(line1[9])))**2+((float(line2[10]))-(float(line1[10])))**2)**(1/2)
       #dist=round(dist)
       return dist
       
def calculinterpolation(dist,pair,eng):
    score = pd.read_csv('score.csv',sep ='\t')
    for i in range (0,len(score[pair])-1):
      if i == int(math.floor(dist)):
       y1=float(score[pair][i])
       y2=float(score[pair][i+1])
       x1=float(math.floor(dist))
       x2=float(math.ceil(dist))
       e = y1 + (dist - x1) * (y2 - y1)/(x2 - x1)
       eng=eng+e
    return eng
    

def main(pdbfile):
#filein= open("4gxy.pdb","r")
 liste=[]
 for line in pdbfile:
     line1=[line[0:6],line[6:11],line[12:16],line[16:17],line[17:20],line[21:22],line[22:26],line[26:27],line[30:38],line[38:46],line[46:54]]
   
    
     if line1[2]==" C3'" and line1[0]== "ATOM  " and (line1[4]=="  A" or line1[4]=="  U" or    line1[4]=="  C" or line1[4]=="  G"):
    
      liste.append(line1)

 ll=[] 
 i=0
 while i<len(liste):
    r=1
    while r<len(liste):
  
 
     liste[i][6]=liste[i][6].replace(" ","")
  
     #print(liste[i][6])
     dd=int(liste[r][6])-int(liste[i][6])
     #print(dd)
     if dd>=4:
      #print(liste[i])
      #print(liste[r])
      ll.append(liste[i])
      ll.append(liste[r])
     r=r+1
    i=i+1
    r=1
 Gibbs_Energy=0   
 for i in range(0,len(ll)-1):
    dd=calculdistance(ll[i],ll[i+1])
    if dd<=20:
     pair=(ll[i][4].replace(" ","")+ll[i+1][4].replace(" ",""))
    
     if pair =='GC':
      pair=pair.replace('GC','CG')
     if pair =='GU':
      pair=pair.replace('GU','UG')
     if pair =='UA':
      pair=pair.replace('UA','AU')
     if pair =='GA':
      pair=pair.replace('GA','AG')
     if pair =='CA':
      pair=pair.replace('CA','AC')
     if pair =='CU':
      pair=pair.replace('CU','UC')

     Gibbs_Energy = calculinterpolation(dd, pair, Gibbs_Energy)

 print("Gibbs energy = ",Gibbs_Energy)       
      
      


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Python script")
    parser.add_argument("file",type=argparse.FileType('r'))
    args=parser.parse_args()
    main(args.file)





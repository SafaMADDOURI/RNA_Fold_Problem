import argparse
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate the distance
def calculdistance(line1,line2):
       dist=((((float(line2[8]))-(float(line1[8])))**2)+((float(line2[9]))-(float(line1[9])))**2+((float(line2[10]))-(float(line1[10])))**2)**(1/2)
       dist=round(dist)
       return dist

def freqobs(df,Nij):
  freq_obs={}
  for elt,t in zip(df.columns,Nij):
    fo=df[elt]/t
    freq_obs[elt] = fo
  freq_obs = pd.DataFrame.from_dict(freq_obs)
  freq_obs=freq_obs.transpose()	
  return freq_obs
  
def freqref(df,Nxx):  
  freq_ref={}
  for elt2,t2 in zip(df,Nxx):
    fr=df[elt2]/t2
    freq_ref[elt2] = fr
  freq_ref = pd.DataFrame.from_dict(freq_ref)	
  return freq_ref
  
def calculscore(freq_obs,freq_ref):   
  score={}
  for fo in freq_obs:
    sc= -1 * (np.log10(freq_obs[fo]/freq_ref[fo]))
    score[fo]=sc
  score = pd.DataFrame.from_dict(score)
  return score
  
def main(pdb_file):
 concatfile=open("concatfile","w") ##create a new file to concatenate all the pdb files
 args = parser.parse_args()
 for f in args.file:
  for line in f:
   concatfile.write(line)
  concatfile.write("\n") #Here, I created a file that englobe all the files that the user has input and the rest of the script will be executed on this file.
 ff=open('concatfile','r')  
 LISTEK=[] 
 for line in ff:
    line1=[line[0:6],line[6:11],line[12:16],line[16:17],line[17:20],line[21:22],line[22:26],line[26:27],line[30:38],line[38:46],line[46:54]]    
    if line1[2]==" C3'" and line1[0]== "ATOM  " and (line1[4]=="  A" or line1[4]=="  U" or    line1[4]=="  C" or line1[4]=="  G"):
      LISTEK.append(line1)   #list containig the lines that we need.
   
 ll=[] 
  #Here, I am going to build another list containing the lines (residues) separated by at least 3 positions,in the list it exisrs all the possible combinations of residues separated by at least 3 positions consecutively:
 #example 
 #line i / line i+1 ==> the 2 first residues separated by at least 3 positions
 #line i+2 / line i+3 ==> the 2 second residues separated by at least 3 positions
 i=0
 while i<len(LISTEK):
  r=1
  while r<len(LISTEK):
   LISTEK[i][6]=LISTEK[i][6].replace(" ","")
   dd=int(LISTEK[r][6])-int(LISTEK[i][6]) 
   if dd>=3:
    ll.append(LISTEK[i])
    ll.append(LISTEK[r])
   r=r+1
  i=i+1
  r=1

################### Calcul distance ############################################
##initialize the dictionnary with 0 distances for each pair
 dictdis={}
 listedistAA=[]
 listedistAU=[]
 listedistAG=[]
 listedistAC=[]
 listedistCG=[]
 listedistCC=[]
 listedistUU=[]
 listedistUG=[]
 listedistUC=[]
 listedistGG=[]
 for j in range(0,21):
  listedistAA.append(0)
  listedistGG.append(0)
  listedistAU.append(0)
  listedistAG.append(0)
  listedistAC.append(0)
  listedistCG.append(0)
  listedistCC.append(0)
  listedistUU.append(0)
  listedistUG.append(0)
  listedistUC.append(0)
 
 
 for i in range(0,len(ll)-1):
  if (ll[i][4] == "  A" and ll[i+1][4] == "  A" and ll[i][5]==ll[i+1][5]):
     dAA=calculdistance(ll[i],ll[i+1])
     if dAA<=21:
       #print(dAA)
       listedistAA[dAA-1]=listedistAA[dAA-1]+1

  if (ll[i][4] == "  G" and ll[i+1][4] == "  G" and ll[i][5]==ll[i+1][5]):
     dGG=calculdistance(ll[i],ll[i+1])
     if dGG<=21:
       listedistGG[dGG-1]=listedistGG[dGG-1]+1    
#print(listedistAA)

  if (ll[i][4] == "  U" and ll[i+1][4] == "  U" and ll[i][5]==ll[i+1][5]):
     dUU=calculdistance(ll[i],ll[i+1])
     if dUU<=21:
        listedistUU[dUU-1]=listedistUU[dUU-1]+1  
    
  if (ll[i][4] == "  C" and ll[i+1][4] == "  C" and ll[i][5]==ll[i+1][5]):
     dCC=calculdistance(ll[i],ll[i+1])
     if dCC<=21:
       listedistCC[dCC-1]=listedistCC[dCC-1]+1   
      
  if ((ll[i][4] == "  A" or ll[i][4] == "  G") and (ll[i+1][4] == "  G" or ll[i+1][4] == "  A") and ll[i][5]==ll[i+1][5]):
    dAG=calculdistance(ll[i],ll[i+1])
    if dAG<=21:
      listedistAG[dAG-1]=listedistAG[dAG-1]+1
      
  if ((ll[i][4] == "  A" or ll[i][4] == "  c") and (ll[i+1][4] == "  c" or ll[i+1][4] == "  A") and ll[i][5]==ll[i+1][5]):
    dAC=calculdistance(ll[i],ll[i+1])
    if dAC<=21:
       listedistAC[dAC-1]=listedistAC[dAC-1]+1
      
  if ((ll[i][4] == "  A" or ll[i][4] == "  U") and (ll[i+1][4] == "  U" or ll[i+1][4] == "  A") and ll[i][5]==ll[i+1][5]):
     dAU=calculdistance(ll[i],ll[i+1])
     if dAU<=21:
       listedistAU[dAU-1]=listedistAU[dAU-1]+1
      
  if ((ll[i][4] == "  U" or ll[i][4] == "  C") and (ll[i+1][4] == "  C" or ll[i+1][4] == "  U") and ll[i][5]==ll[i+1][5]):
     dUC=calculdistance(ll[i],ll[i+1])
     if dUC<=21:
       listedistUC[dUC-1]=listedistUC[dUC-1]+1
      
  if ((ll[i][4] == "  U" or ll[i][4] == "  G") and (ll[i+1][4] == "  G" or ll[i+1][4] == "  U") and ll[i][5]==ll[i+1][5]):
     dUG=calculdistance(ll[i],ll[i+1])
     if dUG<=21:
       listedistUG[dUG-1]=listedistUG[dUG-1]+1
      
  if ((ll[i][4] == "  C" or ll[i][4] == "  G") and (ll[i+1][4] == "  G" or ll[i+1][4] == "  C") and ll[i][5]==ll[i+1][5]):
     dCG=calculdistance(ll[i],ll[i+1])
     if dCG<=21:
       listedistCG[dCG-1]=listedistCG[dCG-1]+1   
  #Dictionnary for the distances
           
 dictdis["AA"]=listedistAA
 dictdis["GG"]=listedistGG
 dictdis["UU"]=listedistUU
 dictdis["CC"]=listedistCC
 dictdis["AG"]=listedistAG
 dictdis["AC"]=listedistAC
 dictdis["AU"]=listedistAU
 dictdis["UC"]=listedistUC
 dictdis["UG"]=listedistUG
 dictdis["CG"]=listedistCG
 #convert the dictopnnary into a dataframe
 df = pd.DataFrame(dictdis)
 dfdistance=df.transpose()
 #sum of all the columns in the dataframe
 dfdistance['sum_l'] = dfdistance.sum(axis=1)
  
 #sum of all rows as a new row in Dataframe
 sum = dfdistance.sum()
 sum.name = 'Sum_c'
 # Assign sum of all rows of DataFrame as a new Row
 dfdistance = dfdistance.append(sum.transpose()) #dataframe of the distances of each combinations with the sums of columns and rows.
  #print(dfdistance)
  
  
#This is another dataframe of the distances without the sum
 dfdistance2 = dfdistance.drop('sum_l', 1)
 dfdistance2=dfdistance2.drop(index='Sum_c')

 dfdistance2=dfdistance2.transpose()


  
  ##################Calcul score, observed freq and ref freq###########################
  
 Nij=[] ##List for the sum in all the lines 
 for i in range(0,len(dfdistance)):
  sum_l=dfdistance.iloc[i,21]
  Nij.append(sum_l)
 Nij.pop() 

 Nxx=[] ##List of the sum in the different columns
 for i in dfdistance.columns:
   sum_c=dfdistance.iloc[10][i]
   Nxx.append(sum_c)
 Nxx.pop()
  
 freq_obs=freqobs(dfdistance2,Nij)
 df3=dfdistance2.transpose()
 freq_ref=freqref(df3,Nxx)
 score=calculscore(freq_obs,freq_ref)
 for s in score.columns:        
    score.loc[score[s] >= 10, s] = 10
 score = score.fillna(10)        
 #print(score)
 score=score.transpose()
 print(score)
 score.to_csv('score.csv', sep='\t') # ==> Here, we will obtain the different scores of the different combinations in the different distances, all organized in a table (a csv file)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Python script")
    parser.add_argument("file",type=argparse.FileType('r'), nargs='+')
    args=parser.parse_args()
    main(args.file)




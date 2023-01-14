
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
def Plot(df,output_dir):
 for i in df.columns:
   if (i!="dist" and i!='0'):
      plt.figure()
      plt.plot(df["dist"],df[i], color='r',linewidth=2)
      plt.xlabel("Distance")
      plt.ylabel("Score")
      plt.title(i)                
      plt.savefig("{0}{1}".format(output_dir,i))
     
def main(score,output_dir):
 score = pd.read_csv(score,sep ='\t')
 score["dist"] = [i for i in range(0,21)]
 Plot(score,output_dir)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="plot")
    parser.add_argument("file",type=argparse.FileType('r'))
    parser.add_argument("-o","--output_dir",type=str,required=True)
    args=parser.parse_args()
    main(args.file,args.output_dir)



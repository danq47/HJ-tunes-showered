# Program to find the outlying st-2 files.

# POWHEG naively computes the arithmetic mean of the total cross sections
# However, if we have one rouge point, this will throw the average off by a lot
# Our solution is to compute the average weighted by 1/error**2
# Then, we can iterate over the list, removing the further outlying files
# until the POWHEG (arithmetic) mean and the weighted mean agree to within each others errors

import pandas as pd
import linecache    # Allows us to read a certain line of the file
import glob
import os

passed = False # if false, remove another file

while not passed:
   seeds=[]
   xsec=[]
   err=[]
   powheg_ave=[]
   running_pwg_err=[]
   weighted_ave=[]
   weighted_err=[]
   for filename in glob.iglob('pwg-st2-*-stat.dat'):
      xsec_and_err=linecache.getline(filename, 5) # this is the line containing the error and xsec
      xsec_and_err=xsec_and_err.split()         # Split on any whitespace 
      xsec.append(float(xsec_and_err[6]))       # The xsec
      err.append(float(xsec_and_err[8]))        # The error
      filename=filename.split("-")
      seeds.append(filename[2])              # The seed (from the filename)   
      
   d = {'Seed' : seeds , 'Xsec' : xsec , 'Error' : err}
   df = pd.DataFrame(d, columns=['Seed','Xsec','Error']) # Make a dataframe out of the these lists 

   y=0.0
   error=0.0
   ave=0.0
   err_pwg=0.0

   for i in range(df['Xsec'].size):    # Work out running averages
      ave = (df['Xsec'][i] + ave*i)/(i+1)
      y = y + df['Xsec'][i]/(df['Error'][i]**2)
      error = error + 1.0/(df['Error'][i]**2)
      weighted_ave.append(y/error)
      powheg_ave.append(ave)
      weighted_err.append(1.0/(error**0.5))  

      err_pwg = err_pwg + df['Error'][i]**2
      running_pwg_err.append((err_pwg/((i+1)**2))**0.5)
   

   df['POWHEG mean']=powheg_ave
   df['POWHEG running error']=running_pwg_err
   df['Weighted Average']=weighted_ave   
   df['Weighted Error']=weighted_err   

   # Check if it fulfils our conditions (i.e. are the two means equal (up to some error)?)
   normal = True     # Check that the POWHEG (arithmetic) mean is greater than the weighted average
   if powheg_ave[-1] > weighted_ave[-1]:
      normal = True
   else:
      normal = False 

   passed = False # Does our set of stat files pass our test?
   if normal:
      if powheg_ave[-1] > weighted_ave[-1] + weighted_err[-1]:
         passed = False
      else:
         passed = True  
   else:
      if powheg_ave[-1] < weighted_ave[-1] - weighted_err[-1]:
         passed = False
      else:
         passed = True
# If passed is true, we will just continue on
# If passed is false, we need to find the worst file, and remove it  

   if not passed:
      files = glob.glob('pwg-*-stat.dat')
      files.extend(glob.glob('pwggrid-*.dat'))
      files.extend(glob.glob('pwgcounters-*.dat'))
      files.extend(glob.glob('NLO-output/pwg-*-NLO.top'))
      seedToRemove=df[df['Error'] == df['Error'].max()]['Seed'].values[0]
      for filename in files:
         if seedToRemove in filename:
            os.rename(filename,filename+'-save')

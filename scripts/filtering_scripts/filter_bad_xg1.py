# Program to find the outlying xg1 files.
# We pick out the outliers so that the weighted mean and the arithmetic mean match

def is_Float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


import pandas as pd
import linecache    # Allows us to read a certain line of the file
import glob
import os

passed=[]
# if false, remove another file
passed.append(False)
passed.append(False)

while not passed[0] and not passed[1]:
   ini=True
   seeds=[]
   seeds_long=[]
   xsec=[]   # 0 will be abs, 1 will be rem
   err=[]
   powheg_ave=[]
   running_pwg_err=[]
   weighted_ave=[]
   weighted_err=[]
   xsec_locations=[]

   loop_list=[xsec,err,powheg_ave,running_pwg_err,weighted_ave,weighted_err]

   for _ in loop_list:   # Make them into a 2d array to have abs and rem
      _.append([])
      _.append([])

   for filename in glob.iglob('run-st1-xg1-*.log'):

# First we need to figure out which line has the numbers on it
      if ini:
         f = open(filename,'r')
         line_num=0
         while True:
            x = f.readline()
            if not x: 
               break
            line_num += 1
            x = x.split()
            if len(x) == 2:
               if is_Float(x[0]) and is_Float(x[1]):
                  xsec_locations.append(line_num)

      ini = False

      abs_xsec_and_err=linecache.getline(filename,xsec_locations[0])
      abs_xsec_and_err=abs_xsec_and_err.split()
      xsec[0].append(float(abs_xsec_and_err[0]))
      err[0].append(float(abs_xsec_and_err[1]))

      rem_xsec_and_err=linecache.getline(filename,xsec_locations[1])
      rem_xsec_and_err=rem_xsec_and_err.split()
      xsec[1].append(float(rem_xsec_and_err[0]))
      err[1].append(float(rem_xsec_and_err[1]))

      filename=filename.split("-")
      filename[3]=filename[3].split(".")
      seed=int(filename[3][0])
      seeds.append(str(seed))
      if seed > 999:
         seed=str(seed)
      elif seed > 99 and seed < 1000:
         seed = '0'+str(seed)
      elif seed > 9 and seed < 100:
         seed = '00'+str(seed)
      else:
         seed = '000'+str(seed)

      seeds_long.append(seed)

   d = {'Seed' : seeds_long, 'xsec_abs' : xsec[0], 'err_abs' : err[0], 'xsec_rem' : xsec[1], 'err_rem' : err[1]}
   df= pd.DataFrame(d,columns=['Seed','xsec_abs','err_abs','xsec_rem','err_rem'])

# Check the averages. the first index will be the absolute cross section, and the second will be the rem

   y=[0.0,0.0]
   error=[0.0,0.0]
   ave=[0.0,0.0]
   err_pwg=[0.0,0.0]
   xsec_ix=['xsec_abs','xsec_rem']
   err_ix=['err_abs','err_rem']


   for j in range(2):                           # work out for abs and rem
      for i in range(df['xsec_abs'].size):      # Work out running averages

         ave[j] = (df[ xsec_ix[j] ][i] + ave[j]*i)/(i+1)
         y[j] = y[j] + df[ xsec_ix[j] ][i]/(df[ err_ix[j] ][i]**2)
         error[j] = error[j] + 1.0/(df[ err_ix[j] ][i]**2)
         weighted_ave[j].append(y[j]/error[j])
         powheg_ave[j].append(ave[j])
         weighted_err[j].append(1.0/(error[j]**0.5))  

         err_pwg[j] = err_pwg[j] + df[ err_ix[j] ][i]**2
         running_pwg_err[j].append((err_pwg[j]/((i+1)**2))**0.5)
   

      df['POWHEG mean '+xsec_ix[j] ]=powheg_ave[j]
      df['POWHEG running '+err_ix[j] ]=running_pwg_err[j]
      df['Weighted Average '+xsec_ix[j] ]=weighted_ave[j]
      df['Weighted ' + err_ix[j] ]=weighted_err[j]

      # Check if it fulfils our conditions (i.e. are the two means equal (up to some error)?)
      normal=[]
      normal.append(True) # Check that the POWHEG (arithmetic) mean is greater than the weighted average
      normal.append(True)
      if powheg_ave[j][-1] > weighted_ave[j][-1]:
         normal[j] = True
      else:
         normal[j] = False 

         passed[0] = False
         passed[1] = False

      if normal[j]:
         if powheg_ave[j][-1] > weighted_ave[j][-1] + 2*weighted_err[j][-1]:
            passed[j] = False
         else:
            passed[j] = True  
      else:
         if powheg_ave[j][-1] < weighted_ave[j][-1] - 2*weighted_err[j][-1]:
            passed[j] = False
         else:
            passed[j] = True

# If passed is true, we will just continue on
# If passed is false, we need to find the worst file, and remove it  

   if not passed[0]: # the abs xsec has failed, remove the btl grid dat and top files
      seedToRemove=df[df['err_abs'] == df['err_abs'].max()]['Seed'].values[0]
      top='pwg-xg1-'+seedToRemove+'-btlgrid.top'
      dat='pwggridinfo-btl-xg1-'+seedToRemove+'.dat'
      log='run-st1-xg1-'+str(int(seedToRemove))+'.log'
      if os.path.isfile(top):
         os.rename(top,top+'-bad-xg1-abs')
      if os.path.isfile(dat):
         os.rename(dat,dat+'-bad-xg1-abs')
      if os.path.isfile(log): # may have already been deleted in the previous section from abs
         os.rename(log,log+'-bad-xg1-abs')
   
   if not passed[1]: # remnant has failed, however, both could fail so we put them in two separate loops
      seedToRemove=df[df['err_rem'] == df['err_rem'].max()]['Seed'].values[0]
      top='pwg-xg1-'+seedToRemove+'-rmngrid.top'
      dat='pwggridinfo-rmn-xg1-'+seedToRemove+'.dat'
      log='run-st1-xg1-'+str(int(seedToRemove))+'.log'
      if os.path.isfile(top):
         os.rename(top,top+'-bad-xg1-rem')
      if os.path.isfile(dat):
         os.rename(dat,dat+'-bad-xg1-rem')
      if os.path.isfile(log): # may have already been deleted in the previous section from abs
         os.rename(log,log+'-bad-xg1-rem')


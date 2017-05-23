#!/bin/bash
#########################################################
#########################################################
#########################################################
#
# Higgs mass to use:
#
hmass=125.5

#
# Square root of S hadronic
#
sroot=8000

#
# Beam particle types (protons, 1, or antiprotons, -1)
#
ih1=1
ih2=1

#
# PDF selection:
# See https://twiki.cern.ch/twiki/bin/view/LHCPhysics/PDF for suggested PDFs.
# MSTW2008nnlo68cl.LHgrid   0  1.0 1000000000.   1.E-06  1. MSTW2008-nnlo(central)
#
pdfset="MSTW2008nnlo68cl.LHgrid"
member=0

#
# Job prefix
#
prefix='HNNLO-LHC8-MSTW'

#
# Number of jobs to run for each mu_R and mu_F choice
#
njobs=10

#
# Number of iterations*points-to-throw in grid finding
# (itmx1*ncall1) & subsequent integration (itmx2*ncall2).
# itmx=15, ncall1=600000,   itmx2=30, ncall2=600000 runs in 31mins.
# itmx=15, ncall1=10000000, itmx2=30, ncall2=10000000 should then
# run in 8 1/2 hours 
#
itmx1=15
ncall1=10000000
itmx2=30
ncall2=8000000

#
# mu_R factors to scan over (mu_R=murfact*0.5*hmass)
#
murfacts="0.5 1.0 2.0"

#
# mu_F factors to scan over (mu_F=muffact*0.5*hmass)
#
muffacts="0.5 1.0 2.0"

#
# Initial random seed value to start looping from
#
initseed=10001

#
# Increment in random seeds between jobs
#
incrementseed=100

#########################################################
#########################################################
#########################################################

#########################################################
# Function to write scripts corresponding to each input #
# file to run on lxplus ...                             #
#########################################################
my_name=`whoami`
write_cmd_file() {
echo  "\

export LHAPATH="$LHAPATH"

touch "$PWD"/"$job_name".log

cd /tmp/
mkdir "$job_name"
cd "$job_name"

echo  \`( df -h 2>/dev/null || exit 0 )\` >> "$PWD"/"$job_name".log
echo  >> "$PWD"/"$job_name".log

cp "$PWD"/hnnlo ./
cp "$PWD"/br.sm1 ./
cp "$PWD"/br.sm2 ./
cp "$PWD"/"$job_name".input ./
mkdir PDFsets
cp "$PWD"/PDFsets/"$pdfset" ./PDFsets/

if [ -e "$PWD"/"$prefix-$smurfact$smuffact-$firstseed"_dvegas_virt.grid ] ; then
cp "$PWD"/"$prefix-$smurfact$smuffact-$firstseed"_dvegas_virt.grid ./
fi
if [ -e "$PWD"/"$prefix-$smurfact$smuffact-$firstseed"_dvegas_real.grid ] ; then
cp "$PWD"/"$prefix-$smurfact$smuffact-$firstseed"_dvegas_real.grid ./
fi

ls    >> "$PWD"/"$job_name".log

echo  >> "$PWD"/"$job_name".log
echo  >> "$PWD"/"$job_name".log

if [ -e "$prefix-$smurfact$smuffact-$firstseed"_dvegas_virt.grid ] ; then
mv "$prefix-$smurfact$smuffact-$firstseed"_dvegas_virt.grid "$job_name"_dvegas_virt.grid
fi
if [ -e "$prefix-$smurfact$smuffact-$firstseed"_dvegas_real.grid ] ; then
mv "$prefix-$smurfact$smuffact-$firstseed"_dvegas_real.grid "$job_name"_dvegas_real.grid
fi

chmod u+rwx hnnlo
chmod u+rwx *.*

ls    >> "$PWD"/"$job_name".log

./hnnlo < "$job_name".input >> "$PWD"/"$job_name".log

cp *.*   "$PWD"/
cd ..

echo \"Listing run directory BEFORE removal:\" >> "$PWD"/"$job_name".log
ls "$job_name" >> "$PWD"/"$job_name".log
echo \"Removing directory "$job_name" \"  >> "$PWD"/"$job_name".log
rm -rf "$job_name" 

# echo \"Listing run directory AFTER removal:\" >> "$PWD"/"$job_name".log
# ls "$job_name" >> "$PWD"/"$job_name".log

" >> $1.sh
}

#########################################################
#########################################################
#########################################################

# Set LHAPDF pdf sets path
LHAPATH=`lhapdf-config --pdfsets-path`

#
# Make template from the standard template ...
#
if [ -e template.input ] ; then rm template.input ; fi
cp HNNLO.input template.input

#
# Set the Higgs mass in the template file.
#
echo ""
echo "Setting Higgs mass to $hmass in template file ..."
sed -i -e "s/.*\! hmass/$hmass   \! hmass/" template.input
# grep hmass template.input

#
# Set the sqrt(S) in the template file.
#
echo ""
echo "Setting sqrt(S) $sroot in template file ..."
sed -i -e "s/.*\! sroot/$sroot   \! sroot/" template.input
# grep sroot template.input

#
# Set the beam particles in the template file.
#
echo ""
echo "Setting ih1=$ih1 and $ih2=$ih2 in template file ..."
sed -i -e "s/.*\! ih1 ih2/$ih1 $ih2  \! ih1 ih2/" template.input
# grep ih1 template.input

#
# Set the PDF set and member.
#
echo ""
echo "Setting PDF set=$pdfset and member=$member in template file ..."
sed -i -e "s/.*\! set, member (LHAPDFs)/\'$pdfset\' $member  \! set, member (LHAPDFs)/" template.input
# grep ih1 template.input

#
# Number of iterations*points-to-throw in grid finding
# (itmx1*ncall1) & subsequent integration (itmx2*ncall2).
#
echo ""
echo "Setting itmx1=$itmx1, ncall1=$ncall1 and "
echo "        itmx2=$itmx2, ncall2=$ncall2 ..."
sed -i -e "s/.*\! itmx1, ncall1/$itmx1 $ncall1  \! itmx1, ncall1/" template.input
sed -i -e "s/.*\! itmx2, ncall2/$itmx2 $ncall2  \! itmx2, ncall2/" template.input
# grep ih1 template.input

#
# Loop over muR and muF values and also the random seeds
# copying the template.input to $prefix-$smurfact$smuffact-$seed.input
# and then editing the result in each iteration of the loop to set
# muR, muF, the random seed and the runstrong appropriately.
#
for murfact in $murfacts
do
  for muffact in $muffacts
  do
    case $murfact in 0.25) smurfact=Q ;; 0.5) smurfact=H ;; 1.0) smurfact=1 ;; 2.0) smurfact=2 ;; esac
    case $muffact in 0.25) smuffact=Q ;; 0.5) smuffact=H ;; 1.0) smuffact=1 ;; 2.0) smuffact=2 ;; esac
    mur=`echo "scale=10; 0.5*$hmass*$murfact" | bc`
    muf=`echo "scale=10; 0.5*$hmass*$muffact" | bc`
#    echo $murfact $mur $muffact $muf
    for ((ixx=1;ixx<=njobs;ixx++))
    do seed=$((ixx*incrementseed+initseed))
      if [ "$ixx" = "1" ] ; then
	firstseed=$seed
	if [ -e $prefix-$smurfact$smuffact-$firstseed\_dvegas_virt.grid ] ; then
	   echo "Found existing virtual grid: "$prefix-$smurfact$smuffact-$firstseed\_dvegas_virt.grid 
           echo "All $prefix-$smurfact$smuffact jobs will use this!"
	fi
	if [ -e $prefix-$smurfact$smuffact-$firstseed\_dvegas_real.grid ] ; then
	   echo "Found existing real grid: "$prefix-$smurfact$smuffact-$firstseed\_dvegas_real.grid 
           echo "All $prefix-$smurfact$smuffact jobs will use this!"
	fi
      fi
#      echo $seed
      job_name=$prefix-$smurfact$smuffact-$seed
      inputfile=$job_name.input
      cp template.input $inputfile
      sed -i -e "s/.*\! mur,muf/$mur $muf \! mur,muf/" $inputfile
      sed -i -e "s/.*\! rseed/$seed \! rseed/" $inputfile
      sed -i -e "s/.*\! runstring/\'$job_name\' \! runstring/" $inputfile
      if [ -e $job_name.sh ] ; then rm $job_name.sh ; fi
      write_cmd_file $job_name
    done
  done
done

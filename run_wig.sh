#!/bin/sh
./cl_wig
TD=` grep NTLOOP input00 | sed -e 's/NTLOOP=//' `
echo $TD
PARAM="NTLOOP"
i=0  # $PARAM=i
PARAM2="TEMP"
j=1.0  # $PARAM2=j

if [ $TD = 3 ]
then
PARAM3="GANMA"
k=0.10 # $PARAM3=k

for k in `cat glist`
do
for i in 3
do
for j in 2.0
do 
  sed -e "s/$PARAM=.*/$PARAM=$i,/g" input00 | sed -e "s/$PARAM2=.*/$PARAM2=$j,/g" | sed -e "s/$PARAM3=.*/$PARAM3=$k,/g" > input
./FRACTION.exe
done
done

rm fort.*
done
mkdir gamma
mv GCR gamma/
fi

if [ $TD = 2 ]
then
mkdir nloop
PARAM3="NA"
k=1000 # $PARAM3=k
PARAM4="NB"
l=1000 # $PARAM4=l
for i in 2
do
for j in `cat tlist`
do 
for k in `cat pplist`
do
for l in `cat pslist`
do
 sed -e "s/$PARAM=.*/$PARAM=$i,/g" input00 | sed -e "s/$PARAM2=.*/$PARAM2=$j,/g" | sed -e "s/$PARAM3=.*/$PARAM3=$k,/g" | sed -e "s/$PARAM4=.*/$PARAM4=$l,/g" > input
./FRACTION.exe
done
cat fort.720 >> fort.600
cat fort.720 >> fort.607
cat fort.720 >> fort.608
cat fort.720 >> fort.605
cat fort.720 >> fort.606
done
mv fort* nloop
mv nloop/fort.600 nloop/CR_K_"$j"
mv nloop/fort.607 nloop/CR_mol_"$j"
mv nloop/fort.608 nloop/CR_raw_"$j"
mv nloop/fort.605 nloop/CB_"$j"
mv nloop/fort.606 nloop/CT_"$j"
done
done

fi

if [ $TD = 1 ]
then

for i in 1
do
for j in 2.0
do 
  sed -e "s/$PARAM=.*/$PARAM=$i,/g" input00 | sed -e "s/$PARAM2=.*/$PARAM2=$j,/g"  > input
./FRACTION.exe
done
done

mkdir data
rm fort.9940
mv fort* data
mv data/fort.600 data/CR
mv data/fort.7100 data/lf
mv data/fort.601 data/PSPD_K
mv data/fort.602 data/PSPD_RB
mv data/fort.603 data/PSPD_BEC
mv data/fort.604 data/PSPD_THE
mv data/fort.605 data/convert_BEC
mv data/fort.606 data/convert_THE
mv data/fort.607 data/PB
cd data
gnuplot <<EOF
load "lf"
EOF

fi

if  [ $TD = 0 ]
then

#mkdir ps_contour

for i in 0
do
for j in `cat tlist`
do 
  sed -e "s/$PARAM=.*/$PARAM=$i,/g" input00 | sed -e "s/$PARAM2=.*/$PARAM2=$j,/g"  > input
./FRACTION.exe
mkdir data_T="$j"
mv fort* data_T="$j"
mv data_T="$j"/fort.3001 data_T="$j"/RPK_"$j"
mv data_T="$j"/fort.3002 data_T="$j"/RPRB_"$j"
mv data_T="$j"/fort.3003 data_T="$j"/RPK_AFT_"$j"
mv data_T="$j"/fort.3004 data_T="$j"/RPRB_AFT_"$j"
mv data_T="$j"/fort.3005 data_T="$j"/RPK_R_"$j"
mv data_T="$j"/fort.3006 data_T="$j"/RPRB_R_"$j"
mv data_T="$j"/fort.3007 data_T="$j"/RPK_R_AFT_"$j"
mv data_T="$j"/fort.3008 data_T="$j"/RPRB_R_AFT_"$j"
mv data_T="$j"/fort.4001 data_T="$j"/KPX_"$j"
mv data_T="$j"/fort.4002 data_T="$j"/KPY_"$j"
mv data_T="$j"/fort.4003 data_T="$j"/KPZ_"$j"
mv data_T="$j"/fort.4101 data_T="$j"/KPX_AFT_"$j"
mv data_T="$j"/fort.4102 data_T="$j"/KPY_AFT_"$j"
mv data_T="$j"/fort.4103 data_T="$j"/KPZ_AFT_"$j"
mv data_T="$j"/fort.5001 data_T="$j"/RBPX_"$j"
mv data_T="$j"/fort.5002 data_T="$j"/RBPY_"$j"
mv data_T="$j"/fort.5003 data_T="$j"/RBPZ_"$j"
mv data_T="$j"/fort.5101 data_T="$j"/RBPX_AFT_"$j"
mv data_T="$j"/fort.5102 data_T="$j"/RBPY_AFT_"$j"
mv data_T="$j"/fort.5103 data_T="$j"/RBPZ_AFT_"$j"
mv data_T="$j"/fort.4004 data_T="$j"/MXK_"$j"
mv data_T="$j"/fort.4005 data_T="$j"/MYK_"$j"
mv data_T="$j"/fort.4006 data_T="$j"/MZK_"$j"
mv data_T="$j"/fort.4104 data_T="$j"/MXK_AFT_"$j"
mv data_T="$j"/fort.4105 data_T="$j"/MYK_AFT_"$j"
mv data_T="$j"/fort.4106 data_T="$j"/MZK_AFT_"$j"
mv data_T="$j"/fort.5004 data_T="$j"/MXRB_"$j"
mv data_T="$j"/fort.5005 data_T="$j"/MYRB_"$j"
mv data_T="$j"/fort.5006 data_T="$j"/MZRB_"$j"
mv data_T="$j"/fort.5104 data_T="$j"/MXRB_AFT_"$j"
mv data_T="$j"/fort.5105 data_T="$j"/MYRB_AFT_"$j"
mv data_T="$j"/fort.5106 data_T="$j"/MZRB_AFT_"$j"
mv data_T="$j"/fort.6001 data_T="$j"/KRBPX_"$j"
mv data_T="$j"/fort.6002 data_T="$j"/KRBPY_"$j"
mv data_T="$j"/fort.6003 data_T="$j"/KRBPZ_"$j"
mv data_T="$j"/fort.6004 data_T="$j"/MXKRB_"$j"
mv data_T="$j"/fort.6005 data_T="$j"/MYKRB_"$j"
mv data_T="$j"/fort.6006 data_T="$j"/MZKRB_"$j"
mv data_T="$j"/fort.7001 data_T="$j"/KEN_re_"$j"
mv data_T="$j"/fort.7003 data_T="$j"/RBEN_re_"$j"
mv data_T="$j"/fort.7002 data_T="$j"/KEN_"$j"
mv data_T="$j"/fort.7004 data_T="$j"/RBEN_"$j"
mv data_T="$j"/fort.9940 data_T="$j"/DRDP_"$j"
mv data_T="$j"/fort.9941 data_T="$j"/DRDIST_"$j"
mv data_T="$j"/fort.9942 data_T="$j"/DPDIST_"$j"
#cd ps_contour
#gnuplot <<EOF
#set contour
#set cntrparam levels 5
#unset surface
#set view 0,0
#set term table
#set output "kps.dat"
#sp "Kps"
#set output "rbps.dat"
#sp "RBps"
#set term postscript eps enhanced color
#set output"contour.eps"
#set xr [-4e-5:4e-5]
#set yr [-7e-28:7e-28]
#set xlabel "z position"
#set ylabel "z momentum"
#plot "kps.dat" u 1:2 w l lt 1 ti "K","rbps.dat" u 1:2 w l lt 2 ti "Rb"
#EOF
#mv contour.eps contour_"$j".eps
#cd ../
done
done

mkdir POSITION
cp data_T=*/*PX_* POSITION/
cp data_T=*/*PY_* POSITION/
cp data_T=*/*PZ_* POSITION/
mkdir MOMENTUM
cp data_T=*/MX*_* MOMENTUM/
cp data_T=*/MY*_* MOMENTUM/
cp data_T=*/MZ*_* MOMENTUM/
mkdir ENERGY
cp data_T=*/*EN*_* ENERGY/
mkdir RPSPACE
cp data_T=*/D*D* RPSPACE/
cp data_T=*/RP* RPSPACE/

fi
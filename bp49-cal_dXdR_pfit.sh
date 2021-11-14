#command: sh cal_dXdRzz_pfit.sh num_atoms -4 -3 -2 -1 0 1 2 3 4
#-----------------------------------------------------------------------#
gfortran -c PolynomialFitting.f90
gfortran -o polyfit.x PolynomialFitting.o -L$LAPACK_LIB_DIR -llapack -lblas
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
nat=$1
h1=$2
h2=$3
h3=$4
h4=$5
h5=$6
h6=$7
h7=$8
h8=$9
shift 
h9=$9
shift
h10=$9
shift
h11=$9
shift
h12=$9
shift
h13=$9
shift
h14=$9
shift
h15=$9
shift
h16=$9
shift
h17=$9


#parameter for unit convert
E=0.001
F_to_eV=25.711057
c=36.3609     #  1 Ry.a.u. = 36.3609 V/Angstrom  USE FOR BERRY PHASE CALCULATION
#c=51.4220632   #  1  H.a.u. = 51.4220632 V/Angstrom USE FOR SAWTOOTH POTENTIAL CALN
c1=14.40036    #  1 e.Angstrom^2/V = 14.40036 Angstrom^3
c0=$(echo "$F_to_eV*$c1/($c^2)" | bc -l )


nat1=$(echo $nat+1 | bc )



#------------------xx
#------------------

pref=../E-Ex
pref2=-Ey0-Ez0
#-------------------------------------------------------------------------
for h in  $h1 $h2 $h3 $h4 $h5 $h6 $h7 $h8 $h9 $h10 $h11 $h12 $h13 $h14 $h15 $h16 $h17; do
e=$(echo $h*$E | bc)
echo "$e" > xval$h
grep -A$nat1 Forces $pref$h$pref2.out | tail -$nat | awk '{printf"  %3.9f     %3.9f     %3.9f\n",$(NF-2),$(NF-1),$(NF)}' > F_Ez$h
done

grep -A$nat1 Forces $pref$h$pref2.out | tail -$nat  | awk '{printf"  %3s\n", $(NF-5) }' > list_spec_index.dat

#--------------------------------------
nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=1;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr1


nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=2;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr2

nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=3;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr3
#--------------------------------------
#echo "1 1" >> dXidR.dat
paste dxdr1 dxdr2 dxdr3  | awk '{printf "  %3.9f     %3.9f     %3.9f \n", $1*2*'$c0', $2*2*'$c0', $3*2*'$c0'}' > xx-tmp.dat
#paste list_spec_index.dat tmp.dat >> dXidR.dat
#--------------------------------------------

rm xval*  dxdr* d2ydx2.out F_Ez* inp list_spec_index.dat





#------------------yy
#------------------


pref=../E-Ex0-Ey
pref2=-Ez0
#-------------------------------------------------------------------------
for h in  $h1 $h2 $h3 $h4 $h5 $h6 $h7 $h8 $h9 $h10 $h11 $h12 $h13 $h14 $h15 $h16 $h17; do
e=$(echo $h*$E | bc)
echo "$e" > xval$h
grep -A$nat1 Forces $pref$h$pref2.out | tail -$nat | awk '{printf"  %3.9f     %3.9f     %3.9f\n",$(NF-2),$(NF-1),$(NF)}' > F_Ez$h
done

grep -A$nat1 Forces $pref$h$pref2.out | tail -$nat  | awk '{printf"  %3s\n", $(NF-5) }' > list_spec_index.dat

#--------------------------------------
nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=1;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr1


nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=2;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr2

nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=3;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr3
#--------------------------------------
#echo "2 2" >> dXidR.dat
paste dxdr1 dxdr2 dxdr3  | awk '{printf "  %3.9f     %3.9f     %3.9f \n", $1*2*'$c0', $2*2*'$c0', $3*2*'$c0'}' > yy-tmp.dat
#paste list_spec_index.dat tmp.dat >> dXidR.dat
#--------------------------------------------

rm xval*  dxdr* d2ydx2.out inp list_spec_index.dat F_Ez*



#------------------zz
#------------------

pref=../E-Ex0-Ey0-Ez
#-------------------------------------------------------------------------
for h in  $h1 $h2 $h3 $h4 $h5 $h6 $h7 $h8 $h9 $h10 $h11 $h12 $h13 $h14 $h15 $h16 $h17; do
e=$(echo $h*$E | bc)
echo "$e" > xval$h
grep -A$nat1 Forces $pref$h.out | tail -$nat | awk '{printf"  %3.9f     %3.9f     %3.9f\n",$(NF-2),$(NF-1),$(NF)}' > F_Ez$h
done

grep -A$nat1 Forces $pref$h.out | tail -$nat  | awk '{printf"  %3s\n", $(NF-5) }' > list_spec_index.dat

#--------------------------------------
nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=1;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr1


nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=2;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr2

nf=$(ls F_Ez* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Ez*   | awk '{ for (i=3;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr3
#--------------------------------------
#echo "3 3" >> dXidR.dat
paste dxdr1 dxdr2 dxdr3  | awk '{printf "  %3.9f     %3.9f     %3.9f \n", $1*2*'$c0', $2*2*'$c0', $3*2*'$c0'}' > zz-tmp.dat
#paste list_spec_index.dat tmp.dat >> dXidR.dat
#--------------------------------------------

rm xval*  dxdr* d2ydx2.out inp list_spec_index.dat F_Ez*




#------------------xy
#------------------

pref=../E-Ex
pref2=-Ey
pref3=-Ez0
#-------------------------------------------------------------------------
for h in  $h1 $h2 $h3 $h4 $h5 $h6 $h7 $h8 $h9 $h10 $h11 $h12 $h13 $h14 $h15 $h16 $h17; do
e=$(echo $h*$E | bc)
echo "$e" > xval$h
grep -A$nat1 Forces $pref$h$pref2$h$pref3.out | tail -$nat | awk '{printf"  %3.9f     %3.9f     %3.9f\n",$(NF-2),$(NF-1),$(NF)}' > F_Exy$h
done

grep -A$nat1 Forces $pref$h$pref2$h$pref3.out | tail -$nat  | awk '{printf"  %3s\n", $(NF-5) }' > list_spec_index.dat

#--------------------------------------
nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=1;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr1


nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=2;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr2

nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=3;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr3

#--------------------------------------
#echo "1 2" >> dXidR.dat
paste dxdr1 dxdr2 dxdr3  | awk '{printf "  %3.9f     %3.9f     %3.9f \n", $1*2*'$c0', $2*2*'$c0', $3*2*'$c0'}' > xy_tmp.dat
#paste list_spec_index.dat tmp.dat >> dXidR.dat
#--------------------------------------------

rm xval* F_Exy* dxdr* d2ydx2.out inp list_spec_index.dat 



#------------------yz
#------------------

pref=../E-Ex0
pref2=-Ey
pref3=-Ez
#-------------------------------------------------------------------------
for h in  $h1 $h2 $h3 $h4 $h5 $h6 $h7 $h8 $h9 $h10 $h11 $h12 $h13 $h14 $h15 $h16 $h17; do
e=$(echo $h*$E | bc)
echo "$e" > xval$h
grep -A$nat1 Forces $pref$pref2$h$pref3$h.out | tail -$nat | awk '{printf"  %3.9f     %3.9f     %3.9f\n",$(NF-2),$(NF-1),$(NF)}' > F_Exy$h
done

grep -A$nat1 Forces $pref$pref2$h$pref3$h.out | tail -$nat  | awk '{printf"  %3s\n", $(NF-5) }' > list_spec_index.dat

#--------------------------------------
nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=1;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr1


nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=2;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr2

nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=3;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr3

#--------------------------------------
#echo "2 3" >> dXidR.dat
paste dxdr1 dxdr2 dxdr3  | awk '{printf "  %3.9f     %3.9f     %3.9f \n", $1*2*'$c0', $2*2*'$c0', $3*2*'$c0'}' > yz_tmp.dat
#paste list_spec_index.dat tmp.dat >> dXidR.dat
#--------------------------------------------

rm xval*  dxdr* d2ydx2.out inp list_spec_index.dat F_Exy*


#------------------zx
#------------------

pref=../E-Ex
pref2=-Ey0
pref3=-Ez
#-------------------------------------------------------------------------
for h in  $h1 $h2 $h3 $h4 $h5 $h6 $h7 $h8 $h9 $h10 $h11 $h12 $h13 $h14 $h15 $h16 $h17; do
e=$(echo $h*$E | bc)
echo "$e" > xval$h
grep -A$nat1 Forces $pref$h$pref2$pref3$h.out | tail -$nat | awk '{printf"  %3.9f     %3.9f     %3.9f\n",$(NF-2),$(NF-1),$(NF)}' > F_Exy$h
done

grep -A$nat1 Forces $pref$h$pref2$pref3$h.out | tail -$nat  | awk '{printf"  %3s\n", $(NF-5) }' > list_spec_index.dat

#--------------------------------------
nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=1;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr1


nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=2;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr2

nf=$(ls F_Exy* | wc -l )
echo "$nf  $nat" > inp
paste xval* >> inp

paste F_Exy*   | awk '{ for (i=3;i<=NF;i+=3) printf"  %3.9f ", $i ; printf "\n"}' >> inp
./polyfit.x
cat d2ydx2.out > dxdr3

#--------------------------------------
#echo "3 1" >> dXidR.dat
paste dxdr1 dxdr2 dxdr3  | awk '{printf "  %3.9f     %3.9f     %3.9f \n", $1*2*'$c0', $2*2*'$c0', $3*2*'$c0'}' > zx_tmp.dat
#paste list_spec_index.dat tmp.dat >> dXidR.dat
#--------------------------------------------

rm xval*  dxdr* d2ydx2.out F_Exy*






#----------------get dXiRi  h=ei=ej
#----------------equation: d2F/deidej = 0.5*(d2F/dh2-d2F/dei2-d2F/dej2)
echo "1 1" > dXidR.dat
paste list_spec_index.dat xx_tmp.dat >> dXidR.dat

echo "1 2" >> dXidR.dat
paste xy_tmp.dat xx_tmp.dat yy_tmp.dat |  awk '{printf "  %3.9f     %3.9f     %3.9f \n", (($1)+(-$4)+(-$7))/2, (($2)+(-$5)+(-$8))/2,(($3)+(-$6)+(-$9))/2}' > tmp1.dat
paste list_spec_index.dat tmp1.dat >> dXidR.dat

echo "1 3" >> dXidR.dat
paste zx_tmp.dat xx_tmp.dat zz_tmp.dat |  awk '{printf "  %3.9f     %3.9f     %3.9f \n", (($1)+(-$4)+(-$7))/2, (($2)+(-$5)+(-$8))/2,(($3)+(-$6)+(-$9))/2}' > tmp1.dat
paste list_spec_index.dat tmp1.dat >> dXidR.dat

echo "2 1" >> dXidR.dat
paste xy_tmp.dat xx_tmp.dat yy_tmp.dat |  awk '{printf "  %3.9f     %3.9f     %3.9f \n", (($1)+(-$4)+(-$7))/2, (($2)+(-$5)+(-$8))/2,(($3)+(-$6)+(-$9))/2}' > tmp1.dat
paste list_spec_index.dat tmp1.dat >> dXidR.dat

echo "2 2" >> dXidR.dat
paste list_spec_index.dat yy_tmp.dat >> dXidR.dat

echo "2 3" >> dXidR.dat
paste yz_tmp.dat yy_tmp.dat zz_tmp.dat |  awk '{printf "  %3.9f     %3.9f     %3.9f \n", (($1)+(-$4)+(-$7))/2, (($2)+(-$5)+(-$8))/2,(($3)+(-$6)+(-$9))/2}' > tmp1.dat
paste list_spec_index.dat tmp1.dat >> dXidR.dat

echo "3 1" >> dXidR.dat
paste zx_tmp.dat xx_tmp.dat zz_tmp.dat |  awk '{printf "  %3.9f     %3.9f     %3.9f \n", (($1)+(-$4)+(-$7))/2, (($2)+(-$5)+(-$8))/2,(($3)+(-$6)+(-$9))/2}' > tmp1.dat
paste list_spec_index.dat tmp1.dat >> dXidR.dat

echo "3 2" >> dXidR.dat
paste yz_tmp.dat yy_tmp.dat zz_tmp.dat |  awk '{printf "  %3.9f     %3.9f     %3.9f \n", (($1)+(-$4)+(-$7))/2, (($2)+(-$5)+(-$8))/2,(($3)+(-$6)+(-$9))/2}' > tmp1.dat
paste list_spec_index.dat tmp1.dat >> dXidR.dat

echo "3 3" >> dXidR.dat
paste list_spec_index.dat zz_tmp.dat >> dXidR.dat


rm list_spec_index.dat tmp1.dat *_tmp.dat





# Non-resonant-Raman-finitedifference
How to use codes?
--------------------
1. This code is based on the method in ref [10.1016/j.diamond.2004.12.007]. So to compute Raman tensor, we need to compute forces on atoms for systems without electric field, with positive electric field, and with negative electric field. Phonon eiegenvectors are needed.
2. Computation of dXidR. A polynomial fitting is performed using the bash file "bp49-cal_dXdR_pfit.sh". Read it and modify accordingly. It compute the second derivative of force with respect to electric field. 
3. Computation of Raman spectra. After step 2, we can get dXidR (full tensor: xx,yy,zz,xy,xz,yz). Then using the code "RAMAN-CODE.f90" to get the Raman spectra. Read the head of the code to get what input files you need.

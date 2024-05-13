These files are from Eric Christy <christy@jlab.org>.

The ep elastic subroutine is called 'elastic.f' and contains several
parameterizations, including John Arrington's fit of the cross section data,
which also describes the RESR data quite well.  The default parameterization is
good for 0.5<Q^2<8 GeV^2.
 
Run the shell script 'makerates' to compile.

Input the current in uA.  Out rate in kHz is for HMS, using 4 cm LH2 target.

The code is good for inclusive scattering so it will first check whether
elastic kinematics are in the acceptance or not and simply average across the
theta acceptance if it is.  It will then integrate the inelastic cross section
across the sprectrometer acceptance and add that contribution as well.  There
is a switch for selecting the spectrometer (the logical sos) in the main code
which will then set the correct effective momentum acceptace (this is currently
set to ?8% for the HMS and should be ~+/-11%).

Regards

-Eric

---

Type: make 
to compile and create the executable from the source FORTRAN code.

To run the program, put your inputs in a .dat file in the order:
  E_Beam E_spec TH_spec Beam_uA
specifying HMS or SHMS in the top line.  You can customize the solid angle 
and momentum bite of the spectrometer in rates.f

Then type: rates
and it will ask you for the name of the input file.  Two example input 
files are given.

A file: rates.out is created.  Move it to filename.out, where filename 
is the same as the input file, to keep things organized.

GH

#!/usr/bin/perl -w
use strict;

#
# This example script generates ini files 
# and condor submission files for CBwaves
# to be submitted to Condor clusters.
#
# Modify it accordin to your needs.
#
# Set $generatesubfileonly = "no" below if 
# you want to run cbwaves on the fly.
#

# Definition of auxiliary variables

my $PI     = 3.141592653589793;
my $SI_c   = 2.99792458E8;
my $SI_G   = 6.67428E-11;
my $SI_Gc2 = 7.42613801611754e-28;
my $SI_ly  = 9460730472580800.0;
my $SI_pc  = 3.26156*$SI_ly;


# - Where is the cbwaves executable
my $cbwbinary="/Users/b/Documents/CppProjects/cbwaves-verziok/cbwaves-1.2.0-2/build/CBwaves";

# - Output file name prefix

my $filenameprefix="test";

# - Mass ratio: 10:1.4
#       m1 = 10*G/c^2*2e30 = 14850 meters
#	m2 = 1.4*G/c^2*2e30 = 2079 meters

my $msun = 1477.;
my $m1 = 10*$msun;
my $m2 = 10*$msun;

my $z = 0.3;
my $M = $m1 + $m2;

# - Initial separation: r = 15*(m1+m2)
my $r0 = 300;
my $r = $r0*($M*(1+$z));
my $ri = 300*($M*(1+$z));

# - Time step
#my $dt = 100000.;
my $dt = 1./4096.;
#my $dt = 1./16384.;

# - Maximum evolution time
my $tmax = 3.0e20;

# - Maximum number of orbits
my $orbitsmax = 1.0e8;
#my $orbitsmax = 203;
#my $orbitsmax = 6;


# - Eccentricity of the orbit
#my $epsilon = 0.6693;
my $epsilon = 0.0;

# - Polar angles
my $iota    = 0.0;
my $phi     = 0.0;
my $theta   = 0.0;
my $varphi  = 0.0;
my $psi     = 0.0;

# - Spin definition
my $s1 = 0;
my $s2 = 0;
#my ($s1x, $s1y, $s1z) = (-sin(0.1020645929),0,cos(0.1020645929));
#my ($s2x, $s2y, $s2z) = (sin(0.4215341828),0,cos(0.4215341828));
my ($s1x, $s1y, $s1z) = (0,0,0);
my ($s2x, $s2y, $s2z) = (0,0,0);

	# $s1x = sin(56.3/180 * $PI);
	# $s1y = cos(56.3/180 * $PI);
	# $s1z = 0;
	# $s2x = sin(56.3/180 * $PI);
	# $s2y = cos(56.3/180 * $PI);
	# $s2z = 0;


# - Output filename definition
my $outfile="${filenameprefix}.dat";
my $ftfile="/tmp/cbwaves.ft";

# - Hterms and corrections
my $hterms= "'Q','P05Q','PQ','P15Q','P15Qtail','PQSO','P15QSO','P2Q','PQSS'";
#my $corrs="'PN','2PN','SO','SS','RR','PNSO','3PN','1RR','2PNSO','RRSO','RRSS'";
my $corrs = "'PN','2PN','3PN','4PN'";
# my $corrs = "'PN','2PN','3PN','4PN','SO','SS','PNSO','2PNSO','1RR','RRSO','RRSS'";

# - Output variables
# my $outvars = "t,orbits,rx,ry,rz,h_+,E_N,E_PNtot,E_tot,E_rad,hp22,hx22,x1,y1,z1,x2,y2,z2,hp22,hp21,hp20,hp2m1,hp2m2,hx22,hx21,hx20,hx2m1,hx2m2,h_x,h";
my $outvars = "t,x1,y1,z1,x2,y2,z2";
#my $outvars = "t,x1,y1,z1,x2,y2,z2,E,E_tot,E_PNtot,h,h_x,h_+";



# - Do we want checkpointing 
my $checkpoint = "no";

# - Do we want initial eccentricity approximation 
my $eccapprox = "no";

# - What is the description of the run - to be used for the checkpoint files
my $description = "debug";

# - We want only to generate condor submission files
my $generatesubfileonly = "no";

# - distance from observer: D = 2*mu = 2*m1*m2/(m1+m2)
my $D = 5.0e8*$SI_pc;
#my $D = 2*$m1*$m2/($m1+$m2);

# - orbital time (defined here, set later)
#my $T = 10;
my $T = 2.*$PI*$r/($SI_c*sqrt(($m1+$m2)/$r));

# - orbit frequency
my $f = 1./$T;

# - Time step
#my $dt = 1e-3*$T;

# - logging level
my $loglevel = 0;

# - simulations stops when r<rmin or r>rmax
my $rmin = 1*($m1 + $m2);
my $rmax = 1e50*($m1 + $m2);

# - Gauge params for RR
my $eta = $m1*$m2/($m1+$m2)/($m1+$m2);
my $galpha = 4;
my $gbeta = 5;
my $gdelta1 = -99./14.+27.*$eta;		
my $gdelta2 = 5.*(1. - 4.*$eta);
my $gdelta3 = 274./7. + 67./21.*$eta;
my $gdelta4 = 5./2.*(1. - $eta);
my $gdelta5 = -1./7.*(292. + 57.*$eta);
my $gdelta6 = 51./28. + 71./14.*$eta;	

#my      $flow=10.0;
#        $r=(($SI_c/2./$PI/$flow)**2.*($m1+$m2) )**(1./3.)*(1-$epsilon);
#        $T = 2.*$PI*sqrt($r**3./($m1+$m2))/$SI_c;
#        $f = 1./$T;
#        $dt = 1./16384.;

# - loop variables
my $i;
my $j;


#
# Subroutine declarations
#

sub Run {
	my ($cmd) = @_;
        printf(STDERR "%s\n", $cmd);
        system($cmd);
}

sub WriteInifile {

open(myFILE,">test.ini");
print myFILE <<EOF;

#
# Initialization file for the CBwaves executable
#       generated by cbwgen.pl
#

[output]
outfile = $outfile              # output filename
ftfile = $ftfile                # file for the Fourier-transformed
outvars = $outvars

[input]
m1        = $m1                         # mass of bigger star
m2        = $m2                         # mass of lighter star
tmax      = $tmax                       # maximum evolution time
orbitsmax = $orbitsmax                  # maximum number of revolutions
T         = $T		                	# orbit time
f         = $f                          # orbiting frequency
dt        = $dt                         # RK4 time step
epsilon   = $epsilon                    # eccentricity
rmin      = $rmin                       # final relative distance
rmax      = $rmax                       # maximal relative distance
r         = $r                          # minimal separation
ri        = $ri                         # initial separation for hyperbolic orbits
D         = $D                          # distance to observer
iota      = $iota                       # polar angle in the source frame
phi       = $phi                        # polar angle in the source fram
theta     = $theta                      # Euler angle between TT and detector frame
varphi    = $varphi                     # Euler angle between TT and detector frame
psi       = $psi                        # Euler angle between TT and detector frame
s1x       = $s1x                        # spin x components of m1
s1y       = $s1y                        # spin y components of m1
s1z       = $s1z                        # spin z components of m1
s2x       = $s2x                        # spin x components of m2
s2y       = $s2y                        # spin y components of m2
s2z       = $s2z                        # spin z components of m2
hterms    = $hterms                     # PN order for the waveform
corrs     = $corrs                      # PN order for the motion
checkpoint= $checkpoint                 # wheter to make or not checkpoint files
eccapprox = $eccapprox                  # approximate or not the initial eccentricity
description = $description              # human readable description of the run
printstep  = 1                       # print variables after every X step
printorbit = 0                          # print variables after every X orbit
loglevel   = $loglevel                  # logging level (0 - 6)
adaptive   = no
adaptive_step = 1000
alpha  = $galpha;                       # newtonian radiation term gauge parameters
beta   = $gbeta;
delta1 = $gdelta1;                      # post-newtonian radiation term gauge parameters
delta2 = $gdelta2;
delta3 = $gdelta3;
delta4 = $gdelta4;
delta5 = $gdelta5;
delta6 = $gdelta6;
EOF
close(myFILE);
}


WriteInifile;
#Run "$cbwbinary ${filenameprefix}.ini";
print "$cbwbinary test.ini \n";


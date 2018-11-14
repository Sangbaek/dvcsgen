# dvcsgen
dvcs/pi0/eta  generator using pdfs and gpds. 
To get command line options `./dvcsgen --help`

Code requires GPD grid files. Set the `CLASDVCS_PDF` variable to point to directory where the
gpd.dat file is located (ex. CLASDVCS_PDF=/scratch/username/dvcsgen )


```
 dvcspi0gen [options]
       option  value      default    comment
      --pi0                    exclusive pi-0 on
      --eta                    exclusive eta on
      --nodvcs                 DVCS off
      --v    verbos_level    0  additional printout
      --trig nevents  10      number of triggers
      --nmax nevents   40k     # of triggers per file
      --beam beam_energy   5.754 e- momentum in GeV
      --ycol P_1 cut        0.005      P_1>ycol_cut
      --x  x-min   x-max   0.1 0.65   min max x_Bj
      --q2 Q2-min Q2-max 1.0  10. min max Q2 in GeV^2
      --y y-min y-max 0.05 0.9    min max for y=nu/E
      --w w^2-min         4.0 min for w^2=M^2+2MyE-Q^2
      --t tmin tmax  0 1.0      t  min/max in GeV^2
      --th thmin thmax  0.2 1 theta min/max \theta rad
      --zpos z-position      0 target z position in cm
      --zwidth z-width 0  width z in cm (zpos+/-zwidth)
       --raster diameter 0.75   raster diameter in cm
      --weight   flat distributions with weight(part12)
      --printgpd               print gpds and exit
      --nont               do not write out the ntuple
      --file              dvcspi0gen   filename
      --rndm               randomize
      --gpd  Igpd 3  GPD model(1-A,2-B,3-C,4-D)>=100VGG
      --scale  scale      1.0   scale the sinphi-mom
      --targ target       proton   deut/neut possible
      --lpol                    Long.pol.target
      --tpol                    Trans.pol.target
      --writef    format      0-lund12, 1-lundgsim
      --mod     write-mod      0-all, 1-cut on events
      --mom                include moments in ntuple
      --proloss                  add proton loss
      --ktcor          FALSE   turn on k_t cor for A_LU
      --radgen                   include radgen
      --nodat               do not write a data file
      --acce16            include e16 acceptance for e-
      --acceg1           include eg1 acceptance for e-
      --acc12         include clas++ acceptance for e-
      --smear                   smear moments
      --A    value         0.006   smear param-A
      --B    value         0.001  smear param-B
      --C    value         0.0008  smear param-C
      --D    value         0.001  smear param-D
      --nmax   value     2000  maximum events per file
      --print nprint     1000   print ev nprint event
      --bh  value      3 BH status:3-All, 1-only BH

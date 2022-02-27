c       The i,i+2 contacts purely repulsive

        program cg
        implicit double precision(a-h,o-z)

        parameter(len=10000) !maximum number of all residues together
        character aseq*3,pdbfile*32,outfile*32,savfile*32,state*2,buf1*4
        character rstfile*64,seqfile*32,arg*32,stafile*32,filname*32
        character buffer*128,buf2*4,paramfile*32,mapfile*32,cmapf*32
        
        logical lconftm,lmedian,lthermo,lforce,lunfold,lvelo,lmass,ldisp
        logical lchiral,lpullrel,lallatom,langle,lj1012,lconstvol,lradii
        logical loscillate,l3rdcn(len),lfrmscrtch,lwarmup,lshear,ldelrst
        logical lfrompdb(len),lwritemap,lsimpang,ldi,lsim,lstartpdb,lpbc
        logical lgln,lmrs,lparam,lchargend,ldet,lpbcx,lpbcy,lpbcz,lcintr
        logical lcoilang,lcoildih,lsawconftm,lpid,lsqpbc,ldynss,lpullfin
        logical lteql,lposcrd,lsslj,lnatend,ldisimp,lkmt,lsselj,lminfbox
        logical lunwrap,lampstrict,ljwal,lenetab,lfcc,lcleanrst,lrepcoul
        logical lnowal,lmj,lbar,lcmap,lwals,lsink,lwrtang,lwal,lcpb,lii4
        logical lobo,lwritego,lcdnat,lcospid,lepid,lrmsmax,lecperm,lsldh
        logical lrst,lcontin,lconect(len),lpdb,ldens,lmaxforce,lwritexyz
        logical lcpot
        
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/parm/f02,f12,f32,f42,f52
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/equil/kteql,ktrest,nratvel,nratveld,kconnecttime
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
        common/mass/rmas(len),rsqmas(len)
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
        common/masscenter/xcm,ycm,zcm,xmcm,ymcm,zmcm
        common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
        common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
        common/restr/bckbmin,bckb2min,sdchnmax,sfact,kmaxbckb,lecperm
        common/restart/delta,work,sep0,rstfile,stafile,filname,klenstr
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/xyforces/xuforce,xdforce,yuforce,ydforce
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
        common/mr/tene(18001,3,3,2),pp,lenetab
        dimension ttab(150) ! table of temperatures
        dimension kbt(len*50),fbt(len*50),fbt2(len*50)
        dimension kut(len*50),fut(len*50),fut2(len*50)
        dimension mtraj(len*20),ncorders(len)
        dimension tfold(501) ! folding time for each trajectory
        dimension works(501),youngmod(4,501)! work and E per oscillation
        dimension aufres(100000),aufres2(100000),aree(100000)
        dimension adfres(100000),adfres2(100000)
        dimension caac(len*50),kaak(len*50)

        pdbfile='1ubq.pdb'         ! PDB INPUT FILEAME
        seqfile='glut.txt'         ! SEQUENCE INPUT FILE
        outfile='saw.out'          ! OUTPUT FILENAME
        savfile='saw.pdb'          ! FILENAME FOR SAVING CONFORMATIONS
        rstfile='saw.rst'          ! RESTART FILE
        paramfile='parameters.txt' ! PARAMETER FILE
        stafile='(a,0,a)'          ! FORMAT STRING FOR OTHER FILES
        ! SIMULATION MODE ----------------------------------------------
        lpdb=.false.      ! USING A PDB FILE FOR GO-MODEL SIMULATIONS
        lforce=.false.    ! USING EXTERNAL FORCE FOR AFM-LIKE STRETCHING
        lvelo=.false.     ! CONST VELO. WITH AFM, MUST BE FALSE IF LWALL
        lwal=.false.      ! EXISTENCE OF THE SIMULATION BOX
        lwals=.false.     ! WALLS IN ALL 3 DIRECTIONS, NOT ONLY Z
        ldens=.false.     ! SYSTEM IS SQUEEZED TO REACH TARGET DENSITY
        lmedian=.false.   ! CALCULATE THE MEDIAN FOLDING TIME
        lthermo=.false.   ! CALCULATE THERMODYNAMIC PROPERTIES
        lunfold=.false.   ! .TRUE. FOR STUDY OF UNFOLDING
        loscillate=.false.! WALL OSCILLATION MEASUREMENTS
        lshear=.false.    ! TRUE FOR SHEARING SIMULATIONS
        ! SIMULATION OPTIONS -------------------------------------------
        lfcc=.false.      ! FCC WALLS
        lstartpdb=.false. ! START FROM PDB FILE
        lnatend=.false.   ! END WHEN ALL NATIVE CONTACTS PRESENT
        lconftm=.false.   ! CALCULATE TIME FOR CONTACT FORMATION
        lsawconftm=.true. ! USE SELF-AVOIDING WALK FOR THAT TASK
        lallatom=.false.  ! GO MODEL WITH ALL-ATOM CONTACT MAP
        lnowal=.false.    ! PBC EVEN FOR WALL DIRECTIONS
        lwarmup=.false.   ! WARM UP THE SYSTEM STARTING FROM NAT
        lteql=.false.     ! USE Teql FOR EQUILIBRATION, THEN USE T
        lrst=.false.      ! RESTART FROM SAVED STATE
        ldelrst=.true.    ! DELETE OLD RESTART FILES DURING SIMULATION
        lfrmscrtch=.false.! RESTART FROM PDB FILE (NO VELOCITIES)
        lminfbox=.false.  ! EXTEND SQUEEZED BOX TO FIND FORCE MINIMUM
        lpullfin=.true.   ! PULL THE BOX WHEN SIMULATION FINISHES
        lobo=.false.      ! ATTACH RESIDUES TO WALL ONE BY ONE
        lwritemap=.false. ! WRITE MAPFILE WITH CONTACTS
        lwritego=.false.  ! WRITE GO CONTACTS IN MAPFILE FORMAT
        lwrtang=.true.    ! IF LWRITEMAP IS TRUE WRITE ALSO ANGLES
        lwritexyz=.false. ! WRITE COORDINATES IN XYZ FILE INSTEAD OF PDB
        lposcrd=.false.   ! SHIFT COORDINATES TO STAY POSITIVE
        lcdnat=.false.    ! USE CUSTOM CUTOFF PARAMETERS FOR GO MAP
        lconstvol=.false. ! KEEP CONSTANT VOLUME DURING STRETCHING
        lcmap=.false.     ! FOR CONTACT MAP LOADED FROM A FILE
        ldisp=.false.     ! FOR DISPLACING ONE PROTEIN FROM ANOTHER
        lpbc=.false.      ! PERIODIC BOUNDARY COND-S IN EVERY DIRECTION
        lpbcx=.false.     ! PERIODIC BOUNDARY CONDITIONS IN X DIRECTION
        lpbcy=.false.     ! PERIODIC BOUNDARY CONDITIONS IN Y DIRECTION
        lpbcz=.false.     ! PERIODIC BOUNDARY CONDITIONS IN Z DIRECTION
        lcpb=.false.      ! PBC ALSO DURING GENERATING CHAINS BY SAW
        lunwrap=.false.   ! UNWRAP INPUT PDB FILE FROM PBC
        lsqpbc=.true.     ! PBC ALSO DURING SQUEEZING
        lampstrict=.true. ! AMPLITUDE OF OSCIL. IS SET TO ampstrict
        lmaxforce=.false. ! GATHER DATA ABOUT FMAX BETWEEN 2 THRESHOLDS
        lpullrel=.false.  ! PULL WITH CONSTANT VEL THEN RELEASE (mpull)
        ldet=.false.      ! FOR WRITING DETAILS IN THE CONTACT MAP
        lkmt=.true.       ! USE KMT ALGORITHM TO COMPUTE KNOT POSITIONS
        lrmsmax=.false.   ! STOP SIMULATION IF RMS IS BIGGER THAN RMSMAX
        ! SIMULATION ITSELF --------------------------------------------
        lparam=.true.     ! IF NO PARAMFILE IS SPECIFIED, LSIM=.TRUE.
        lcleanrst=.true.  ! IF CONTACTS SHOULD BE PURGED EVERY KRST TAU
        lsim=.false.      ! FOR MODEL WITH SIMPLIFIED CONTACT POTENTIAL
        ldisimp=.false.   ! HARMONIC DIHEDRAL POTENTIAL
        lmass=.false.     ! FOR TAKING INTO ACCOUNT DIFF. MASSES
        lchargend=.false. ! SHOULD PROTEIN ENDS BE CHARGED
        lsimpang=.false.  ! FOR MODEL WITH SIMPLIFIED ANGLE POT.
        lcoilang=.false.  ! RAND. COIL ANGLE POTENTIAL FOR GOMODEL
        lcoildih=.false.  ! RAND. COIL ANG. POT. FOR DIHEDRAL POTENTIAL
        lrepcoul=.false.  ! ONLY REPULSIVE PART OF COULOMB POT, NO ADIAB
        lecperm=.true.    ! USE CONSTANT REL. PERMITTIVITY FOR COULOMB
        lcpot=.true.   ! CUSTOM POTENTIAL NOT BASED ON NATIVE STRUCTURE
        lcintr=.true.  ! CUSTOM POTENTIAL ALSO FOR INTRACHAIN CONTACTS
        lpid=.false.   ! PSEUDO IMPROPER DIHEDRAL POTENTIAL
        lcospid=.false.! COSINE FUNCTION FOR THE PID POTENTIAL
        lepid=.false.  ! USE PID POTENTIAL TO CALCULATE ELECTROSTATICS
        lbar=.false.   ! SOLVATION BARRIER FOR PID POTENTIAL
        lmj=.false.    ! CONTACT DEPTH SCALED BY MIYAZAWA-J-LIKE MATRIX
        lchiral=.false.! POTENTIAL OF CHIRALITY
        langle=.true.  ! FOR MODEL WITH ANGLE POTENTIALS
        ldi=.true.     ! FOR DIHEDRAL TERM IN POTENTIALS, IF LANGLE IS T
        lii4=.true.    ! i+4 CONTACTS ARE INCLUDED BY DEFAULT, EXCEPT BB
        lsink=.false.  ! FOR SINK-LIKE POTENTIAL FOR NON-NATIVE CONTACTS
        lradii=.false. ! IF AMINO ACIDS SHOULD HAVE RADII FROM PARAMFILE
        lgln=.false.   ! FOR COUNTING ONLY GLNs IN CONTACTS
        lmrs=.false.   ! MORSE POTENTIAL FOR SSBONDS (DEFAULT HARMONIC)
        ldynss=.false. ! DYNAMIC FORMING OF DISULFIDE BRIDGES
        lsslj=.false.  ! TREAT SS BONDS AS NORMAL L-J CONTACTS
        lsselj=.false. ! TREAT SS BONDS AS EXCLUSIVE L-J CONTACTS
        lj1012=.false. ! FOR L-J 10-12 POTENTIALS
        ljwal=.true.   ! LENARD-JONES POTENTIAL FOR WALL INTERACTIONS
        lenetab=.false.! TABULARIZED VALUES FOR LOCAL POTENTIALS
        lsldh=.false.  ! SLOWLY TURN ON DIHEDRAL POTENTIAL IN THE START
        ! SIMULATION PARAMETERS - GENERAL ------------------------------
        c216=2.0d0**(1.d0/6.d0)
        Pi=dacos(-1.d0)! PI
        twopi=2.d0*Pi  ! 2PI
        unit=5.d0      ! LENGTH UNIT
        iseed=448      ! RANDOM SEED
        ntraj=1        ! NUMBER OF TRAJECTORIES
        mskip=0        ! SKIPPING STEPS [tau]
        af=4.d0        ! UNIT CELL SIZE FOR FCC WALLS [Angstrem]
        mstep=3000000  ! TOTAL SIMULATION TIME [tau]
        ktrest=10000   ! TIME AFTER SQUEEZING, BEFORE PULLING [TAU]
        kteql=0        ! EQUILIBRATION TIME BEFORE STRETCHING ETC [TAU]
        kwrite=100     ! HOW OFTEN TO PRINT ENERGY, 0: DON'T PRINT
        ksave= 1000    ! HOW OFTEN TO SAVE CONFORMATION, 0: DON'T SAVE
        kksave=1000    ! KSAVE BEFORE COOLDOWN (for kteql equilibration)
        krst=0         ! RESTART
        tstart=0.35    ! STARTING TEMPERATURE
        tend=0.35      ! ENDING TEMPERATURE
        tstep=0.05     ! TEMPERATURE STEP
        nen1=1         ! FIRST SIMULATED INDEX (INDICES<nen1 ARE FROZEN)
        rmsmax=10.0    ! STOP THE SIMULATION IF RMSD IS BIGGER THAN THIS
        ! SIMULATION PARAMETERS - DISULFIDE BRIDGES
        nssb=0         ! NUMBER OF NATIVE DISULFIDE BRIDGE, 0: NO BRIDGE
c        ksb(1,1)=56   ! FOR SSBONDS PUT RES NUMBER AS IS IN THE PROGRAM
c        ksb(2,1)=268  ! THIS ADDS 1 SSBOND BETWEEN RESIDUES 56 AND 268
        disul=1.d0     ! 5.d0 STRENGTH OF A DISULFIDE BRIDGE
        dislj=4.d0     ! HOW MANY TIMES L-J SSBOND IS STRONGER THAN L-J
        amrs=0.5       ! REVERSE WIDTH OF MORSE POTENTIAL [Angstrem^-1]
        smorse=15.0    ! SPRING CONSTANT EQUIVALENT TO MORSE [eps/A^2]
        rmrs=5.9       ! MORSE POT. MINIMUM FOR C-C PAIR [Angstrems]
        neimin=0       ! MAX. NUM. OF NEIGHBORS TO INCREASE COORD NUMBER
        neimaxdisul=9  ! MIN. NUM. OF NEIGHBORS TO PREVENT SS OXIDATION
        ! SIMULATION PARAMETERS - ANGLE POTENTIALS
        echi=1.0       ! COEFF. FOR CHIRALITY POTENTIAL [epsilon]
        CBA=30.d0!20.0COEFF FOR BOND ANGLE POTENTIALS,   20 for Clementi
        CDA=0.66 !1.0 COEFF FOR DIHED. ANGLE POT. (K1),   1 for Clementi
        CDB=0.66 !0.5 COEFF FOR DIHED. ANGLE POT. (K3), 0.5 for Clementi
        CDH=3.33 !1.0 COEFF FOR DIHED. ANGLE POT. (HARMONIC), 1 for Clem
        ! SIMULATION PARAMETERS - PULLING AND STRETCHING
        dar=1.d-10     ! DISPLACEMENT FOR COMPUTING THE FORCES
        coef=0.01      ! COEFFICIENT FOR CONST FORCE  [epsilon/(a*tau)]
        velo=0.005     ! PULLING VELOCITY IN ANGSTROMS/TAU
        veldist=12.d0  ! DISTANCE BEFORE ACQUIRING FULL VELOCITY
        mpull=100      ! TIME OF PULLING BEFORE RELEASE [tau]
        HH1=30.0       ! THE PULLING SPRING [epsilon/A^2]
        HH2=0.         ! THE PULLING SPRING [epsilon/A^2]
        cofp=1.0       ! COEFFICIENT OF THE PULLING FORCE
        dnaver = 0.0   ! LENGTH OVER WHICH FORCE IS AVERAGED (Angstrem)
        naver= 100     ! AVERAGING INSTATANEOUS F OVER MD STEPS (tau)
        kwforce=naver  ! 100     HOW OFTEN TO RECORD THE FORCE [tau]
        !       100000 is maximum simulation time divided by kwforce
        ! SIMULATION PARAMETERS - SIMULATION BOX
        tdens=0.001    ! TARGET DENSITY (IN RESIDUES/ANSTREM^3)
        sdens=0.0001   ! INITIAL DENSITY (THE BOX WILL SQUEEZE TO TDENS)
        densvelo=0.02  ! VELOCITY OF SQUEEZING FROM EVERY DIRECTION
        sepmin=10.0    ! FINAL DISTANCE BETWEEN WALLS [Angstrems]
        kbperiodmax=6  ! NUMBER OF OSCILLATIONS AFTER WHICH IS PULLING
        amp=0.1        ! AMPLITUDE OF OSCILLATIONS, RELATIVE TO WALLDIST
        ampstrict=10.0 ! AMPLITUDE OF OSCILLATIONS, STRICT [Angstrem]
        omega=0.0001   ! ANGULAR FREQUENCY OF OSCILLATIONS [1/tau]
c       period around 18000 tau means vmax=0.01 A/tau for large systems
        period=twopi/omega   ! TIME PERIOD OF OSCILLATIONS [tau]
        walmindst=5.d0 ! MINIMAL DISTANCE FOR STICKING TO WALL [A]
        kconnecttime=9 ! WHEN TO CONNECT BEADS TO WALL (1,3,5 OR 7)
        ipwn=-1        ! NR OF RES. STICKED TO WALL (NEGATIVE-AUTO)
        fwal = 4.0     ! WALL POTENTIAL COEFFICIENT
        ! SIMULATION PARAMETERS - TECHNICAL
        alphacos(1)=6.4! INVERSE WIDTH OF THE 1ST BB PID POTENTIAL 
        alphacos(2)=6.0!3.6! INVERSE WIDTH OF THE 2ND BB PID POTENTIAL
        alphacos(3)=1.2! INVERSE WIDTH OF THE SS PID POTENTIAL
        rbb(1)=5.6     ! R_MIN OF THE 1ST BB PID POTENTIAL 
        rbb(2)=6.2     ! R_MIN OF THE 2ND BB PID POTENTIAL 
        psi0ss=-0.23   ! ANGLE MINIMUM OF THE 1ST BB PID POTENTIAL
        psi0bb(1)=1.05 ! ANGLE MINIMUM OF THE 2ND BB PID POTENTIAL
        psi0bb(2)=-1.44! ANGLE MINIMUM OF THE SS PID POTENTIAL
        epsbb=0.2      ! AMPLITUDE OF BB INTERACTIONS in PID [epsilon]
        rnei=7.5       ! CUT-OFF DISTANCE FOR NEIGHBOUR RES. [Angstrem]
        ad=2000        ! TIME FOR ADIABATIC POTENTIAL TURNOFF [steps]
        potcoeff=1.d0  ! SCALE OF LCPOT FORCES
        cntfct=1.3d0   ! MULTIPLIES SIGMA TO CHECK IF CONTACT PRESENT
c        factor=2.0    ! HOW MANY TIMES BACKB-BACKB. TOLERANCE IS BIGGER
        confcut=4.56   ! MIN. DISTANCE FOR CHAINS GENERATED BY SAW [A]
        bond=3.8       ! BOND LENGTH [Angstrem]
        ethi=1.0       ! THIRUMALAI OVERLAP FUNCTION [Angstrem] unused
        dnat=0.0       ! CUT-OFF DISTANCE FOR NATIVE CONTACTS [Angstrem]
        delta=0.005    ! INTEGRATION TIME STEP [tau]
        gamma=2.0      ! LANGEVIN PARAMETER [m/tau]
        H1 = 50.0      ! HARMONIC COEFF. [epsilon/A^2]
        H2 = 0.0       ! ANHARMONIC COEFF. [epsilon/A^4]
        verlcut=10.0   ! CUT-OFF DISTANCE FOR VERLET LIST[Angstrem]
        tolerance=0.d0 ! HOW FAR FROM THE MINIMUM CONTACT TURNS ON [%]
        sigma1(1)=1.5  ! FOR LEGACY (UNLESS LSIM IS T SHOULD NOT MATTER)
        sigma1(2)=1.5  ! FOR LEGACY (UNLESS LSIM IS T SHOULD NOT MATTER)
        sigma1(3)=1.5  ! FOR LEGACY (UNLESS LSIM IS T SHOULD NOT MATTER)
        sigma1(4)=5.0  ! L-J MINIMUM FOR BACKBONE CONTACTS [Angstrems]
        sigma1(5)=7.5  ! L-J MINIMUM FOR SIDECHAIN CONTACTS [Angstrems]
        sigma1(6)=6.6  ! L-J MINIMUM FOR BCKB-SDCH CONTACTS [Angstrems]
        sigma1(7)=6.6  ! L-J MINIMUM FOR SDCH-BCKB CONTACTS [Angstrems]
        sigma1(8)=6.1  ! L-J MINIMUM FOR i,i+4 CONTACTS [Angstrems]
        screend=10.0   ! ELECTROSTATIC SCREENING LENGTH [Angstrems]
        coul=85.0      ! CONSTANT FOR COULOMBIC INTERACTION [eps*A*A]
        if(lecperm) coul=2.63 ! 210 IF REL. PERMITTIVITY=1 [eps*A]
        cut=5.0        ! CUT-OFF DISTANCE FOR REPULSIVE TERM [Angstrem]
        rcut=18.d0     ! CUT-OFF DISTANCE FOR THE REST [Angstrem]
        kmaxbckb=2     ! MAXIMAL NUMBER OF BACKBONE CONTACTS
        bckbmin=0.75   ! MINIMUM VALUE OF BACKBONE ANGLE
        bckb2min=0.92  ! MINIMUM VALUE OF BACKBONE ANGLE
        sdchnmax=0.5   ! MAXMUM VALUE OF SIDECHAIN ANGLE
        ndomain=1      ! NUMBER OF DOMAIN (FOR TITIN ONLY)
        ! READING VARIABLES --------------------------------------------
        if(iargc().gt.0) then ! reading parameters from external file
        call getarg(1,arg) ! for more than 1 args, do i_arg=1,iargc()
        open(7,file=arg,status='old') ! no spaces in filenames allowed
        write(*,*) 'RUNNING INPUTFILE ',arg
11      read(7,'(a)',end=12, err=12) buffer
        if(buffer(1:4).eq.'velo') then
            read(buffer(5:),*) velo
        elseif(buffer(1:7).eq.'pdbfile') then
            read(buffer(8:),*) pdbfile 
        elseif(buffer(1:7).eq.'seqfile') then
            read(buffer(8:),*) seqfile 
        elseif(buffer(1:7).eq.'outfile') then
            read(buffer(8:),*) outfile 
        elseif(buffer(1:7).eq.'savfile') then
            read(buffer(8:),*) savfile
        elseif(buffer(1:7).eq.'rstfile') then
            read(buffer(8:),*) rstfile
        elseif(buffer(1:7).eq.'mapfile') then
            read(buffer(8:),*) mapfile
            lwritemap=.true.
        elseif(buffer(1:7).eq.'verlcut') then
            read(buffer(8:),*) verlcut
        elseif(buffer(1:6).eq.'dnaver') then
            read(buffer(7:),*) dnaver
        elseif(buffer(1:6).eq.'factor') then
            read(buffer(7:),*) factor
        elseif(buffer(1:9).eq.'lconstvol') then
            read(buffer(10:),*) lconstvol
        elseif(buffer(1:9).eq.'lstartpdb') then
            read(buffer(10:),*) lstartpdb
        elseif(buffer(1:8).eq.'lslowdih') then
            read(buffer(9:),*) lsldh
        elseif(buffer(1:5).eq.'lsldh') then
            read(buffer(6:),*) lsldh
        elseif(buffer(1:5).eq.'lcmap') then
            read(buffer(6:),*) lcmap
        elseif(buffer(1:5).eq.'lsink') then
            read(buffer(6:),*) lsink
        elseif(buffer(1:4).eq.'lpid') then
            read(buffer(5:),*) lpid
        elseif(buffer(1:4).eq.'lfcc') then
            read(buffer(5:),*) lfcc
        elseif(buffer(1:4).eq.'lii4') then
            read(buffer(5:),*) lii4
        elseif(buffer(1:4).eq.'lbar') then
            read(buffer(5:),*) lbar
        elseif(buffer(1:4).eq.'lrst') then
            read(buffer(5:),*) lrst
        elseif(buffer(1:6).eq.'lradii') then
            read(buffer(7:),*) lradii
        elseif(buffer(1:6).eq.'lcdnat') then
            read(buffer(7:),*) lcdnat
        elseif(buffer(1:5).eq.'ljwal') then
            read(buffer(6:),*) ljwal
        elseif(buffer(1:5).eq.'lepid') then
            read(buffer(6:),*) lepid
        elseif(buffer(1:5).eq.'disul') then
            read(buffer(6:),*) disul
        elseif(buffer(1:5).eq.'dislj') then
            read(buffer(6:),*) dislj
        elseif(buffer(1:5).eq.'epsbb') then
            read(buffer(6:),*) epsbb
        elseif(buffer(1:5).eq.'gamma') then
            read(buffer(6:),*) gamma
        elseif(buffer(1:5).eq.'cmapf') then
            read(buffer(6:),*) cmapf
            lcmap=.true.
        elseif(buffer(1:9).eq.'tolerance') then
            read(buffer(10:),*) tolerance
        elseif(buffer(1:9).eq.'lcleanrst') then
            read(buffer(10:),*) lcleanrst
        elseif(buffer(1:9).eq.'lwritemap') then
            read(buffer(10:),*) lwritemap
        elseif(buffer(1:9).eq.'lwritexyz') then
            read(buffer(10:),*) lwritexyz
        elseif(buffer(1:7).eq.'lwrtang') then
            read(buffer(8:),*) lwrtang
        elseif(buffer(1:7).eq.'lcospid') then
            read(buffer(8:),*) lcospid
        elseif(buffer(1:7).eq.'lecperm') then
            read(buffer(8:),*) lecperm
            coul=210.0
        elseif(buffer(1:7).eq.'lrmsmax') then
            read(buffer(8:),*) lrmsmax
        elseif(buffer(1:11).eq.'wallmindist') then
            read(buffer(12:),*) walmindst
        elseif(buffer(1:12).eq.'lfromscratch') then
            read(buffer(13:),*) lfrmscrtch
        elseif(buffer(1:11).eq.'lwritegomap') then
            read(buffer(12:),*) lwritego
        elseif(buffer(1:8).eq.'lwritego') then
            read(buffer(9:),*) lwritego
        elseif(buffer(1:7).eq.'lunwrap') then
            read(buffer(8:),*) lunwrap
        elseif(buffer(1:10).eq.'lfrmscrtch') then
            read(buffer(11:),*) lfrmscrtch
        elseif(buffer(1:10).eq.'lampstrict') then
            read(buffer(11:),*) lampstrict
        elseif(buffer(1:12).eq.'kconnecttime') then
            read(buffer(13:),*) kconnecttime
        elseif(buffer(1:8).eq.'displace') then
            read(buffer(9:),*) iprota,jprotb,away
            ldisp=.true.
        elseif(buffer(1:8).eq.'lpullfin') then
            read(buffer(9:),*) lpullfin
        elseif(buffer(1:8).eq.'lrepcoul') then
            read(buffer(9:),*) lrepcoul
        elseif(buffer(1:8).eq.'kmaxbckb') then
            read(buffer(9:),*) kmaxbckb
        elseif(buffer(1:7).eq.'screend') then
            read(buffer(8:),*) screend
        elseif(buffer(1:7).eq.'bckbmin') then
            read(buffer(8:),*) bckbmin
        elseif(buffer(1:7).eq.'bckbmax') then
            read(buffer(8:),*) bckbmax
        elseif(buffer(1:8).eq.'sdchnmin') then
            read(buffer(9:),*) sdchnmin
        elseif(buffer(1:8).eq.'sdchnmax') then
            read(buffer(9:),*) sdchnmax
        elseif(buffer(1:8).eq.'bckb2min') then
            read(buffer(9:),*) bckb2min
        elseif(buffer(1:8).eq.'bckb2max') then
            read(buffer(9:),*) bckb2max
        elseif(buffer(1:8).eq.'lrestart') then
            read(buffer(9:),*) lrst
        elseif(buffer(1:8).eq.'lsimpang') then
            read(buffer(9:),*) lsimpang
        elseif(buffer(1:8).eq.'lcoilang') then
            read(buffer(9:),*) lcoilang
        elseif(buffer(1:8).eq.'lcoildih') then
            read(buffer(9:),*) lcoildih
        elseif(buffer(1:8).eq.'lallatom') then
            read(buffer(9:),*) lallatom
        elseif(buffer(1:9).eq.'lchargend') then
            read(buffer(10:),*) lchargend
        elseif(buffer(1:9).eq.'ampstrict') then
            read(buffer(10:),*) ampstrict
        elseif(buffer(1:8).eq.'potcoeff') then
            read(buffer(9:),*) potcoeff
        elseif(buffer(1:8).eq.'densvelo') then
            read(buffer(9:),*) densvelo
        elseif(buffer(1:6).eq.'cntfct') then
            read(buffer(7:),*) cntfct
        elseif(buffer(1:6).eq.'lsqpbc') then
            read(buffer(7:),*) lsqpbc
        elseif(buffer(1:5).eq.'lpbcx') then
            read(buffer(6:),*) lpbcx
        elseif(buffer(1:5).eq.'lpbcy') then
            read(buffer(6:),*) lpbcy
        elseif(buffer(1:5).eq.'lpbcz') then
            read(buffer(6:),*) lpbcz
        elseif(buffer(1:4).eq.'lobo') then
            read(buffer(5:),*) lobo
        elseif(buffer(1:4).eq.'lsim') then
            read(buffer(5:),*) lsim
        elseif(buffer(1:4).eq.'ldet') then
            read(buffer(5:),*) ldet
        elseif(buffer(1:4).eq.'lpbc') then
            read(buffer(5:),*) lpbc
        elseif(buffer(1:4).eq.'lcpb') then
            read(buffer(5:),*) lcpb
        elseif(buffer(1:4).eq.'lkmt') then
            read(buffer(5:),*) lkmt
        elseif(buffer(1:4).eq.'ldih') then
            read(buffer(5:),*) ldi
        elseif(buffer(1:7).eq.'ldisimp') then
            read(buffer(8:),*) ldisimp
        elseif(buffer(1:7).eq.'lmedian') then
            read(buffer(8:),*) lmedian
        elseif(buffer(1:7).eq.'lenetab') then
            read(buffer(8:),*) lenetab
        elseif(buffer(1:3).eq.'CBA') then
            read(buffer(4:),*) CBA
        elseif(buffer(1:3).eq.'CDA') then
            read(buffer(4:),*) CDA
        elseif(buffer(1:3).eq.'CDB') then
            read(buffer(4:),*) CDB
        elseif(buffer(1:3).eq.'lmj') then
            read(buffer(4:),*) lmj
        elseif(buffer(1:4).eq.'coul') then
            read(buffer(5:),*) coul
        elseif(buffer(1:4).eq.'nen1') then
            read(buffer(5:),*) nen1
        elseif(buffer(1:4).eq.'fwal') then
            read(buffer(5:),*) fwal
        elseif(buffer(1:5).eq.'adiab') then
            read(buffer(6:),*) ad
        elseif(buffer(1:5).eq.'lvelo') then
            read(buffer(6:),*) lvelo 
        elseif(buffer(1:5).eq.'lmass') then
            read(buffer(6:),*) lmass
        elseif(buffer(1:6).eq.'langle') then
            read(buffer(7:),*) langle        
        elseif(buffer(1:7).eq.'lchiral') then
            read(buffer(8:),*) lchiral
        elseif(buffer(1:7).eq.'lconftm') then
            read(buffer(8:),*) lconftm
        elseif(buffer(1:7).eq.'confcut') then
            read(buffer(8:),*) confcut
        elseif(buffer(1:7).eq.'lnatend') then
            read(buffer(8:),*) lnatend
        elseif(buffer(1:7).eq.'lposcrd') then
            read(buffer(8:),*) lposcrd
        elseif(buffer(1:7).eq.'lthermo') then
            read(buffer(8:),*) lthermo
        elseif(buffer(1:7).eq.'ldelrst') then
            read(buffer(8:),*) ldelrst
        elseif(buffer(1:6).eq.'lnowal') then
            read(buffer(7:),*) lnowal
        elseif(buffer(1:10).eq.'lsawconftm') then
            read(buffer(11:),*) lsawconftm
        elseif(buffer(1:4).eq.'lmrs') then
            read(buffer(5:),*) lmrs
        elseif(buffer(1:6).eq.'lmorse') then
            read(buffer(7:),*) lmrs
        elseif(buffer(1:6).eq.'lparam') then
            read(buffer(7:),*) lparam
        elseif(buffer(1:6).eq.'ldynss') then
            read(buffer(7:),*) ldynss
        elseif(buffer(1:6).eq.'lsselj') then
            read(buffer(7:),*) lsselj
        elseif(buffer(1:6).eq.'lcintr') then
            read(buffer(7:),*) lcintr
        elseif(buffer(1:5).eq.'lsslj') then
            read(buffer(6:),*) lsslj
        elseif(buffer(1:5).eq.'ldens') then
            read(buffer(6:),*) ldens
        elseif(buffer(1:5).eq.'tdens') then
            read(buffer(6:),*) tdens
        elseif(buffer(1:5).eq.'sdens') then
            read(buffer(6:),*) sdens
        elseif(buffer(1:4).eq.'krst') then
            read(buffer(5:),*) krst
        elseif(buffer(1:5).eq.'kteql') then
            read(buffer(6:),*) kteql
        elseif(buffer(1:5).eq.'ksave') then
            read(buffer(6:),*) ksave
        elseif(buffer(1:6).eq.'kksave') then
            read(buffer(7:),*) kksave
        elseif(buffer(1:6).eq.'ktrest') then
            read(buffer(7:),*) ktrest
        elseif(buffer(1:6).eq.'kwrite') then
            read(buffer(7:),*) kwrite
        elseif(buffer(1:6).eq.'sepmin') then
            read(buffer(7:),*) sepmin
        elseif(buffer(1:6).eq.'lwall ') then
            read(buffer(7:),*) lwal 
        elseif(buffer(1:6).eq.'lwalls') then
            read(buffer(7:),*) lwals
        elseif(buffer(1:6).eq.'lshear') then
            read(buffer(7:),*) lshear
        elseif(buffer(1:9).eq.'paramfile') then
            read(buffer(10:),*) paramfile
            lparam=.true.
        elseif(buffer(1:10).eq.'loscillate') then
            read(buffer(11:),*) loscillate
        elseif(buffer(1:11).eq.'kbperiodmax') then
            read(buffer(12:),*) kbperiodmax
        elseif(buffer(1:5).eq.'omega') then
            read(buffer(6:),*) omega
            period=twopi/omega
        elseif(buffer(1:6).eq.'period') then
            read(buffer(7:),*) period
            omega=twopi/period
        elseif(buffer(1:5).eq.'iseed') then
            read(buffer(6:),*) iseed 
        elseif(buffer(1:5).eq.'ntraj') then
            read(buffer(6:),*) ntraj
        elseif(buffer(1:5).eq.'mstep') then
            read(buffer(6:),*) mstep
        elseif(buffer(1:5).eq.'lcpot') then
            read(buffer(6:),*) lcpot
        elseif(buffer(1:4).eq.'lpdb') then
            read(buffer(5:),*) lpdb
        elseif(buffer(1:4).eq.'ipwn') then
            read(buffer(5:),*) ipwn
        elseif(buffer(1:4).eq.'bond') then
            read(buffer(5:),*) bond
        elseif(buffer(1:4).eq.'dnat') then
            read(buffer(5:),*) dnat
        elseif(buffer(1:4).eq.'rcut') then
            read(buffer(5:),*) rcut
        elseif(buffer(1:4).eq.'cofp') then
            read(buffer(5:),*) cofp
        elseif(buffer(1:4).eq.'rnei') then
            read(buffer(5:),*) rnei
        elseif(buffer(1:6).eq.'neimin') then
            read(buffer(7:),*) neimin
        elseif(buffer(1:11).eq.'neimaxdisul') then
            read(buffer(7:),*) neimaxdisul
        elseif(buffer(1:3).eq.'HH1') then
            read(buffer(4:),*) HH1
        elseif(buffer(1:3).eq.'amp') then
            read(buffer(4:),*) amp
        elseif(buffer(1:3).eq.'cut') then
            read(buffer(4:),*) cut
        elseif(buffer(1:5).eq.'acos1') then
            read(buffer(6:),*) alphacos(1)
        elseif(buffer(1:5).eq.'acos2') then
            read(buffer(6:),*) alphacos(2)
        elseif(buffer(1:5).eq.'acos3') then
            read(buffer(6:),*) alphacos(3)
        elseif(buffer(1:6).eq.'psi0ss') then
            read(buffer(7:),*) psi0ss
        elseif(buffer(1:7).eq.'psi0bb1') then
            read(buffer(8:),*) psi0bb(1)
        elseif(buffer(1:7).eq.'psi0bb2') then
            read(buffer(8:),*) psi0bb(2)
        elseif(buffer(1:4).eq.'rbb1') then
            read(buffer(5:),*) rbb(1)
        elseif(buffer(1:4).eq.'rbb2') then
            read(buffer(5:),*) rbb(2)
        elseif(buffer(1:4).eq.'bbrm') then
            read(buffer(5:),*) sigma1(4)
        elseif(buffer(1:4).eq.'ssrm') then
            read(buffer(5:),*) sigma1(5)
        elseif(buffer(1:4).eq.'bsrm') then
            read(buffer(5:),*) sigma1(6)
            sigma1(7)=sigma1(6)
        elseif(buffer(1:4).eq.'i4rm') then
            read(buffer(5:),*) sigma1(8)
        elseif(buffer(1:4).eq.'temp') then
            read(buffer(5:),*) tstart
            tend=tstart
        elseif(buffer(1:4).eq.'teql') then
            read(buffer(5:),*) teql
            lteql=.true.
        elseif(buffer(1:4).eq.'tend') then
            read(buffer(5:),*) tend
        elseif(buffer(1:5).eq.'tstep') then
            read(buffer(6:),*) tstep
        elseif(buffer(1:6).eq.'tstart') then
            read(buffer(7:),*) tstart
        elseif(buffer(1:7).eq.'klenstr') then
            read(buffer(8:),*) klenstr
        elseif(buffer(1:4).eq.'file') then
            read(buffer(5:),*) filname
            write(stafile,*) '(a',klenstr,',a)'
            write(outfile,stafile) filname,'.out'
            write(mapfile,stafile) filname,'.map'
            write(savfile,stafile) filname,'.pdb'
        else ! writing to console, unless file indexes 5 or 6 are in use
            write(*,*) 'UNRECOGNIZED OPTION: ',buffer
        endif
        goto 11
12      close(7)
        endif
        if(lparam) call load_paramfile(paramfile)
        do i=1,len      ! LFROMPDB MUST BE ZEROED BEFORE LOADING CMAPS
            the0(i)=-1.0
            phi0(i)=0.0
            lfrompdb(i)=lsimpang
        enddo
        if(lpbc) then
            lpbcx=.true.
            lpbcy=.true.
            lpbcz=.true.
        endif
        klont=0
        if(lwal.and..not.lnowal) lpbcz=.false. ! no PBC in Z
        if(lsselj) lmrs=.true. ! lmorse is overriden by lsselj
        if(dnaver.gt.0.0) naver=nint(dnaver/velo)
        nratvel=nint(veldist/velo) ! NUMBER OF STEPS TO PULLING VELOCITY
        nratveld=nint(veldist/densvelo) ! N OF STEPS TO SQUEEZE VELOCITY
        open(1,file=outfile,status='unknown')
        if(lwritemap) open(22,file=mapfile,status='unknown')
        if(ksave.ne.0) open(2,file=savfile,status='unknown')

        write(1,*)'#I,I+2 CONTACTS PURELY REPULSIVE'
        
        ! SCALE LENGTHS
        cut=cut/unit  ! REPULSIVE INTERACTIONS CUT-OFF
        sfact=(1.d0+tolerance)*c216 !2.0d0**(1.d0/6.d0)
        sigma0=cut/c216 !*0.5d0**(1.d/6.d)! REPULSIVE INTERACTIONS SIGMA
c       dfold=dfold/unit
        rcut = rcut/unit ! OTHER POTENTIALS CUT-OFF
        rnei=rnei/unit
        rcutsq=rcut*rcut
        cutsq=cut*cut
        rneisq=rnei*rnei
c        bckbmin=bckbmin*bckbmin
c        bckb2min=bckb2min*bckb2min
        sdchnmax=sdchnmax*sdchnmax
        verlcut=verlcut/unit
        confcut=confcut/unit
        vrcut2sq=verlcut*verlcut/4.d0
        dnat=dnat/unit
        H1=H1*unit*unit
        H2=H2*unit**4
        HH1=HH1*unit*unit
        HH2=HH2*unit**4
        ethi=ethi/unit
        bond=bond/unit
        af=af/unit
        rbb(1)=rbb(1)/unit*0.5d0**(1.d0/6.d0)
        rbb(2)=rbb(2)/unit*0.5d0**(1.d0/6.d0)
        if(.not.lcospid) then
            alphacos(1)=alphacos(1)/Pi
            alphacos(2)=alphacos(2)/Pi
            alphacos(3)=alphacos(3)/Pi
            vmp=1.d0
        else
            vmp=Pi
        endif
        ampstrict=ampstrict/unit
        sepmin=sepmin/unit
        sigma1(9)=walmindst ! for wall interactions via L-J potential
        walmindst=walmindst/unit
        tdens=tdens*unit**3
        sdens=sdens*unit**3
        vpull=velo/unit
        vpull=vpull*delta ! corresponding infinitesimal displacement
        vpulld=densvelo/unit*delta
        vtarget=vpull ! target velocity acquired after nratvel steps
        ! SETUP TABLE OF TEMPERATURE
        nt=nint(abs(tstart-tend)/tstep)+1
        if(nt.gt.150) nt=150
        ttstep=tstep
        if(tstart.gt.tend) ttstep=-tstep
        do i=1,nt
        ttab(i)=tstart+(i-1)*ttstep
        enddo
c        bckbmax=bckbmax-bckbmin
c        bckb2max=bckb2max-bckb2min
c        sdchnmin=sdchnmax-sdchnmin
        
        ! LOAD PROTEIN CONFORMATION
        if(lpdb.or.lstartpdb) then
            write(1,'(/,a,2x,a,/)')'#PDB FILE =',pdbfile
            call load_protein(pdbfile,lunwrap)
            do i=1,men
            x0(i)=xn(i)
            y0(i)=yn(i)
            z0(i)=zn(i)
            enddo
        else
            write(1,'(/,a,2x,a,/)')'#SEQ FILE =',seqfile
            call load_sequence(seqfile)
            call confstart(sdens,confcut)
            do i=1,men
            xn(i)=x0(i)
            yn(i)=y0(i)
            zn(i)=z0(i)
            enddo
        endif
        
        if(lallatom) then
            write(1,'(a)')'#CONSTRUCT CONTACT-MAP BASED ON ALL-ATOM'
            call compute_contact_map(pdbfile)
        elseif(klont.eq.0) then
        write(1,'(a,f6.2,a)')'#CONSTRUCT CONTACTMAP BASED ON CUT-OFF =',
     +      dnat*unit,' ANGSTROM'
            call compute_cmap(dnat,lcdnat)
        else
            write(1,'(a)')'#CONTACT-MAP BASED ON THE SEQUENCE FILE'
        endif
        
        if(lcmap) then
          write(1,'(a,a)')'#CONSTRUCT CONTACT-MAP BASED ON FILE ',cmapf
          call load_cmap(cmapf)
        endif
        if(ipwn.lt.0) then 
        ipwn=(men-nint((men**(1./3)-walmindst*tdens**(1./3))**3))/3
        endif
        ipwn=2*(ipwn/2)
        ! BUILD TITIN
        if(ndomain.gt.1) then
            write(1,'(a,i10)')'#NUMBER OF DOMAINS',ndomain
            call build_titin(ndomain)
            write(1,'(a)')'#DO NOT ALLOW CONTACTS BETWEEN DOMAINS'
            call interdomain(ndomain)
        endif
        
        ip1=1
        ip2=men
        if(lmass) then
            write(1,'(a)')'#CONSIDERING AMINO ACID MASSES'
            call amino_acid_mass
        else
            do i=1,men
                rmas(i)=1.d0
            enddo
        endif
        jq=1
        part=men-nen1+1
        mcmr=men*(men-1)/2-(men-1)
        mchi=(men*men-5*men+6)/2
        targetvolume=men/tdens
        reemax=3.0*targetvolume**(1./3)
        xmin=x0(1)
        ymin=y0(1)
        zmin=z0(1)
        xmax=x0(1)
        ymax=y0(1)
        zmax=z0(1)
        if(.not.lcpb) then
            do ib=2,men
                if(x0(ib).lt.xmin) xmin=x0(ib)
                if(y0(ib).lt.ymin) ymin=y0(ib)
                if(z0(ib).lt.zmin) zmin=z0(ib)
                if(x0(ib).gt.xmax) xmax=x0(ib)
                if(y0(ib).gt.ymax) ymax=y0(ib)
                if(z0(ib).gt.zmax) zmax=z0(ib)
            enddo
            xdown=xmin-2*bond
            xup=xmax+2*bond
            ydown=ymin-2*bond
            yup=ymax+2*bond
            zdown=zmin-2*bond
            zup=zmax+2*bond
            xsep=xup-xdown
            ysep=yup-ydown
            zsep=zup-zdown
            xinv=1.d0/xsep
            yinv=1.d0/ysep
            zinv=1.d0/zsep
        endif
        if(lobo) kconnecttime=3
        call update_verlet_list(verlcut,nen1)
        corder=0
        icor=0
        do k=1,kront
C            if(krist(3,k).ne.0) then
            corder=corder + abs(krist(1,k)-krist(2,k))
            icor=icor+1
C            endif
        enddo
        if(icor.gt.0) then
            corder=corder/icor
            corder=corder/men
            write(1,'(a,f8.4)') '#RELATIVE TOTAL CONTACT ORDER', corder
        else
            write(1,'(a)') '#NO NON-REPULSIVE CONTACTS'
        endif
        if(klont.gt.0) then
            flush(1)
            corder=0
            do k=1,klont
            corder=corder + abs(klist(1,k)-klist(2,k))
            enddo
            corder=corder/klont
            corder=corder/men
            write(1,'(a,f8.4)') '#RELATIVE NATIVE CONTACT ORDER', corder
        else
            write(1,'(a)') '#NO NATIVE CONTACTS'
        endif
        
        do 47 ks=1,nssb
        i1=ksb(1,ks)
        i2=ksb(2,ks)
        if(aseq(i1).ne.'CYS') then
        write(1,*)'#RESIDUE ',ksb(ks,1),' IS NOT A CYSTEINE. PLS CHECK!' 
        stop
        endif
        if(aseq(i2).ne.'CYS') then
        write(1,*)'#RESIDUE ',ksb(ks,2),' IS NOT A CYSTEINE. PLS CHECK!' 
        stop
        endif
        icheck=0
        if(.not.(lsslj.or.lsselj)) then
            do k=1,klont
                ki1=klist(1,k)
                ki2=klist(2,k)
        if((ki1.eq.i1.and.ki2.eq.i2).or.(ki1.eq.i2.and.ki2.eq.i1)) then
        klist(3,k)=sign(631,klist(3,k)) ! SS BONDS HAVE VALUE +-631
                icheck=1
                write(1,'(a,i4,5x,2(a3,i4,3x))')'#NATIVE SS BOND',ks,
     +          aseq(i1),iseq(i1),aseq(i2),iseq(i2)
        endif
           enddo
        if(icheck.eq.0) then
        write(1,'(a,i4,5x,2(a3,i4,3x))')'#SS BOND COULD NOT BE MADE',ks,
     +   aseq(i1),iseq(i1),aseq(i2),iseq(i2)
        endif
        endif
47      continue
c        write(1,'(a,f7.1)')'SS BOND STRENGTH',disul
        if(lcpot) then
        write(1,'(a)')'#USING CUSTOM ATTRACTIVE L-J CONTACT POTENTIAL!'
        if(lsselj) then
        write(1,'(a)')'#USING L-J POTENTIAL FOR NON-NATIVE SS BONDS!'
        write(1,'(a,f6.2,a,f6.2)')'#R_LJ=',rmrs,' DEPTH=',dislj
        else
        write(1,'(a)')'#USING MORSE POTENTIAL FOR NON-NATIVE SS BONDS!'
        write(1,'(a,f6.2,a,f6.2)')'#R_MORSE=',rmrs,' A_MORSE=',amrs
        endif
        if(lsim.or.(lcpot.and..not.lparam)) then
            lsim=.true.
        write(1,'(a,f8.2)')'#SS CONTACT EQUILIBRIUM DISTANCE=',sigma1(5)
        else
          write(1,'(a,a)')'#CONTACTS BASED ON DATA FROM FILE ',paramfile
        endif
        write(1,'(a,f6.2)') '#USING DEBYE SCREENING LENGTH',screend
        endif
        
        ! PUT CUSTOM POTENTIALS HERE
        
        smorse=smorse*unit*unit
        amrs=amrs*unit
        screend=screend/unit
        coul=coul/unit
        if(.not.lecperm) coul=coul/unit
        dmrs=smorse/amrs**2
        rmrs=rmrs/unit
        sigma1(88)=rmrs*0.5d0**(1.d0/6.d0)
        sigss=sigma1(88)
        do i=4,9
            sigma1(i)=sigma1(i)/unit*0.5d0**(1.d0/6.d0)
        enddo
c        sigma1(2)=dalQQ*0.5d0**(1.d0/6.d0)
        
        write(1,'(a)')'#USING HARMONIC POTENTIALS FOR NATIVE SS BONDS!'
        olno=ran2(iseed) ! FOR LEGACY (COMPARE OLD VERSION, SAME SEED)
        write(1,'(a,2(a,f6.2))')'#USING ANHARMONIC POTENTIAL',
     +  '   H1 =',H1/unit/unit,'   H2 =',H2/unit**4

        if(lchiral) then
        write(1,'(a)')'#USING CHIRALITY POTENTIALS'
        call model_chirality
        endif

        if(langle) then
        if(lsimpang) ldi=.false.
        if(ldi) then
          write(1,'(a)')'#USING POTENTIALS FOR BOND AND DIHEDRAL ANGLES'
        else
          write(1,'(a)')'#USING POTENTIALS ONLY FOR BOND ANGLES'
        endif
        endif
        
        if(lpdb) call compute_native_angles
        
        write(1,'(a)')'#DISABLE NATIVE CONTACTS (I,I+2)' ! AND (I,I+3)'
        km=0
        do k=1,klont
        i=klist(1,k)
        j=klist(2,k)
        ijdiff=j-i
        if(ldynss .and. inameseq(i).eq.4 .and. inameseq(j).eq.4 
     +  .and. abs(klist(3,k)).eq.631) then
        ijdiff=0
        write(1,'(a,2i4,5x,2(a3,i4,3x))')'#DYNAMIC SS BOND',i,j,
     +   aseq(i),iseq(i),aseq(j),iseq(j)
        endif
        if(.not.lconect(i)) ijdiff=5
        if(ijdiff.ge.3) then ! 4) then
        km=km+1
        klist(1,km)=i
        klist(2,km)=j
        klist(3,km)=klist(3,k)
        endif
        enddo
        klont=km
        
        ngln=0
        do ib=1,men
            x0(ib)=xn(ib)     ! assign native coordinates to actual ones
            y0(ib)=yn(ib)
            z0(ib)=zn(ib)
            adia(ib)=0       ! zero the table for adiabatic turning off
            z0temp(ib)=zn(ib) ! initialize the table to be sorted by z
            if(ljwal) z0temp(ib) = 0
            ksorted(ib)=ib    ! indexes to be sorted
            xpul(ib)=0.0      ! initialize ref. values of pulled resid.
            ypul(ib)=0.0
            zpul(ib)=0.0
            nei(1,ib)=0       ! set neighbour counter to zero
            nei(2,ib)=0
            if(aseq(ib).eq.'GLN') ngln=ngln+1 ! count glutamines
            ksdchns(ib)=ksdchn(inameseq(ib),1) ! type of sidechain
            !khbful(3,ib)=ksdchn(inameseq(ib),2) !nr of potential hbonds
        enddo
        if(lchargend) then
          do ic=1,nchains
           if(ksdchns(menchain(ic)+1).eq.4) then
            ksdchns(menchain(ic)+1)=2 !N-terminal is zwitterion
           else
            ksdchns(menchain(ic)+1)=5 !N-terminal of protein is positive
           endif
           if(ksdchns(menchain(ic+1)).eq.5) then
            ksdchns(menchain(ic+1))=2 ! C-terminal is zwitterion
           else
            ksdchns(menchain(ic+1))=4 ! C-terminal is negative
           endif
          enddo
        endif
        
        bond=0.d0 ! length of the Ca-Ca bond is averaged over all bonds
        do ic=1,nchains
            do i=menchain(ic)+1,menchain(ic+1)-1 
                dx=xn(i)-xn(i+1)
                dy=yn(i)-yn(i+1)
                dz=zn(i)-zn(i+1)
                dal=dx*dx+dy*dy+dz*dz
                dal=dsqrt(dal)
                b(i)=dal
                bond = bond + dal
            enddo
        enddo
        bond=bond/(men-1) ! this averaging previously was in gopotential
        
        write(1,'(a)')'#USING THE GO-LIKE 6-12 LJ POTENTIALS'
        if(lpdb) call gopotential(asigma)
        if(lwritego) then
            call print_cmap(22,0)
            close(1)
            close(2)
            close(22)
            stop
        endif
        
        call prepare(edsg)
        call evalgo(edsg,chi)
        if(lwal.and..not.lnowal) call evalwall(edsg)
        if(lpid) then
            call evalimproper(edsg,lcospid,epsbb)
        else
            call evalcpot(edsg)
        endif
        if(langle.or.lwritemap) call evalangles(edsg,lsldh,1.d0)
        if(lchiral) then
            call eval_chirality(enechi)
            edsg=edsg+enechi
        endif
        if(lwritemap) call print_map(22,0)
        if(lpullrel) lvelo=.TRUE.

        write(1,'(/,a,i10)') '#TOTAL PROTEIN LENGTH      ',men
        write(1,'(a,i10)')  '#NUMBER OF NATIVE CONTACTS ',klont
        if(lpdb) then
          write(1,'(a,f10.4)')'#AVERAGE LENGTH OF CONTACTS',asigma*unit
        else
          write(1,'(a)')'#NO PDB FILE USED FOR GO MODEL CONSTRUCTION'
        endif
        write(1,'(a,f10.2)')'#ENERGY OF NATIVE STATE    ',edsg
        write(1,*)
        if(lforce) write(1,'(a,f7.4)')'#USING AFM FORCE  ',coef
        if(lvelo.or.lforce) write(1,'(a,2f10.2)')
     +  '#PULLING SPRING CONSTANTS  ',HH1/unit/unit,HH2/unit**4
        if(naver.ne.0) 
     +  write(1,'(a,i8,a)') '#FORCE AVERAGED OVER ',naver,' TAU'
        write(1,'(a,f7.2)') '#VERLET LIST CUTOFF ',verlcut*unit
        if(loscillate) then
            write(1,'(a,f8.6,a,f9.1)')
     +      '#ANGULAR FREQUENCY ',omega,' PERIOD ',period
        else if(lvelo .OR. lwal) then
            write(1,'(a,f7.4)') '#CONSTANT VELOCITY',velo
        endif
        if(ldens) then
            write(1,'(a,f8.5)')'#SQEEZING VELOCITY',densvelo
            write(1,'(a,f8.5,a)')'#TARGET DENSITY ',
     +      tdens/(unit**3),' RESIDUES/A^3'
        endif
        if(lvelo) write(1,'(a,i10,a)') '#KWFORCE ',kwforce,' TAU'
        if(lwal) write(1,'(a,f7.4)')'#WALL POTENTIAL COEFFICIENT ',fwal
        if(lpullrel) then
        write(1,'(a)')'#STUDY PULLING AT CONSTANT VELOCITY THEN STOP'
        write(1,'(a)')'#PULLING AT CERTAIN DISTANCE'
        write(1,'(a,i10)')'#PULLING TIME [tau]        ',mpull
        else
        if(lmedian) write(1,'(a)')'#COMPUTING MEDIAN FOLDING TIMES'
        if(lthermo) write(1,'(a)')'#COMPUTING THERMODYNAMIC PROPERTIES'
        if(lunfold) write(1,'(a)')'#STUDYING UNFOLDING'
        endif
        if(lconftm) then
        if(lunfold) then
        write(1,'(a)')'#COMPUTING AVERAGED TIMES FOR CONTACT BREAKING'
        write(1,'(a)')'#AND AVERAGED LIFE TIMES OF CONTACTS'
        else
        write(1,'(a)')'#COMPUTING AVERAGED TIMES FOR CONTACT FORMATION'
        endif
        endif
        write(1,'(/,a,f10.3)')'#DELTA    =',delta
        write(1,'(a,f10.3)')'#GAMMA    =',gamma
        write(1,'(/,a,i10)') '#NUMBER OF TRAJECTORIES ',ntraj
        write(1,'(a,i10)') '#SIMULATION TIME        ',mstep
        write(1,'(a,i10)') '#SKIPPING STEPS         ',mskip
        write(1,'(a,i10)') '#EQUILIBRATION TIME     ',kteql
        write(1,'(a,i10)')'#RANDOM SEED            ',iseed
        write(1,'(/,a,7x,3f7.3)')'#TSTART TEND TSTEP',tstart,tend,tstep
        write(1,'(a,i6)')'#NUMBER OF TEMPERATURE STEPS',nt

        ! RESCALE TIME PARAMETERS
        kunit=nint(1.d0/delta) ! if delta = 0.005, this is 200
        naver=naver*kunit
        krst=krst*kunit
        mskip=mskip*kunit
        mstep=mstep*kunit
        ktrest=ktrest*kunit
        kteql=kteql*kunit
        nratvel=nratvel*kunit
        nratveld=nratveld*kunit
        kwrite=kwrite*kunit
        ksave=ksave*kunit
        kksave=kksave*kunit
        mpull=mpull*kunit
        kwforce=kwforce*kunit
        kwquarterperiod=nint(0.25*period)*kunit ! 1/4 of period
        omega=omega*delta
        ! SCALE FACTORS FOR VELOCITIES DURING EQUILIBRATION
        delsq=delta*delta
        deltsq=0.5d0*delsq
        ! SET PARAMETERS IN PREDICTOR-CORRECTOR METHOD
        f02=dble(3./16.)
        f12=dble(251./360.)
        f32=dble(11./18.)
        f42=dble(1./6.)
        f52=dble(1./60.)
        dfr=ran2(iseed) ! dfr is not used anywhere, just here (legacy)
        
C ===============================================
        ! LOOP OVER TEMPERATURES
        do 2000 it=1,nt
        if(lteql  .and. kteql.gt.0) then
            tr = teql
        write(1,'(/,a,f7.3,a,i7,a,f7.3)')
     +  '#Teql',tr,' for the first ',kteql/kunit,' Tau, then ',ttab(it)
        else
            tr = ttab(it)
            write(1,'(/,a,f7.3)')'#TEMPERATURE ',tr
        endif

        ! LANGEVIN PARAMETERS 
        gamma2=gamma/delta
        const2=2.d0*tr*gamma*delta     ! assume xmas=1
        const2=dsqrt(const2)*delta
        aheat=delsq*part*3*tr

        if(lconftm) then
        do i=1,klont
        fbt(i)=0.d0
        fbt2(i)=0.d0
        mtraj(i)=0
        enddo
        endif

        inot=0

        if(lthermo) then
        acv=0.d0
        acv2=0.d0
        apnat=0.d0
        apnat2=0.d0
        achi=0.d0
        achi2=0.d0
        endif

        nfm=mstep/kwforce
        if(lvelo.or.lwal) then
        do i=1,nfm ! averages over multiple trajectories
        aufres(i)=0.d0
        aufres2(i)=0.d0
        aree(i)=0.d0
        adfres(i)=0.d0
        adfres2(i)=0.d0
        enddo
        afmax1=0
        afmax2=0
        atmax1=0
        atmax2=0
        endif

C =========================================
        ! LOOP OVER STARTING CONFIGURATIONS
        iterate=0
        do 1000 itraj=1,ntraj
        
        if(lteql) then
            tr = teql
            const2=2.d0*tr*gamma*delta     ! assume xmas=1
            const2=dsqrt(const2)*delta
            aheat=delsq*part*3*tr
        endif

        ! LANGEVIN PARAMETERS 
        
        kbwal(1)=0 ! TIME COUNTER AFTER SQUEEZING
        kbwal(2)=0 ! T COUNTER BEFORE REACHING VELOCITY TO FIND MINIMUM
        kbwal(3)=0 ! T COUNTER AFTER REACHING MINIMUM
        kbwal(4)=0 ! T COUNTER BEFORE REACHING SQUEEZING VELOCITY
        kbwal(5)=0 ! T COUNTER AFTER REACHING MAXIMUM AMPLITUDE
        kbwal(6)=-2 ! T COUNTER OF OSCILLATIONS
        kbwal(7)=0 ! T COUNTER AFTER OSCILLATIONS
        kbwal(8)=0 ! T COUNTER BEFORE REACHING PULLING VELOCITY
        kbwal(9)=0 ! ALWAYS 0
        lcontin=.true.  !true: simulation proceeds; false: it stops
        vtarget=abs(vtarget)
        fresist = 0.d0
        afresist = 0.d0
        bfresist = 0.d0
        fresistperp = 0.d0
        afresistperp = 0.d0
        bfresistperp = 0.d0
        aufresist = 0.d0
        adfresist = 0.d0
        axufresist = 0.d0
        axdfresist = 0.d0
        ayufresist = 0.d0
        aydfresist = 0.d0
        bufresist=0.0
        bdfresist=0.0
        bxufresist = 0.d0
        bxdfresist = 0.d0
        byufresist = 0.d0
        bydfresist = 0.d0
        work=0.d0  !    WORK DONE IN ONE OSCILLATION CYCLE
        kbperiod=1 !    NUMBER OF OSCILLATION CYCLE
        intrsc=nssb!    NUMBER OF INTRACHAIN DISULFIDE BONDS
        intesc=0   !    NUMBER OF INTERCHAIN DISULFIDE BONDS
        icnss=0    !    NUMBER OF NATIVE DISULFIDE BONDS
        icdss=0    !    NUMBER OF NON-NATIVE DISULFIDE BONDS
        cofdih=0.0 !    COEFFICIENT FOR SLOWLY TURNING ON DIHEDRAL POT.
        shear=0.0  !    THE SYSTEM IS NOT SHEARED AT THE BEGINNING
        menw=0     !    NUMBER OF FCC WALL BEADS IS 0 AT THE BEGINNING
        kfccw=0    !    LENGTH OF THE FCC WALL CONTACT LIST IS ALSO 0
        kqont=0    !    NUMBER OF NON-NATIVE CONTACTS IS ALS0 0
        jq=1       !    VERLET LIST HAS 2 COPIES, JQ=1 AND JQ=2
        if(lpullrel) then
        lvelo=.TRUE.
        lunfold=.TRUE.
        lmedian=.FALSE.
        endif

        if(lmedian.and.(.not.lconftm).and.(.not.lthermo)) then
        if(inot.gt.ntraj/2+1) goto 1001
        endif
        iterate=iterate+1

        do i=1,klont
        kbt(i)=0
        kut(i)=0
        imap(i)=0
        enddo
        
        do i=1,4 ! first row of E is preparation
        youngmod(i,1)=0.d0
        enddo    ! work is the work required to place polymers together
        
        if(lthermo) then
c       initialize averages over the trajectories
        ave=0.d0
        ave2=0.d0
        pnat=0.d0
        bchi=0.d0
        bchi2=0.d0
        endif

        ! STARTING CONFORMATION
        if(lconftm.or.lmedian) then   ! A STRAIGHT-LINE 
            do i=1,men
                x0(i)=0.d0
                y0(i)=0.d0
                z0(i)=(i-1)*bond
            enddo
            if(lsawconftm) call confstart(sdens,confcut)
        else if(lwarmup) then
            if(it.eq.1.and.itraj.eq.1) then
                do i=1,men
                    x0(i)=xn(i)
                    y0(i)=yn(i)
                    z0(i)=zn(i)
                enddo
            endif
        else if(lpdb.or.lstartpdb) then        ! THE NATIVE STATE
            do i=1,men
                x0(i)=xn(i)
                y0(i)=yn(i)
                z0(i)=zn(i)
            enddo
        else
            call confstart(sdens,confcut)
            do i=1,men
                xn(i)=x0(i)
                yn(i)=y0(i)
                zn(i)=z0(i)
            enddo
        endif
        
        if(ldisp) call displace(iprota,jprotb,away)
        
        if(.not.lcpb) then
            startvolume=(men)/(sdens) !*unit**3)
            startboxsize=0.5*startvolume**(1.0/3.0)
            xmin=-startboxsize
            ymin=-startboxsize
            zmin=-startboxsize
            xmax=startboxsize
            ymax=startboxsize
            zmax=startboxsize
            do ib=1,men
                if(x0(ib)-2*bond.lt.xmin) xmin=x0(ib)-2*bond
                if(y0(ib)-2*bond.lt.ymin) ymin=y0(ib)-2*bond
                if(z0(ib)-2*bond.lt.zmin) zmin=z0(ib)-2*bond
                if(x0(ib)+2*bond.gt.xmax) xmax=x0(ib)+2*bond
                if(y0(ib)+2*bond.gt.ymax) ymax=y0(ib)+2*bond
                if(z0(ib)+2*bond.gt.zmax) zmax=z0(ib)+2*bond
            enddo
            xdown=xmin
            xup=xmax
            ydown=ymin
            yup=ymax
            zdown=zmin
            zup=zmax
            xsep=xup-xdown
            ysep=yup-ydown
            zsep=zup-zdown
            xinv=1.d0/xsep
            yinv=1.d0/ysep
            zinv=1.d0/zsep
        endif
        if(lwal) then ! CALCULATE WALL ENDS
          ip1=0
          ip2=0
          vpull=0.d0
          do ib=1,men
            z0temp(ib)=z0(ib)
            if(ljwal) z0temp(ib)=0
            ksorted(ib)=ib
            ipw(1,ib)=0
            ipw(2,ib)=ib
          enddo
          if(ldens.and..not.lcpb) then
            startvolume=(men)/(sdens) !*unit**3)
            startboxsize=0.5*startvolume**(1.0/3.0)
            xyzmin=-startboxsize
            xyzmax=startboxsize
            do ib=1,men ! to form a cube, wals must have same distance
                if(x0(ib)-2*bond.lt.xyzmin) xyzmin=x0(ib)-2*bond
                if(y0(ib)-2*bond.lt.xyzmin) xyzmin=y0(ib)-2*bond
                if(z0(ib)-2*bond.lt.xyzmin) xyzmin=z0(ib)-2*bond
                if(x0(ib)+2*bond.gt.xyzmax) xyzmax=x0(ib)+2*bond
                if(y0(ib)+2*bond.gt.xyzmax) xyzmax=y0(ib)+2*bond
                if(z0(ib)+2*bond.gt.xyzmax) xyzmax=z0(ib)+2*bond
            enddo
            xdown=xyzmin
            xup=xyzmax
            ydown=xyzmin
            yup=xyzmax
            zdown=xyzmin
            zup=xyzmax
            xsep=xup-xdown
            ysep=yup-ydown
            zsep=zup-zdown
            xinv=1.d0/xsep
            yinv=1.d0/ysep
            zinv=1.d0/zsep
          endif
        endif
        oldxup=xup
        oldxdown=xdown
        oldyup=yup
        oldydown=ydown
        oldzup=zup
        oldzdown=zdown
        
        if(lforce.or.lvelo) then ! COMPUTE THE DIRECTION OF FORCE
        afx=xn(ip2)-xn(ip1)
        afy=yn(ip2)-yn(ip1)
        afz=zn(ip2)-zn(ip1)
        aff=sqrt(afx*afx+afy*afy+afz*afz)
        afx=afx/aff
        afy=afy/aff
        afz=afz/aff
        vpulx=vpull*afx
        vpuly=vpull*afy
        vpulz=vpull*afz
        endif

        if(lvelo) then
        reemax=0.975*bond*men ! maximum length, after it trajectory ends
        xpul(1)=xn(ip2)
        ypul(1)=yn(ip2)
        zpul(1)=zn(ip2)
        if(lmaxforce) then
        fmax1=0
        fmax2=0
        tresh1=70  ! threshold for measuring maximum force for 1tit hard
        tresh2=160 ! numbers are in angstrems?
        tresh1=tresh1/velo  ! thresh1 and 2 are in taus
        tresh2=tresh2/velo
        endif
        endif

        ! LOAD INITIAL VELOCITIES OF PARTICLES
        call intvel3d(aheat,part,nen1)
        kb0=0
        sep0=0.0
        do i=1,men ! ZERO THE CONNECTION TABLES
        l3rdcn(i)=.false. ! true if residue i forms a disulfide bond
        knct3rd(i)=0        ! index of residue that is bonded to i
        khbful(1,i)=0
        khbful(2,i)=0
        khbful(3,i)=0
        khbful(4,i)=0
        nei(1,ib)=0
        nei(2,ib)=0
        enddo
        if(.not.ldynss) then
            do ks=1,nssb
                i1=ksb(1,ks)
                i2=ksb(2,ks)
                l3rdcn(i1)=.true.
                l3rdcn(i2)=.true.
                knct3rd(i1)=i2
                knct3rd(i2)=i1
            enddo
        endif
        if(itraj.gt.1) lrst=.false.
        if(lrst) then
            open(21,file=rstfile,status='unknown')
            read(21,*)time,work
            kb0=int(time/delta)
            mstep=mstep+kb0
            read(21,*)zup,zdown,yup,ydown,xup,xdown,sep0,shear
            if(lampstrict) amp=ampstrict/sep0
            read(21,*)(kbwal(i),i=1,8)
            state='E ' ! EQUILIBRATION
            if(kbwal(1).gt.0) vtarget=-1.0*vtarget
            if(loscillate.and.kbwal(3).gt.0) vtarget=-1.0*vtarget
            if(kbwal(kconnecttime).gt.0) then 
                if(lfcc) call make_fcc()
                if(kbwal(kconnecttime+1).eq.0) state='B ' ! AFTER BEADS
            endif
            if(loscillate.and.kbwal(6).gt.-2) then
                state='O '
                vtarget=-1.0*vtarget
            endif
            if(kbwal(8).gt.0) state='P '
            if(lfrmscrtch) then
                kteql=kb0+ktrest  ! EQUILIBRATION AFTER RESTART
c                call update_verlet_list(verlcut,nen1) ! only for ssbond
c                call compute_ssbonds()                ! now obsolete
                do i=1,men ! this is OK, positions are from the PDB file
                    x0(i)=x0(i)+xdown
                    y0(i)=y0(i)+ydown
                    z0(i)=z0(i)+zdown
                enddo
            else
                read(21,*)intrsc,intesc
                do k=1,2*(intrsc+intesc)
                    read(21,*)i,j
                    l3rdcn(i)=.true.
                    l3rdcn(j)=.true.
                    knct3rd(i)=j
                    knct3rd(j)=i
                enddo
                read(21,*)(ksorted(i),i=1,men)
                j=1
c                do k=1,ipwn/2
c                    i=ksorted(k)
c                    ipw(1,j)=0
c                    ipw(2,j)=i
c                    j=j+1
c                    i=ksorted(men+1-k)
c                    ipw(1,j)=1
c                    ipw(2,j)=i
c                    j=j+1
c                enddo
                read(21,*)(xpul(i),ypul(i),zpul(i),i=1,men)
                read(21,*)(x0(i),y0(i),z0(i),i=1,men)
                read(21,*)(x1(i),y1(i),z1(i),i=1,men)
                read(21,*)icnss,icdss,ip1,ip2
                if(lwal) then
                    read(21,*)(ipw(1,i),i=1,men)
                    read(21,*)(ipw(2,i),i=1,men)
                endif
                read(21,*)kbperiod
            endif
            close(21)
        endif
        
        xsep=xup-xdown
        ysep=yup-ydown
        zsep=zup-zdown
        xinv=1.d0/xsep
        yinv=1.d0/ysep
        zinv=1.d0/zsep
        call update_verlet_list(verlcut,nen1) ! set up Verlet list
        ! ASSIGN INITIAL ACCELERATION BASED ON INITIAL POSITIONS
        !if(lj1012) then
        !call evalgo_1012(epot)
        !else
        call prepare(epot)
        call evalgo(epot,chi)
        if(lwal.and..not.lnowal) call evalwall(epot)
        if(lpid) then
            call evalimproper(epot,lcospid,epsbb)
        else
            call evalcpot(epot)
        endif
        !endif
        if(lchiral) then
        call eval_chirality(enechi)
        epot=epot+enechi
        endif
        if(langle.or.lwritemap) call evalangles(epot,lsldh,0.d0)

        ! SCALE ACCELERATIONS
        do 530 i=1,men
            x2(i)=fx(i)*deltsq
            y2(i)=fy(i)*deltsq
            z2(i)=fz(i)*deltsq
530     continue

        ! TABLE HEADING
        if(ksave.ne.0) write(2,'(a,i9)')'MODEL',iterate
        if(kwrite.ne.0) then
        write(1,'(//,a,i4)')'#TRAJECTORY',iterate
        if(lforce) then
        write(1,'(a,a)')'#    TIME   AFORCE      EPOT      ETOT',
     +  '   ICN     RG    RMSD  D(1,N)   VEL'
        else if(lvelo) then
        write(1,'(a,a)')'#    TIME      EPOT      ETOT',
     +  '   ICN     RG    RMSD   D(1,N)   FORCE'
        else if(lwal) then
        if(lwals) then
        write(1,'(a,a,a,a)') '#S     TIME        EPOT        ETOT',
     +  ' INTRHC INTEHC INTRSC INTESC    ICN     RMSD     FZ_UP',
     +  '    FZDOWN       |F|       SEP    FPULLZ    FPULLX',
     +  '         W NCORD     FX_UP    FXDOWN     FY_UP    FYDOWN'
c        else if(lwals) then
c       write(1,'(a,a,a)')'#S    TIME      EPOT      ETOT   ICN     RG',
c     +  '    RMSD   FZ_UP  FZDOWN     |F|      SEP',
c     +  '   FX_UP  FXDOWN   FY_UP  FYDOWN   FPULL'
        else
        write(1,'(a,a,a,a)') '#S     TIME      EPOT        ETOT',
     +  ' INTRHC INTEHC INTRSC INTESC    ICN     RMSD     FZ_UP',
     +  '    FZDOWN       |F|       SEP    FPULLZ    FPULLX',
     +  '         W NCORD     XCM     YCM     ZCM   ICW'
        endif
        else if(ldynss) then
        write(1,'(a,a,a)')'#    TIME          EPOT          ETOT',
     +  '   ICN ICNss ICDss',
     +  '      RG       L    RMSD NCORD     W CORDR KNOTS KNOTE'
        else
        write(1,'(a,a,a)')'#    TIME          EPOT          ETOT',
     +  '   ICN B1-B2 S1-S2 B1-S2 B1-B1 S1-S1 B1-S1',
     +  '      RG       L    RMSD NCORD     W CORDR KNOTS KNOTE'
        endif
        endif

        if(lpullrel) then
        write(1,'(a)')'#START PULLING WITH CONSTANT VELOCITY'
        call flush(1)
        do i=1,mpull
        call lang(twopi,gamma2,const2,nen1)
        call corr(deltsq,nen1)
        call predct(nen1)
        call prepare(epot)
        call evalgo(epot,chi)
        if(lwal.and..not.lnowal) call evalwall(epot)
        if(lpid) then
            call evalimproper(epot,lcospid,epsbb)
        else
            call evalcpot(epot)
        endif
        if(lchiral) then
        call eval_chirality(enechi)
        epot=epot+enechi
        endif
        if(langle.or.lwritemap) then
        call evalangles(epot,lsldh,min(i*1.d0/ad,1.d0))
        endif
        call vafm(fresist,fresistperp,epot)
        enddo
        write(1,'(a,f12.2)')'#STOP PULLING AT D= ',mpull*delta*velo*unit
        lmedian=.TRUE.
        lunfold=.FALSE.
        lvelo=.FALSE.
        lforce=.FALSE.
        kb0=mpull
        endif

c -----------------------------------------
c        ENTER MAIN LOOP OF SIMULATION
c -----------------------------------------

        kb=kb0
533     continue
        kb=kb+1

        if(lmass) then
            call lang_mass(twopi,gamma2,const2,nen1)
        else
            call lang(twopi,gamma2,const2,nen1)
        endif
        call corr(deltsq,nen1)
        call predct(nen1)
        !if(lj1012) then
        !    call evalgo_1012(epot)
        !else
        call prepare(epot)
        call evalgo(epot,chi)
        if(lwal.and..not.lnowal) call evalwall(epot)
        if(lpid) then
            call evalimproper(epot,lcospid,epsbb)
        else
            call evalcpot(epot)
        endif
        !endif
        if(lchiral) then
            call eval_chirality(enechi)
            epot=epot+enechi
        endif
        if(langle.or.lwritemap) then
            if(lsldh.and.kb.le.ad) cofdih=kb*1.d0/ad
            call evalangles(epot,lsldh,cofdih)
        endif
        if(lmass) then
            do ib=1,men
                rma=rmas(ib)
                fx(ib)=fx(ib)/rma
                fy(ib)=fy(ib)/rma
                fz(ib)=fz(ib)/rma
            enddo
        endif
        if(kb.gt.kteql) then ! START OF A GLOBAL EQUILIBRATION CONDITION
        
        if(lwal) then
            totforce=zuforce+zdforce+fresist ! fresist<0, zforce>0
            btotforce=bufresist+bdfresist+bfresist !averaged total force
            if(ldens) then
                !---------PHASE 1: SQUEEZING PROTEINS TOGETHER----------
                if(kbwal(1).eq.0) then
                    if(xsep*ysep*zsep.le.2.0*targetvolume) then
                        if((kb-kteql).lt.nratveld) then !REACHING VPULL
                            vpull=vtarget*float(kb-kteql)/nratveld
                            state='S ' ! SQUEEZING
                        else
                            vpull=vtarget !SQUEEZING AND CHECKING VOLUME
                        endif
                        if(xsep*ysep*zsep.le.targetvolume) then
                            call print_restart(kb,itraj)
                            kbwal(1)=1
                            vtarget=-1.0*vtarget
                            vpull=0.d0
                            if(kconnecttime.eq.1) then
                                if(.not.ljwal) call connect_to_wal()
                                if(lfcc) call make_fcc()
                                state='B ' ! AFTER ATTACHING THE BEADS
                            else
                                state='E ' ! EQUILIBRATION
                            endif
                        endif
                    else
                        if((kb-kteql).lt.nratveld) then !REACHING VPULLD
                            vpull=vpulld*float(kb-kteql)/nratveld
                            state='S ' ! SQUEEZING
                        else
                            vpull=vpulld ! SQUEEZING AND CHECKING VOLUME
                        endif
                    endif
                !---------PHASE 2: RESTING AFTER SQUEEZING PROTEINS-----
                else if(kbwal(1).lt.ktrest) then
                    kbwal(1)=kbwal(1)+1
                !---------PHASE 3: EXTENDING THE BOX TO FIND MINIMUM----
                else if(kbwal(3).eq.0) then
                    bceilf=sign(min(abs(0.5*btotforce),1.d0),btotforce)
                    if(kbwal(2).lt.nratvel) then
                        vpull=vtarget*bceilf*float(kbwal(2))/nratvel
                        kbwal(2)=kbwal(2)+1
                        state='M ' ! SEARCHING FOR MINIMUM (F=0)
                    else ! RELAXATION SPEED IS PROPORTIONAL TO FORCE
                        vpull=vtarget*bceilf
                    endif
                    if(abs(btotforce).lt.(0.05+0.00005*men).or.
     +              .not.lminfbox) then
                        vpull=0.0
                        icheck=0
                        if(lobo) then
                            state='B ' ! AFTER ATTACHING THE BEADS
                        if(.not.ljwal) call connect_to_wal_one_by_one()
                            if(ip1+ip2.ge.ipwn) icheck=1
                        else
                            icheck=1
                        endif
                        if(icheck.eq.1) then
                         call print_restart(kb,itraj)
                         kbwal(3)=1 
                         if(lshear) then
                             sep0=xsep
                         else
                             sep0=zsep
                         endif
                         if(lampstrict) amp=ampstrict/sep0
                         vpull=0.d0
                         if(loscillate) vtarget=-1.0*vtarget
                         if(.not.lobo) then
                             if(kconnecttime.eq.3) then
                                 if(.not.ljwal) call connect_to_wal()
                                 if(lfcc) call make_fcc()
                                 state='B ' ! AFTER ATTACHING THE BEADS
                             else
                                 if(.not.lminfbox) kbwal(3)=ktrest
                                 state='E ' ! EQUILIBRATION
                             endif
                         endif
                        endif
                    endif
                !---------PHASE 4: RESTING AFTER FINDING MINIMUM--------
                else if(kbwal(3).lt.ktrest) then
                        kbwal(3)=kbwal(3)+1
                !---------PHASE 5: SQUEEZING BOX TO FULL AMPLITUDE------
                else if(loscillate.and.kbwal(5).eq.0) then
                    if(((.not.lshear).and.zsep.gt.sep0*(1.0-amp))
     +              .or.(lshear.and.shear.lt.sep0*amp*0.5)) then
                        if(kbwal(4).lt.nratvel) then
                            kbwal(4)=kbwal(4)+1
                            vpull=vtarget*float(kbwal(4))/nratvel
                            state='A ' ! REACHING MAXIMAL AMPLITUDE
                        else
                            vpull=vtarget
                        endif
                    else
                        call print_restart(kb,itraj)
                        kbwal(5)=1
                        vpull=0.0
                        if(kconnecttime.eq.5) then
                            if(.not.ljwal) call connect_to_wal()
                            if(lfcc) call make_fcc()
                            state='B ' ! AFTER ATTACHING THE BEADS
                        else
                            state='E ' ! EQUILIBRATION
                        endif
                    endif
                !---------PHASE 6: RESTING AFTER FINDING MAX AMPLITUDE--
                else if(loscillate.and.kbwal(5).lt.ktrest) then
                    kbwal(5)=kbwal(5)+1
                !---------PHASE 7: OSCILLATIONS-------------------------
                else if(loscillate.and.kbwal(6).eq.-2) then
                    kbwal(6)=-1
                    state='O ' ! OSCILLATIONS
                    vtarget=-1.0*vtarget
                else if(loscillate.and.kbperiod.le.kbperiodmax) then
                    kbwal(6)=kbwal(6)+1
                    vpull=-1.d0*sep0*amp*omega*sin(omega*kbwal(6))
                    if(mod(kbwal(6),kwquarterperiod).eq.0) then
                        kremainder=mod(kbwal(6)/kwquarterperiod,4)
                        if(kremainder.eq.0) then
                            works(kbperiod)=work
                            k=kbperiod
                            write(1,'(a,4f13.9,f10.2)') '#',
     +                      (youngmod(j,k)/unit**3,j=1,4),works(k)
                            work=0.d0
                            kbperiod=kbperiod+1
                            call print_restart(kb,itraj)
                            kbwal(6)=kbwal(6)+1 ! for the last check
                        endif
                        emd=btotforce/(amp*xsep*ysep)
                        youngmod(kremainder+1,kbperiod)=emd
                    endif
                !---------PHASE 8: RESTING AFTER OSCILLATIONS-----------
                else if(kbwal(7).lt.ktrest) then
                    klastquarter=mod(kbwal(6),kwquarterperiod)
                    if(loscillate.and.klastquarter.ne.0) then
                        kbwal(6)=kbwal(6)+1 !last check return to sep0
                        vpull=-1.d0*sep0*amp*omega*sin(omega*kbwal(6))
                        if(klastquarter.eq.0) then
                            emd=btotforce/(amp*xsep*ysep)
                            write(1,'(a,4f13.9,f10.2)') '#',emd/unit**3,
     +                      0.0,0.0,0.0,work
                            work=0.d0
                        endif
                    else
                        if(kconnecttime.eq.8 .and. kbwal(7).eq.0) then
                            if(.not.ljwal) call connect_to_wal()
                            if(lfcc) call make_fcc()
                            state='B ' ! AFTER ATTACHING THE BEADS
                        else
                            state='E ' ! EQUILIBRATION
                        endif
                        kbwal(7)=kbwal(7)+1
                        vpull=0.d0
                    endif
                !---------PHASE 9: PULLING------------------------------
                else if(vpull.gt.vtarget.and.lpullfin) then
                    vpull=vtarget*float(kbwal(8))/nratvel
                    kbwal(8)=kbwal(8)+1
                    state='P ' ! PULLING
                    if(lshear) sep0=zsep
                    kbwal(4)=0
                else if(lpullfin) then
                    vpull=vtarget
                else
                    vpull=0.0
                endif
            else !TODO
                if(vpull.lt.abs(vtarget)) then ! REACH VTARGET
                    if(vtarget.gt.0) then
                        state='S ' ! SQUEEZING
                        vpull=vtarget*float(kb-kteql)/nratvel
                    else if(kbwal(3).le.ktrest) then !EQUILIBRATION
                        state='R ' ! RESTING
                        kbwal(3)=kbwal(3)+1
                        vpull=0.d0
                    else
                        state='P ' ! PULLING
                        vpull=vtarget*float(kbwal(8))/nratvel
                        kbwal(8)=kbwal(8)+1
                    endif
                else
                    vpull=vtarget
                endif
            endif
            
            if(lshear.and.kbwal(4).gt.0) then
                shear=shear+vpull*0.5
                work=work+fresistperp*vpull
            else
                zup=zup-vpull*0.5
                zdown=zdown+vpull*0.5
                work=work+totforce*vpull ! if box shrinks, vpull > 0
                if(ldens.and.kbwal(1).eq.0) then
                    xup=xup-vpull*0.5
                    xdown=xdown+vpull*0.5
                    work=work+(yuforce+ydforce+xuforce+xdforce)*vpull
                    yup=yup-vpull*0.5
                    ydown=ydown+vpull*0.5 ! bad indent in the next line
          else if(lconstvol.and.kbwal(8).eq.0.and.kbwal(3).gt.0) then
                    vpullcv=-xsep*(1.d0-1.d0/sqrt(1.d0+vpull/zsep))
                    xup=xup-vpullcv*0.5
                    xdown=xdown+vpullcv*0.5
                    yup=yup-vpullcv*0.5
                    ydown=ydown+vpullcv*0.5
                    work=work+(yuforce+ydforce+xuforce+xdforce)*vpullcv
                endif
                xsep=xup-xdown
                ysep=yup-ydown
                zsep=zup-zdown
                xinv=1.d0/xsep
                yinv=1.d0/ysep
                zinv=1.d0/zsep
            endif
            lcontin=zsep.gt.sepmin
            if(kbwal(8).gt.0) then
                lcontin=lcontin.and.zsep.lt.reemax
            endif
            if((.not.ldens) .AND. (vtarget.gt.0)) then
                if((icw(1).ge.ipwn/2).and.(icw(2).ge.ipwn/2)) then
                    vtarget=-1.0*vtarget
                    if(.not.ljwal) call connect_to_wal()
                    if(lfcc) call make_fcc()
                endif
            endif
        endif
        
        do kwal=1,menw/2
            z0(men+kwal)=zdown
            z0(men+menw/2+kwal)=zup
        enddo
        
        if(lforce) then ! APPLY THE AFM FORCE
           aforce=coef*kb*delta
           call afm(aforce,epot)
        else if(lvelo .or. kbwal(kconnecttime).gt.0) then ! APPLY FORCE
           call vafm(fresist,fresistperp,epot)
           if(naver.ne.0) then 
            afresist=afresist+fresist
            afresistperp=afresistperp+fresistperp
            if(mod(kb,naver).eq.0) then
             bfresist=afresist/naver
             afresist=0.
             bfresistperp=afresistperp/naver
             afresistperp=0.
            endif
           endif
        endif
        
        else ! Continuation of equlibration condition
            state='E '
            if(kb.eq.kteql .and. kteql.gt.0) then ! last step of equil.
                tr = ttab(it) ! matters only if lteql is true
                const2=2.d0*tr*gamma*delta     ! assume xmas=1
                const2=dsqrt(const2)*delta
                aheat=delsq*part*3*tr
                if(lforce.or.lvelo) then
                    afx=x0(ip2)-x0(ip1)
                    afy=y0(ip2)-y0(ip1)
                    afz=z0(ip2)-z0(ip1)
                    aff=sqrt(afx*afx+afy*afy+afz*afz)
                    afx=afx/aff
                    afy=afy/aff
                    afz=afz/aff
                    vpulx=vpull*afx
                    vpuly=vpull*afy
                    vpulz=vpull*afz
                    if(lvelo) then
                        xpul(1)=x0(ip2)
                        ypul(1)=y0(ip2)
                        zpul(1)=z0(ip2)
                    endif
                    do i=1,men
                        xn(i)=x0(i)
                        yn(i)=y0(i)
                        zn(i)=z0(i)
                    enddo
                endif
            endif
        endif ! end of equilibration condition
        
        if(lwal) then ! MEASURE THE FORCE TO WALLS (A AND D)
           if(naver.ne.0) then
            aufresist=aufresist+zuforce
            adfresist=adfresist+zdforce
            if(mod(kb,naver).eq.0) then
             bufresist=aufresist/naver
             aufresist=0.
             bdfresist=adfresist/naver
             adfresist=0.
            endif
           if(lwals.or.(ldens.and.kbwal(1).eq.0.and..not.lsqpbc))then
             axufresist=axufresist+xuforce
             axdfresist=axdfresist+xdforce
             ayufresist=ayufresist+yuforce
             aydfresist=aydfresist+ydforce
             if(mod(kb,naver).eq.0) then
              bxufresist=axufresist/naver
              axufresist=0.
              bxdfresist=axdfresist/naver
              axdfresist=0.
              byufresist=ayufresist/naver
              ayufresist=0.
              bydfresist=aydfresist/naver
              aydfresist=0.
             endif
            endif
           endif
        endif
        
        if(lnatend.and. .not. lmedian) then
            lcontin=(icn.lt.klont)
            if(ldynss) lcontin=((icn+icnss).lt.(klont+nssb))
        endif
        if(kb.lt.mskip) goto 533             ! SKIPPING STEPS

        if(lthermo.and.(kb.eq.mskip)) then
c       initialize things for one trajectory after skipping
        endif

        ! CALCULATE KINETIC ENERGY AND MEAN COORDINATION NUMBER
        sumvel=0.d0
        do 540 i=1,men
            sumvel=sumvel+(x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i))*rmas(i)
540     continue
        ekin=sumvel/(2*delsq)           ! kinetic energy
        etot=epot+ekin                  ! total energy
        if(lgln) then
            ancord=(1.0*ncord)/ngln+2.0 ! (mean coord. num.) for GLN
        else
            ancord=(1.0*ncord)/men+2.0  ! (mean coord. num.) for all AAs
        endif
        ! FIND CONTACTS WHICH APPEAR FOR THE FIRST TIME DURING FOLDING
        if(lconftm) then
        do k=1,klont
        i=klist(1,k)
        j=klist(2,k)
        if(lunfold.and.imap(k).eq.1) kut(k)=kb ! CONTACT LAST APPERANCE
        if(kbt(k).eq.0) then
        if(lunfold) then
        if(imap(k).eq.0) kbt(k)=kb    ! CONTACT DISAPPEARS
        else
        if(imap(k).eq.1) kbt(k)=kb    ! CONTACT APPEARS
        endif
        endif
        enddo
        endif
        
        ! UPDATE VERLET LIST IF CERTAIN DISTANCE IS CROSSED
        vrcut2sqlim = vrcut2sq
        if(lpbcx) vrcut2sqlim = vrcut2sqlim-max(xdown-oldxdown,0.0)
        if(lpbcx) vrcut2sqlim = vrcut2sqlim-max(oldxup-xup,0.0)
        if(lpbcy) vrcut2sqlim = vrcut2sqlim-max(ydown-oldydown,0.0)
        if(lpbcy) vrcut2sqlim = vrcut2sqlim-max(oldyup-yup,0.0)
        if(lpbcz) vrcut2sqlim = vrcut2sqlim-max(zdown-oldzdown,0.0)
        if(lpbcz) vrcut2sqlim = vrcut2sqlim-max(oldzup-zup,0.0)
        do 1525 i=1,men
        if((x0(i)-oxv(1,i))**2+(y0(i)-oxv(2,i))**2+(z0(i)-oxv(3,i))**2
     +  .gt.vrcut2sqlim) then
            call update_verlet_list(verlcut,nen1)
            oldxup=xup
            oldxdown=xdown
            oldyup=yup
            oldydown=ydown
            oldzup=zup
            oldzdown=zdown
            goto 1526
        endif
1525    continue
1526    continue
        
        if(lvelo.and.mod(kb,kwforce).eq.0) then
        ifn=kb/kwforce
        dx=x0(ip1)-x0(ip2) ! END TO END DISTANCE
        dy=y0(ip1)-y0(ip2)
        dz=z0(ip1)-z0(ip2)
        ree=sqrt(dx*dx+dy*dy+dz*dz)
        aree(ifn)=aree(ifn)+ree
        if(naver.eq.0) then
        aufres(ifn)=aufres(ifn)+fresist
        aufres2(ifn)=aufres2(ifn)+fresist*fresist
        else
        aufres(ifn)=aufres(ifn)+bfresist
        aufres2(ifn)=aufres2(ifn)+bfresist*bfresist
        endif
        endif
        
        if(lwal.and.mod(kb,kwforce).eq.0) then
        ifn=kb/kwforce
        ! END TO END DISTANCE
        aree(ifn)=aree(ifn)+zsep ! aree measures avg distance of ends
        if(naver.eq.0) then      ! in this context avg are not relevant?
        aufres(ifn)=aufres(ifn)+zuforce
        aufres2(ifn)=aufres2(ifn)+zuforce*zuforce
        adfres(ifn)=adfres(ifn)+zdforce
        adfres2(ifn)=adfres2(ifn)+zdforce*zdforce
        else       ! the averages for aufres and adfres are over kwforce
        aufres(ifn)=aufres(ifn)+bufresist
        aufres2(ifn)=aufres2(ifn)+bufresist*bufresist
        adfres(ifn)=adfres(ifn)+bdfresist
        adfres2(ifn)=adfres2(ifn)+bdfresist*bdfresist
        endif ! force for output file (resist) is averaged over naver
        endif ! force for after-simulation statistis (res) over kwforce

        ! PRINTING TO FILE EVERY KWRITE STEPS AND AT THE END
        if(kwrite.ne.0) then
        if(kb.eq.1.or.mod(kb,kwrite).eq.0 .or. .not. lcontin) then
        ktime=nint(kb*delta)
        call gyration(rg)               ! RADIUS OF GYRATION
        call cgyration()                ! Rg and other parameters like W
        call compute_rmsd(rms)          ! RMSD

        if(.not.lwal) then
            dx=x0(ip1)-x0(ip2) ! END TO END DISTANCE
            dy=y0(ip1)-y0(ip2)
            dz=z0(ip1)-z0(ip2)
            ree=sqrt(dx*dx+dy*dy+dz*dz)
        endif

        icor=0
        corder=0.d0
        do k=1,kqont
            if(abs(kqist(4,k,jq)).gt.1) then
                corder=corder + abs(kqist(1,k,jq)-kqist(2,k,jq))
                icor=icor+1
            endif
        enddo
        if(icor.gt.0) corder=corder/icor
        corder=corder/men
        
        if(lvelo) then
        lcontin=ree.lt.reemax
c       projection of the restoring force on the direction of the force
        if(naver.eq.0) then
        write(1,'(i9,2f11.2,i6,f7.2,f8.3,f8.2,f9.2)')
     +  ktime,epot,etot,icn,rg*unit,rms*unit,ree*unit,fresist/unit
        else
        write(1,'(i9,2f11.2,i6,f7.2,f8.3,f8.2,f9.2)')
     +  ktime,epot,etot,icn,rg*unit,rms*unit,ree*unit,bfresist/unit
        if(lmaxforce) then
            if(ktime.le.nint(tresh1).and.bfresist.gt.fmax1) then
                fmax1=bfresist
                tmax1=dble(ktime)
            endif
            if(ktime.gt.nint(tresh1).and.ktime.le.nint(tresh2).and.
     +      bfresist.gt.fmax2) then
                fmax2=bfresist
                tmax2=dble(ktime)
            endif
        endif
        endif 
        go to 9997
        endif

        if(lwal) then
        ubr=bufresist/unit
        dbr=bdfresist/unit
        fbr=-1.d0*bfresist/unit
        fbrp=-1.d0*bfresistperp/unit
        if(lshear.and.kbwal(4).gt.0) then
            zse=shear*2.0*unit ! shear is applied to both wals, so 2x
        else
            zse=zsep*unit ! 0.5vpull is applied to both wals, so 1x
        endif
        if(lwals) then
            xubr=bxufresist/unit
            xdbr=bxdfresist/unit
            yubr=byufresist/unit
            ydbr=bydfresist/unit
            abr=sqrt((xubr-xdbr)**2+(yubr-ydbr)**2+(ubr-dbr)**2)
c'#S     TIME       EPOT       ETOT INTRHC INTEHC INTRSC INTESC       RG
c     RMSD     FZ_UP    FZDOWN       |F|       SEP    FPULLZ    FPULLX
c         W NCORD     FX_UP    FXDOWN     FY_UP    FYDOWN'
        write(1,'(a2,i9,2f11.2,5i7,f9.2,3f10.3,4f10.2,f6.2,4f10.2)')
     +      state,ktime,epot,etot,intrhc,intehc,intrsc,intesc,icn,
     +      rms*unit,ubr,dbr,abr,zse,fbr,fbrp,work,ancord,
     +      xubr,xdbr,yubr,ydbr
        else
            abr=(ubr+dbr)/2
        write(1,'(a2,i9,2f11.2,5i7,f9.2,3f10.3,4f10.2,f6.2,3f8.2,i6)')
     +      state,ktime,epot,etot,intrhc,intehc,intrsc,intesc,icn,
     +      rms*unit,ubr,dbr,abr,zse,fbr,fbrp,work,ancord,
     +      xmcm*unit,ymcm*unit,zmcm*unit,icw(1)+icw(2)
        endif
        go to 9997
        endif

        if(lforce) then
        dx=x0(ip1)-x0(ip2) ! END TO END DISTANCE
        dy=y0(ip1)-y0(ip2)
        dz=z0(ip1)-z0(ip2)
        ree=sqrt(dx*dx+dy*dy+dz*dz)
c       projection of the velocity on the direction of the force
        vl=x1(men)*afx + y1(men)*afy + z1(men)*afz
        vl=vl*unit/delta
        write(1,'(i9,f9.3,2f10.3,i6,f7.2,f8.3,f8.2,f8.2)')
     +  ktime,aforce,epot,etot,icn,rg*unit,rms*unit,ree*unit,vl
        go to 9997
        endif
        
        if(ldynss) then
        write(1,'(i9,2f14.3,3i6,2f8.2,f8.3,3f6.2,2i6)') ktime,epot,etot,
     +  icn,icnss,icdss,
     +  rg*unit,ree*unit,rms*unit,ancord,w(1),corder,knts(1,1),knts(2,1)
        go to 9997
        endif
        
        write(1,'(i9,2f14.3,7i6,2f8.2,f8.3,3f6.2,2i6)') ktime,epot,etot,
     +  icn,icnt(1),icnt(2),icnt(3),icnt(4),icnt(5),icnt(6),
     +  rg*unit,ree*unit,rms*unit,ancord,w(1),corder,knts(1,1),knts(2,1)

9997    continue
        call flush(1)

        endif
        endif
        if(ksave.ne.0.and.((kb.ge.(kteql).and.(mod(kb,ksave).eq.0)).or.
     +  (kb.eq.1.or.(kb.lt.(kteql).and.mod(kb,kksave).eq.0)))) then
! at first kteql time steps, sequence can be saved more often
            time=kb*delta
            if(mod(ksave,kwrite).ne.0) then
                call compute_rmsd(rms)   ! RMSD
                call cgyration() ! RG FOR INDIVIDUAL CHAINS
            endif
            if(lwritexyz) then
                call print_conf_xyz(2,time,epot,rms)
            else
            call print_conformation(2,time,epot,rms)
            endif
            if(lwritemap) call print_map(22,nint(kb*delta))
            call flush(2)
        endif

        if(krst.ne.0.and.mod(kb,krst).eq.0) then
            if(lcleanrst) then
                intrsc=0
                intesc=0
                icnss=0
                icdss=0
                do ib=1,men
                    khbful(1,ib)=0
                    khbful(2,ib)=0
                    khbful(3,ib)=0
                    khbful(4,ib)=0
                    l3rdcn(ib)=.false.
                    knct3rd(ib)=0
                enddo
                jq2=3-jq
                do k=1,kqont
                    kqist(3,k,jq)=0
                    kqist(4,k,jq)=sign(1,kqist(4,k,jq))
                    kqist(3,k,jq2)=0
                    kqist(4,k,jq2)=sign(1,kqist(4,k,jq2))
                enddo
            endif
            call print_restart(kb,itraj)
            if(ldelrst) then
              time=(kb-2*krst)*delta
              if(time.gt.0.d0) then  
                write(stafile,*) '(a',klenstr,',i2.2,i10.10,a)'
                write(rstfile,stafile) filname,itraj,nint(time),'.rst'
                open(38,iostat=irstat,file=rstfile,status='unknown')
                if (irstat.eq.0) close(38, status='DELETE')
              endif
            endif
        endif

        if(lthermo.and.kb.gt.kteql) then
            ave=ave+etot
            ave2=ave2+etot*etot
            if(icn.eq.klont) pnat=pnat+1.0
            bchi=bchi+chi
            bchi2=bchi2+chi*chi
        endif
        
        if(abs(etot).gt.99999999.9) lcontin=.false. ! stop if it blew up
        if(lrmsmax.and.rms*unit.gt.rmsmax) lcontin=.false. ! stop if rms
        
        if(lmedian) then
           if(icn.lt.klont.and.kb.lt.mstep+mskip) goto 533
        else
           if(kb.lt.mstep+mskip .AND. lcontin) goto 533
        endif

c -----------------------------------------
        ! END LOOP OF SIMULATION
c -----------------------------------------
544     continue

        if(ksave.ne.0.and.mod(kb,ksave).ne.0) then
        time=kb*delta
        call compute_rmsd(rms)       ! RMSD
        call cgyration()
        if(lwritexyz) then
            call print_conf_xyz(2,time,epot,rms)
        else
            call print_conformation(2,time,epot,rms)
        endif
        call flush(2)
        endif

        ! ACCUMULATE CONTACT BREAKING TIMES
        if(lconftm) then
        do k=1,klont
        if(kbt(k).ne.0) then
        aa=kbt(k)*delta
        fbt(k)=fbt(k)+aa
        fbt2(k)=fbt2(k)+aa*aa
        bb=kut(k)*delta
        fut(k)=fut(k)+bb
        fut2(k)=fut2(k)+bb*bb
        mtraj(k)=mtraj(k)+1
        endif
        enddo
        endif

        ! THE FOLDING TIME
        if(lmedian) then
        time=kb*delta
        tfold(iterate)=time
        if(icn.ne.klont) inot=inot+1
        if(kwrite.ne.0.and.kb.lt.mstep+mskip) then
        etot=epot+ekin                  ! total energy
        call gyration(rg)               ! RADIUS OF GYRATION
        call compute_rmsd(rms)          ! RMSD
        if(lforce) then
c       projection of the velocity on the direction of the force
        vl=x1(men)*afx + y1(men)*afy + z1(men)*afz
        vl=vl*unit/delta
        write(1,'(f9.1,f9.3,2f10.3,i6,f7.2,f8.3,f8.2,f8.2)')
     +  time,aforce,epot,etot,icn,rg*unit,rms*unit,ree*unit,vl
        else
        write(1,'(f9.1,2f10.3,i6,f7.2,f8.3)')
     +  time,epot,etot,icn,rg*unit,rms*unit
        endif
        call flush(1)
        endif
        endif
        
        if(loscillate) then
        write(1,'(//,a,x,5(9x,a))') '#','E1','E2','E3','E4','W'
        do k=1,kbperiod-1
        write(1,'(a,4f11.5,f10.2)') 
     +  '# ',(youngmod(j,k)/unit**3,j=1,4),works(k)
        enddo
        write(1,'(a,f11.2)') '# ',work
        endif
        
        if(lthermo) then
        tav=kb-kteql-mskip
        ave=ave/tav
        ave2=ave2/tav
        pnat=pnat/tav
        cv=(ave2-ave*ave)/(tr*tr)       ! SPECIFIC HEAT
        apnat=apnat+pnat
        apnat2=apnat2+pnat*pnat
        acv=acv+cv
        acv2=acv2+cv*cv
        bchi=bchi/tav
        bchi2=bchi2/tav
        dchi=bchi2-bchi*bchi     ! STRUCTURAL SUSCEPTIBILITY
        achi=achi+dchi
        achi2=achi2+dchi*dchi
        endif
        ! nfm is maximal 'time' for which all trajectories are recorded
        if(lvelo.and.(kb/kwforce).lt.nfm) nfm = kb/kwforce
        if(lvelo.and.lmaxforce) then
        write(1,*)'#the first peak '
        write(1,*)'#max force ',fmax1/unit,'  displacement ',tmax1*velo
        write(1,*)'#the second peak '
        write(1,*)'#max force ',fmax2/unit,'  displacement ',tmax2*velo
        afmax1=afmax1+fmax1/unit
        afmax2=afmax2+fmax2/unit
        atmax1=atmax1+tmax1*velo
        atmax2=atmax2+tmax2*velo
        endif
c       print distances in the contacts

1000    continue
        ! END LOOP OVER CONFIGURATIONS

1001    continue

        if(lconftm) then
        if(lunfold) then
        write(1,'(/,a)')'#AVERAGE TIME NEEDED FOR BREAKING EACH CONTACT'
        write(1,'(2a)')'# ICN    I    J  J-I       t0    DISP.',
     +  '       t1    DISP.'
        else
        write(1,'(/,a)')'#AVERAGE TIME NEEDED FOR FORMING EACH CONTACT'
        write(1,'(a)')'# ICN    I    J  J-I       t0    DISP.'
        endif
        do k=1,klont
        mntraj=mtraj(k)
        if(mntraj.ne.0) then
        fbt(k)=fbt(k)/mntraj
        fbt2(k)=fbt2(k)/mntraj
        fut(k)=fut(k)/mntraj
        fut2(k)=fut2(k)/mntraj
        endif
        aa=fbt(k)
        bb=sqrt(abs(fbt2(k)-aa*aa))/2
        i=klist(1,k)
        j=klist(2,k)
        if(lunfold) then
        cc=fut(k)
        dd=sqrt(abs(fut2(k)-cc*cc))/2
        write(1,'(4i5,4f11.2)')k,iseq(i),iseq(j),j-i,aa,bb,cc,dd
        else
c        write(1,'(4i5,2f11.2)')k,iseq(i),iseq(j),j-i,aa,bb
        caac(k)=aa ! time of contact break
        kaak(k)=k  ! index of contact
        endif
        enddo
            call sort2(klont,caac,kaak)
            do k=1,klont
            fcaa=caac(k)
            if(fcaa.gt.0) then
            ll=kaak(k)
            i=klist(1,ll)
            j=klist(2,ll)
            ise=iseq(i)
            jse=iseq(j)
            if(fcaa.gt.0.01) write(1,'(4i5,2f11.2)')ll,ise,jse,j-i,fcaa
            endif
            enddo
        endif

        ! THE MEDIAN FOLDING TIME
        if(lmedian) then
        call sort(iterate,tfold)
        tmed=tfold(ntraj/2+1)
        write(1,'(/,a)')'#  TEMP     TMEDIAN  NTRAJ  INOT'
        write(1,'(a,f6.2,f12.2,i7,i6)')'#',tr,tmed,iterate,inot
        endif

        if(lthermo) then
        pnat=apnat/iterate
        pnat2=apnat2/iterate
        pnat2=sqrt(abs(pnat2-pnat*pnat))/2.d0
        cv=acv/iterate
        cv2=acv2/iterate
        cv2=sqrt(abs(cv2-cv*cv))/2.d0
        dchi=achi/iterate
        dchi2=achi2/iterate
        dchi2=sqrt(abs(dchi2-dchi*dchi))/2.d0
        write(1,'(/,a,a)')
     +  '#  TEMP      PNAT    DISP/2           CV       DISP/2',
     +  '         CHI      DISP/2'
        write(1,'(a,f6.2,2f10.5,2f13.3,2f12.4)')
     +  '#',tr,pnat,pnat2,cv,cv2,dchi,dchi2
        endif


        if(lvelo.and.iterate.gt.1) then
        write(1,'(//,a)')
     +  '#AVERAGED FORCE ON STRETCHING WITH CONSTANT VELOCITY'
        write(1,'(a)') '#    TIME  <D(1,N)>  <FORCE>  DISP./2' 
        do i=1,nfm
        af1=aufres(i)/iterate
        af2=aufres2(i)/iterate-af1*af1
        af2=sqrt(abs(af2))*0.5d0
        ree=aree(i)/iterate*unit
        write(1,'(f9.1,f10.2,2f9.2)')i*kwforce*delta,ree,af1,af2
        enddo
        if(lmaxforce) then
        write(1,*)' '
        write(1,*)'#average over trajectories '
        at=1./ntraj
        write(1,*)'#the first peak '
        write(1,*)'#max force ',afmax1*at,'  displacement ',atmax1*at
        write(1,*)'#the second peak '
        write(1,*)'#max force ',afmax2*at,'  displacement ',atmax2*at
        endif
        endif

2000    continue
        ! END LOOP OVER TEMPERATURES
C ===============================================

        close(1)
c        write(2,'(a)') 'ENDMDL'
        close(2)
        close(22)
c ----------------------------
c         END OF EXECUTION
c ----------------------------
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine predct(nen1)
        implicit double precision(a-h,o-z)

        ! USE FIFTH-ORDER TAYLOR SERIES TO PREDICT POSITIONS & THEIR
        ! DERIVATIVES AT NEXT TIME-STEP

        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
!$omp parallel default(shared), private(i)
!$omp do
        do 300 i=nen1,men
        x0(i)=x0(i)+x1(i)+x2(i)+x3(i)+x4(i)+x5(i)
        y0(i)=y0(i)+y1(i)+y2(i)+y3(i)+y4(i)+y5(i)
        z0(i)=z0(i)+z1(i)+z2(i)+z3(i)+z4(i)+z5(i)
        x1(i)=x1(i)+2.d0*x2(i)+3.d0*x3(i)+4.d0*x4(i)+5.d0*x5(i)
        y1(i)=y1(i)+2.d0*y2(i)+3.d0*y3(i)+4.d0*y4(i)+5.d0*y5(i)
        z1(i)=z1(i)+2.d0*z2(i)+3.d0*z3(i)+4.d0*z4(i)+5.d0*z5(i)
        x2(i)=x2(i)+3.d0*x3(i)+6.d0*x4(i)+10.d0*x5(i)
        y2(i)=y2(i)+3.d0*y3(i)+6.d0*y4(i)+10.d0*y5(i)
        z2(i)=z2(i)+3.d0*z3(i)+6.d0*z4(i)+10.d0*z5(i)
        x3(i)=x3(i)+4.d0*x4(i)+10.d0*x5(i)
        y3(i)=y3(i)+4.d0*y4(i)+10.d0*y5(i)
        z3(i)=z3(i)+4.d0*z4(i)+10.d0*z5(i)
        x4(i)=x4(i)+5.d0*x5(i)
        y4(i)=y4(i)+5.d0*y5(i)
        z4(i)=z4(i)+5.d0*z5(i)
300     continue
!$omp enddo nowait
!$omp end parallel 
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine corr(deltsq,nen1)
        implicit double precision(a-h,o-z)

        ! CORRECT PREDICTED POSITIONS AND THEIR DERIVATIVES

        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/parm/f02,f12,f32,f42,f52
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
!$omp parallel default(shared), private(i,xerr,yerr,zerr)
!$omp do
        do 670 i=nen1,men
          xerr=x2(i)-deltsq*fx(i)
          yerr=y2(i)-deltsq*fy(i)
          zerr=z2(i)-deltsq*fz(i)
          x0(i)=x0(i)-xerr*f02
          x1(i)=x1(i)-xerr*f12
          x2(i)=x2(i)-xerr
          x3(i)=x3(i)-xerr*f32
          x4(i)=x4(i)-xerr*f42
          x5(i)=x5(i)-xerr*f52
          y0(i)=y0(i)-yerr*f02
          y1(i)=y1(i)-yerr*f12
          y2(i)=y2(i)-yerr
          y3(i)=y3(i)-yerr*f32
          y4(i)=y4(i)-yerr*f42
          y5(i)=y5(i)-yerr*f52
          z0(i)=z0(i)-zerr*f02
          z1(i)=z1(i)-zerr*f12
          z2(i)=z2(i)-zerr
          z3(i)=z3(i)-zerr*f32
          z4(i)=z4(i)-zerr*f42
          z5(i)=z5(i)-zerr*f52
670       continue
!$omp enddo nowait
!$omp end parallel 
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE PREPARES VARIABLES FOR COMPUTING ENERGY AND FORCE

        subroutine prepare(epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/xyforces/xuforce,xdforce,yuforce,ydforce
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
        
        epot=0.d0
        ncord=0
        icn=0
        !ncnt=0
        intrhc=0
        intehc=0
        icw(1)=0
        icw(2)=0
        do i=1,6
            icnt(i)=0
        enddo
        
        zuforce=0.d0
        zdforce=0.d0
        yuforce=0.d0
        ydforce=0.d0
        xuforce=0.d0
        xdforce=0.d0
!$omp parallel default(shared), private(i,j)
!$omp do
        do i=1,men
            nei(2,i)=nei(1,i)
            nei(1,i)=0
            fx(i) = 0.d0
            fy(i) = 0.d0
            fz(i) = 0.d0
            j=i+1
            v(1,i)=x0(j)-x0(i)
            v(2,i)=y0(j)-y0(i)
            v(3,i)=z0(j)-z0(i)
            j=i+2
            v(4,i)=x0(i)-x0(j)
            v(5,i)=y0(i)-y0(j)
            v(6,i)=z0(i)-z0(j)
        enddo
!$omp enddo nowait
!$omp end parallel
!$omp parallel default(shared), private(i,j,vxvnrm)
!$omp do
        do i=1,men-1
            j=i+1
            vxv(4,j)=v(3,j)*v(5,i)-v(2,j)*v(6,i)
            vxv(5,j)=v(1,j)*v(6,i)-v(3,j)*v(4,i)
            vxv(6,j)=v(2,j)*v(4,i)-v(1,j)*v(5,i)
            vxv(1,j)=v(3,i)*v(2,j)-v(2,i)*v(3,j)
            vxv(2,j)=v(1,i)*v(3,j)-v(3,i)*v(1,j)
            vxv(3,j)=v(2,i)*v(1,j)-v(1,i)*v(2,j)
            vxvnrm=vxv(1,j)*vxv(1,j)+vxv(2,j)*vxv(2,j)+vxv(3,j)*vxv(3,j)
            vxvnrm=dsqrt(vxvnrm)
            vxv(1,j)=vxv(1,j)/vxvnrm
            vxv(2,j)=vxv(2,j)/vxvnrm
            vxv(3,j)=vxv(3,j)/vxvnrm
            vnrm(j)=vxvnrm
        enddo
!$omp enddo nowait
!$omp end parallel
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  THIS SUBROUTINE COMPUTE THE ENERGY AND FORCE OF THE WALL INTERACTIONS

        subroutine evalwall(epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconftm,lwal,lwals,ldens,lcpot,ljwal
        logical lpbcx,lpbcy,lpbcz,lwalup,lwaldown,lsqpbc,lfcc
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/equil/kteql,ktrest,nratvel,nratveld,kconnecttime
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
        common/xyforces/xuforce,xdforce,yuforce,ydforce
        

        if(ljwal) then
            do i=1,men
                zi=z0(i) ! force from the Z plate(s)
                dz=(zi-zdown)
                uz=(zup-zi) 
                lwalup=uz.le.walmindst
                lwaldown=dz.le.walmindst
                if(lwalup.or.lwaldown) then
                    if(lwaldown) r=dz
                    if(lwalup) r=-uz
                    rsi=sigma1(9)/r
                    r6=rsi**6
                    if(z0temp(i).lt.1 .and..not.lfcc) then
                        coefpul = fwal
                        if(kbwal(kconnecttime).ne.0) then
                            k=1
737                             continue
                            if(k.le.ip1+ip2) then
                              if(ipw(1,k).ne.0) then
                                k=k+1
                                goto 737
                              endif
                            else
                             if(lwaldown.and.ip1.lt.men/2) ip1=ip1+1
                               if(lwalup.and.ip2.lt.men/2) ip2=ip2+1
                            endif
                            ipw(2,k)=i
                            if(lwaldown) ipw(1,k)=-1
                            if(lwalup) ipw(1,k)=1
                            z0temp(i)=2
                            xpul(k)=x0(i)
                            ypul(k)=y0(i)
                            zpul(k)=0.d0
                        endif
                    else
                        if(lwaldown) icw(1)=icw(1)+1
                        if(lwalup) icw(2)=icw(2)+1
                        coefpul = fwal !*(ad-z0temp(i))/ad
                    endif
                    epot=epot+coefpul*(4.d0*r6*(r6-1.d0)+1.d0)
                    fce= coefpul*24.d0*r6*(1.d0-2.d0*r6)/r ! sign OK
                    !if(fce.gt.1.d+3) fce=1.d+3
                    !if(fce.lt.-1.d+3) fce=-1.d+3
                    fz(i) = fz(i) - fce ! instead of +(-fce*dz/r)
                    if(lwaldown) zdforce = zdforce - fce
                    if(lwalup) zuforce = zuforce + fce
                endif
            enddo
        else
            if(kbwal(kconnecttime).eq.0) then
                kipwn1=0
                kipwn2=0
            else
                kipwn1=ip1
                kipwn2=ip2
            endif
            do ii=1+kipwn1,men-kipwn2
                if(.not.lpbcz) then
                    i=ksorted(ii)
                    zi=z0(i) ! force from the Z plate(s)
                    dz=(zi-zdown)
                    uz=(zup-zi) 
                    if(dz.le.walmindst) icw(1)=icw(1)+1
                    if(uz.le.walmindst) icw(2)=icw(2)+1
                    zzdowni=fwal*dz**(-10.0)
                    zdforce=zdforce + zzdowni
                    zzupi=fwal*uz**(-10.0)
                    zuforce=zuforce+zzupi
                    fz(i) = fz(i) + zzdowni - zzupi
                endif
      if(lwals.or.(ldens.and.kbwal(1).eq.0.and..not.lsqpbc)) then
                    xi=x0(i) ! force from the X plate(s)
                    xxdowni=fwal*(xi-xdown)**(-10.0)
                    xdforce=xdforce + xxdowni
                    xxupi=fwal*(xup-xi)**(-10.0)
                    xuforce=xuforce+xxupi
                    fx(i) = fx(i) + xxdowni - xxupi
                    yi=y0(i) ! force from the Y plate(s)
                    yydowni=fwal*(yi-ydown)**(-10.0)
                    ydforce=ydforce + yydowni
                    yyupi=fwal*(yup-yi)**(-10.0)
                    yuforce=yuforce+yyupi
                    fy(i) = fy(i) + yydowni - yyupi
                endif
            enddo
        endif
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTE THE ENERGY AND FORCE OF THE CUSTOM POTENTIAL

        subroutine evalcpot(epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lconftm,lwal,lwals,ldens,lcpot,lelectr,ljwal,lecperm
        logical lsimpang,lfrompdb(len),lconect(len),l3rdcn(len),lsim
        logical lwritemap,lsink,lel,lmrs,ldet,lki,lkj,lpid,lsqpbc,ltr
        logical ldynss,lss,lpbcx,lpbcy,lpbcz,lssn,lsselj,lrepcoul
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/restr/bckbmin,bckb2min,sdchnmax,sfact,kmaxbckb,lecperm
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
        
!$omp parallel default(shared), private(i,j,k,l,dx,dy,dz,r,rsq,icm,rsig,
!$omp& ene,fce,repx,repy,repz,kqistabs,iconttype,bbir,bbjr,bbij,lelectr,
!$omp& vxvinrm,vxvjnrm,sdchni,sdchnj,i0,j0,sdir,sdirnrm,sdjrnrm,screenr,
!$omp& lel,adiab,sdjr,r6,enec,expc,coulpotcoeff,totcoeff,phi,rb,rb2,rsi,
!$omp& ksdchnsi,ksdchnsj,kmaxbckbi,kmaxbckbj,iconttype2,kremainder,it2,
!$omp& jt2,itype,jtype,lki,lkj,icheck,lssn,kss,lss,kqadabs)
!$omp do reduction(+:epot,intehc,intrhc,intesc,intrsc,ncord,icnss,icdss)
        do 466 k=1,kqont ! non-native attractive contacts
            lss=.false.
            kqistabs=abs(kqist(4,k,jq))
            kqadabs=abs(kqist(3,k,jq))
            kremainder=kqistabs/kqist(4,k,jq) ! 1 for same protein
            if(kqistabs.gt.8) then
                iconttype=5
            else
                iconttype=kqistabs ! 4-8 for backbone contact
            endif
            rsig=sigma1(kqistabs) ! wrong if kqistabs<4, but w no effect
c            if(abs(j-i).eq.4) rsig=sigma1(8)
            iconttype2=1 ! type of contact after checking conditions
            i=kqist(1,k,jq) ! kqist 1 or 2 = residue numbers
            j=kqist(2,k,jq) ! jq = number of verlet list (1 or 2)
            itype=ksdchns(i)
            jtype=ksdchns(j)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.lt.rneisq) then
                nei(1,i)=nei(1,i)+1
                nei(1,j)=nei(1,j)+1
            endif
            if(rsq.gt.rcutsq) then
            if(kqadabs.ne.0) kqist(3,k,jq)=sign(kqadabs-1,kqist(3,k,jq))
                goto 465
            endif
            r = sqrt(rsq)
            lelectr=(lcpot.and.itype+jtype.gt.7)
            if(iconttype.lt.2) then
            ! ---------- contact not formed ----------------------------
            if(r.lt.sfact*(max(sigma1(inameseq(i)*21+inameseq(j)),
     +      sigma1(6)))) then ! this if is only for speed optimization
            bbir=0.0 ! bb = backbone, ir = criterion for r and residue i
            bbjr=0.0
            kmaxbckbi=kmaxbckb
            if(inameseq(i).eq.2) kmaxbckbi=1  ! pro has only 1 hbond
            if(khbful(1,i).lt.kmaxbckbi) then ! cos(n,r)
                bbir=abs((vxv(1,i)*dx+vxv(2,i)*dy+vxv(3,i)*dz)/r)
            endif
            kmaxbckbj=kmaxbckb
            if(inameseq(j).eq.2) kmaxbckbj=1  ! pro has only 1 hbond
            if(khbful(1,j).lt.kmaxbckbj) then ! cos(n,r)
                bbjr=abs((vxv(1,j)*dx+vxv(2,j)*dy+vxv(3,j)*dz)/r)
            endif
           ! if(bbir+bbjr.gt.2.d0*bckb2min) then
            if(bbir.gt.bckb2min .and. bbjr.gt.bckb2min) then
             if((abs(j-i).ne.4).or.kremainder.ne.1) then
              bbij=vxv(1,i)*vxv(1,j)+vxv(2,i)*vxv(2,j)+vxv(3,i)*vxv(3,j)
              if(abs(bbij).gt.bckbmin) iconttype2=4 ! backbone-backbone
             endif
            endif
           !   if(abs(bbij)+bbir+bbjr.gt.3*bckbmin) iconttype2=4
           ! endif
            if(iconttype2.lt.2) then
              sdchni=1.0 ! sdchn or sd = sidechain
              ksi=khbful(2,i)
              if(nei(2,i).lt.neimin) ksi=ksi+1
              if((i.gt.1).and.ksi.lt.ksdchn(inameseq(i),2))then
                i0=i-1 ! khbful 1-bb contacts,2-sd cont.,3-max sd cont
        sdir=(v(1,i0)-v(1,i))*dx+(v(2,i0)-v(2,i))*dy+(v(3,i0)-v(3,i))*dz
                sdirnrm=0.0
                if(sdir.lt.0.0) then ! if <0, no need to compute norm
                    sdchni=0.0
                else
                    do l=1,3
                        sdirnrm=sdirnrm+(v(l,i)-v(l,i0))**2
                    enddo
                    sdir=sdir**2/(sdirnrm*rsq) ! cos^2(n,r)
                    sdchni=sdir!min(1.0,max(sdchnmax-sdir,0.0)/sdchnmin)
                endif
              endif
              if(bbjr.gt.bckb2min .and. sdchni.lt.sdchnmax) then
                iconttype2=7 ! sidechain-backbone
              else
                sdchnj=1.0
                ksj=khbful(2,j)
                if(nei(2,j).lt.neimin) ksj=ksj+1 
                if((j.gt.1).and.ksj.lt.ksdchn(inameseq(j),2)) then
                  j0=j-1
        sdjr=(v(1,j)-v(1,j0))*dx+(v(2,j)-v(2,j0))*dy+(v(3,j)-v(3,j0))*dz
                  sdjrnrm=0.0
                  if(sdjr.lt.0.0) then
                    sdchnj=0.0
                  else
                    do l=1,3
                        sdjrnrm=sdjrnrm+(v(l,j)-v(l,j0))**2
                    enddo
                    sdjr=sdjr**2/(sdjrnrm*rsq)
                    sdchnj=sdjr!min(1.0,max(sdchnmax-sdjr,0.0)/sdchnmin)
                  endif
                endif
                if(bbir.gt.bckb2min .and. sdchnj.lt.sdchnmax) then
                  iconttype2=6 ! backbone-sidechain
                endif
              endif
            endif
            if(iconttype2.lt.2) then
              if(sdchni.lt.sdchnmax.and.sdchnj.lt.sdchnmax) then
                it2=2+min(itype,2)
                jt2=2+min(jtype,2)
                lki=khbful(jt2,i).lt.ksdchn(inameseq(i),jt2)
                lkj=khbful(it2,j).lt.ksdchn(inameseq(j),it2)
                lel=(.not.(lelectr.and.itype.eq.jtype))
         ltr=(inameseq(i).eq.20.and.inameseq(j).eq.20.and.abs(j-i).eq.3)
                !ksdchnsi=min(2,ksdchns(i))
                !ksdchnsj=min(2,ksdchns(j))
                if(lel.and.lki.and.lkj.and..not.ltr) iconttype2=5 !sd-sd
                !if(ksdchnsi.eq.ksdchnsj) iconttype2=5 ! sidechain-sd
              endif
            endif
            if(lsim.or.iconttype2.ne.5) then
                kqistabs=iconttype2
            else
                kqistabs=inameseq(i)*21+inameseq(j)
            endif
            endif ! end of the speed optimization if
            if(r.le.sigma1(kqistabs)*sfact.and.iconttype2.gt.1) then 
                kqist(4,k,jq)=sign(kqistabs,kqist(4,k,jq))
                if(kqadabs.lt.ad) kqist(3,k,jq)=kqadabs+1
                if(iconttype2.eq.5) then ! sidechain-sidechain
                  khbful(2,i)=khbful(2,i)+1
                  khbful(2,j)=khbful(2,j)+1
                  khbful(jt2,i)=khbful(jt2,i)+1
                  khbful(it2,j)=khbful(it2,j)+1
                  if(inameseq(i).eq.4 .and.inameseq(j).eq.4) then
                    lssn=(nei(2,i)+nei(2,j)).lt.neimaxdisul
                    if(lssn.and..not.l3rdcn(i).and..not.l3rdcn(j)) then
                      l3rdcn(i)=.true.
                      l3rdcn(j)=.true.
                      knct3rd(i)=j
                      knct3rd(j)=i
                      if(kqist(4,k,jq).gt.0) then    ! are the same type
                        intrsc=intrsc+1
                      else
                        intesc=intesc+1
                      endif
                      if(ldynss) then
                        icheck=0
                        do kss=1,nssb
                          if(ksb(1,kss).eq.i.and.ksb(2,kss).eq.j) then
                            icheck=1
                            icnss=icnss+1
                          endif
                          if(ksb(2,kss).eq.i.and.ksb(1,kss).eq.j) then
                            icheck=1
                            icnss=icnss+1
                          endif
                        enddo
                        if(icheck.eq.0) icdss=icdss+1
                      endif
                    endif
                  endif
                else if(iconttype2.eq.7) then ! sidechain-backbone
                    khbful(2,i)=khbful(2,i)+1
                    khbful(1,j)=khbful(1,j)+1
                    !khbful(4,j)=khbful(4,j)+1
                else if(iconttype2.eq.6) then ! backbone-sidechain
                    khbful(2,j)=khbful(2,j)+1
                    khbful(1,i)=khbful(1,i)+1
                    !khbful(4,i)=khbful(4,i)+1
                else if(iconttype2.eq.4) then ! backbone-backbone
                    khbful(1,i)=khbful(1,i)+1
                    khbful(1,j)=khbful(1,j)+1
                endif ! kqist 4 = type of contact
            else
              iconttype2=1
            if(kqadabs.ne.0) kqist(3,k,jq)=sign(kqadabs-1,kqist(3,k,jq))
            endif
            else
            ! ---------- contact already formed ------------------------
                if(r.le.rsig*cntfct) then
                   ncord=ncord+2
                   !if(lwritemap) then
                   !    ncnt=ncnt+1
                   !    nli(1,ncnt)=i
                   !    nli(2,ncnt)=j
                   !    nli(3,ncnt)=0
                   !    nli(4,ncnt)=iconttype
                   !    nli(5,ncnt)=kqist(4,k,jq)
                   !endif
                   iconttype2=min(iconttype,6) ! bb-sd (6) and sd-bb (7)
                   if(kremainder.eq.1) then    ! are the same type
                       intrhc=intrhc+1
                       icnt(iconttype2)=icnt(iconttype2)+1
                   else
                       intehc=intehc+1
                       icnt(iconttype2-3)=icnt(iconttype2-3)+1
                   endif
                   if(kqadabs.lt.ad) kqist(3,k,jq)=kqadabs+1
                else
                    if(kqist(3,k,jq).gt.0) then ! adiabatic turning off
                        kqist(3,k,jq)=1-kqist(3,k,jq)
                    elseif(kqist(3,k,jq).lt.0) then
                        kqist(3,k,jq)=1+kqist(3,k,jq)
                    endif
                endif
            endif
            ! END OF DEFINING CONTACT TYPES, NOW COMPUTING FORCES
            if(ldynss.and.inameseq(i).eq.4 .and. inameseq(j).eq.4
     +          .and.iconttype.eq.5 .and. knct3rd(i).eq.j) lss=.true.
            if(lss.or.(lcpot.and.(iconttype.gt.1))) then !attractive pot
                if(inameseq(i).eq.4 .and. inameseq(j).eq.4
     +          .and.iconttype.eq.5 .and. knct3rd(i).eq.j) then
                    if(kqist(3,k,jq).lt.0) then
                        totcoeff=-kqist(3,k,jq)/ad
                    else
                        totcoeff=1.d0
                    endif
                    rb = r - 6.d0/unit
                    rb2 = rb*rb
                    ene = totcoeff*(H1*rb2 + H2*rb2*rb2)
                    fce = totcoeff*((2*H1+4*H2*rb2)*rb)
                    lssn=(nei(2,i)+nei(2,j)).lt.neimaxdisul
                    if(lssn.and.rb2.gt.0.1) then ! DESTROYING SS BOND
c                      kqist(4,k,jq)=kremainder ! contact destroyed
c                      l3rdcn(i)=.false.
c                      l3rdcn(j)=.false.
c                      knct3rd(i)=0
c                      knct3rd(j)=0
                      if(kqist(4,k,jq).gt.0) then  ! are the same type
                          intrsc=intrsc-1
                      else
                          intesc=intesc-1
                      endif
                      if(kqist(3,k,jq).gt.0) then ! adiab. turning off
                          kqist(3,k,jq)=1-kqist(3,k,jq)
                      elseif(kqist(3,k,jq).lt.0) then
                          kqist(3,k,jq)=1+kqist(3,k,jq)
                      endif
                      if(ldynss) then
                        icheck=0
                        do kss=1,nssb
                          if(ksb(1,kss).eq.i.and.ksb(2,kss).eq.j) then
                            icheck=1
                            icnss=icnss-1
                          endif
                          if(ksb(2,kss).eq.i.and.ksb(1,kss).eq.j) then
                            icheck=1
                            icnss=icnss-1
                          endif
                        enddo
                        if(icheck.eq.0) icdss=icdss-1
                      endif
                    endif
                else
                    if(lsink.and.r.le.rsig*c216) then
                        if(r.le.cut) then
                            rsig=cut/c216 !*0.5d0**(1.d0/6.d0)
                            totcoeff=potcoeff ! repulsion always on
                        else
                            totcoeff=0.d0
                        endif
                    else
                if(lsink.and.(.not.lelectr).and.kqadabs.eq.0) goto 465
                      totcoeff=potcoeff*(kqadabs)/ad ! adiabatic
                    endif
                    if(iconttype.eq.4) totcoeff=2.d0*totcoeff
                    rsi=rsig/r
                    r6=rsi**6
                    ene=totcoeff*4.d0*r6*(r6-1.d0)
                    fce=totcoeff*24.d0*r6*(1.d0-2.d0*r6)/r
                    !r10=rsi**10.
                    !r12=r10*rsi*rsi
                    !ene=5.d0*r12-6.d0*r10
                    !fce=60.d0*(r10-r12)/r
                    if(r.lt.cut.and..not.lsink) then!repulsion always on
                        rsi=sigma0/r
                        r6=rsi**6
                        ene=ene+4.d0*r6*(r6-1.d0)+1.d0
                        fce=fce+ 24.d0*r6*(1.d0-2.d0*r6)/r
                    endif
                endif
            else
                if(r.gt.cut) then
                    if(.not.lelectr) goto 465
                    ene=0.d0
                    fce=0.d0
                else
                    rsi=sigma0/r
                    r6=rsi**6
                    ene=4.d0*r6*(r6-1.d0)+1.d0
                    fce= 24.d0*r6*(1.d0-2.d0*r6)/r
                endif
            endif
            
            ! ----------computing electrostatic forces------------------
            if(lcpot.and.lelectr) then
                if(itype.ne.jtype.and.iconttype.eq.5)then
                    
                    if(iconttype2.gt.1) then
                        adia(i)=min(adia(i)+1,ad) ! ELECTR. TURNED OFF
                        adia(j)=min(adia(j)+1,ad)
                    else
                        adia(i)=max(adia(i)-1,0.0) ! ELECTR. TURNED ON
                        adia(j)=max(adia(j)-1,0.0)
                    endif
                endif
!                if(lwritemap.and.iconttype.lt.3) then
!                    ncnt=ncnt+1
!                    nli(1,ncnt)=i
!                    nli(2,ncnt)=j
!                    nli(3,ncnt)=i*men+j !each cont. has uniq nr
!                    nli(4,ncnt)=1
!                    nli(5,ncnt)=1
!                    if(ldet) call compute_details(i,j)
!                endif
                adiab=max(adia(i),adia(j))
                if(lrepcoul) adiab=0.d0
                if(adiab.lt.ad) then
                    screenr=-1.d0*r/screend
                    expc=exp(screenr)
                    if(itype.eq.jtype) then
                        coulpotcoeff=expc*coul*(ad-adiab)/ad
                    else ! opposite signs
                        coulpotcoeff=expc*coul*(adiab-ad)/ad
                        if(lrepcoul) coulpotcoeff=0.d0
                    endif
                    if(lecperm) then
                        enec=coulpotcoeff/r
                        fce=fce+enec*(screenr-1.d0)/r
                    else
                        enec=coulpotcoeff/rsq
                        fce=fce+enec*(screenr-2.d0)/r
                    endif
                    ene=ene+enec ! ene and fce are initialised earlier
                    !fce=fce-2.d0*enec/r
                endif
            endif
            
            epot=epot+ene
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
465         continue
            if(kqist(3,k,jq).eq.0 .and. iconttype.gt.2) then
                kqist(4,k,jq)=kremainder ! contact destroyed
                if(iconttype.eq.5) then ! sidechain-sidechain
                    if(knct3rd(i).eq.j) then
                        l3rdcn(i)=.false.
                        l3rdcn(j)=.false.
                        knct3rd(i)=0
                        knct3rd(j)=0
                    endif
                    khbful(2,i)=khbful(2,i)-1
                    khbful(2,j)=khbful(2,j)-1
                    it2=2+min(itype,2)
                    jt2=2+min(jtype,2)
                    khbful(jt2,i)=khbful(jt2,i)-1
                    khbful(it2,j)=khbful(it2,j)-1
                else if(iconttype.eq.7) then !sidechain-backbone
                    khbful(2,i)=khbful(2,i)-1
                    khbful(1,j)=khbful(1,j)-1
                else if(iconttype.eq.6) then !backbone-sidechain
                    khbful(1,i)=khbful(1,i)-1
                    khbful(2,j)=khbful(2,j)-1
                else if(iconttype.eq.4) then ! backbone-backbone
                    khbful(1,i)=khbful(1,i)-1
                    khbful(1,j)=khbful(1,j)-1
                endif
            endif
466     continue
!$omp enddo nowait
!$omp end parallel
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE ENERGY AND FORCE OF THE PID MODEL

        subroutine evalimproper(epot,lcospid,epsbb)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lconftm,lwal,lwals,ldens,lcpot,lelectr,ljwal,lmj,lecperm
        logical lsimpang,lfrompdb(len),lconect(len),l3rdcn(len),lsim,lss
        logical lwritemap,lsink,lel,lmrs,ldet,lki,lkj,lpid,lsqpbc,ltr
        logical ldynss,lpbcx,lpbcy,lpbcz,lssn,lsselj,lbar,lcospid,lepid
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/restr/bckbmin,bckb2min,sdchnmax,sfact,kmaxbckb,lecperm
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
        dimension rij(2,3),rkj(2,3),rkl(2,3),bblbd(2),sslbd(2),dvdp(2)
        dimension fi(2,3),fj(2,3),fk(2,3),fl(2,3),rm(2,3),rn(2,3),a(2,4)
        
        do kk=1,kqont ! non-native attractive contacts
            i=kqist(1,kk,jq) ! kqist 1 or 2 = residue numbers
            j=kqist(2,kk,jq) ! jq = number of verlet list (1 or 2)
            if(lmj) then
                epsmj=emj(abs(kqist(4,kk,jq)))
            else
                epsmj=1.d0
            endif
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.lt.rneisq) then
                nei(1,i)=nei(1,i)+1
                nei(1,j)=nei(1,j)+1
            endif
            if(rsq.gt.rcutsq) goto 465
            if(inameseq(i).eq.2) goto 465
            if(inameseq(j).eq.2) goto 465
            do m=1,2 ! i is first, j is second
               bblbd(m)=0.d0
               sslbd(m)=0.d0
               dvdp(m)=0.d0
               psi=0.d0 
               i1=kqist(m,kk,jq)
               i2=kqist(m,kk,jq)+1
               i3=kqist(m,kk,jq)-1
               i4=kqist(3-m,kk,jq)
c               rij(m,1)=x0(i1)-x0(i2)
               rij(m,1)=-v(1,i1)
c               rij(m,2)=y0(i1)-y0(i2)
               rij(m,2)=-v(2,i1)
c               rij(m,3)=z0(i1)-z0(i2)
               rij(m,3)=-v(3,i1)
c               rkj(m,1)=x0(i3)-x0(i2)
               rkj(m,1)=v(4,i3)
c               rkj(m,2)=y0(i3)-y0(i2)
               rkj(m,2)=v(5,i3)
c               rkj(m,3)=z0(i3)-z0(i2)
               rkj(m,3)=v(6,i3)
               rkl(m,1)=x0(i3)-x0(i4)
               rkl(m,2)=y0(i3)-y0(i4)
               rkl(m,3)=z0(i3)-z0(i4)
c               rm(m,1)=rij(m,2)*rkj(m,3)-rij(m,3)*rkj(m,2)
               rm(m,1)=vxv(4,i1)
c               rm(m,2)=rij(m,3)*rkj(m,1)-rij(m,1)*rkj(m,3)
               rm(m,2)=vxv(5,i1)
c               rm(m,3)=rij(m,1)*rkj(m,2)-rij(m,2)*rkj(m,1)
               rm(m,3)=vxv(6,i1)
               rn(m,1)=rkj(m,2)*rkl(m,3)-rkj(m,3)*rkl(m,2)
               rn(m,2)=rkj(m,3)*rkl(m,1)-rkj(m,1)*rkl(m,3)
               rn(m,3)=rkj(m,1)*rkl(m,2)-rkj(m,2)*rkl(m,1)
               d2n=0.d0
               d2m=0.d0
               d2rkj=0.d0
               do k=1,3
                   d2n=d2n+rn(m,k)*rn(m,k)
                   d2m=d2m+rm(m,k)*rm(m,k)
                   d2rkj=d2rkj+rkj(m,k)*rkj(m,k)
               enddo
               if (abs(d2n).lt.0.01d0 .or. abs(d2m).lt.0.01d0) goto 465
               dmn=0.d0
               do k=1,3
                 dmn=dmn+rn(m,k)*rm(m,k)
               enddo
               cospsi=dmn/sqrt(d2n*d2m)
               cospsi=min(cospsi,1.d0)
               cospsi=max(cospsi,-1.d0)
               rijn=0.d0
               do k=1,3
                rijn=rijn+rij(m,k)*rn(m,k)
               enddo
               if (rijn.lt.0.d0) then
                psi=-1.d0*acos(cospsi)
               else
                psi=acos(cospsi)
               endif
               
               a(m,1)=alphacos(3)*(psi-psi0ss)
               a(m,3)=alphacos(3)
               if(a(m,1).lt.vmp .and.a(m,1).gt.-vmp) then
                if(lcospid) then
                 sslbd(m)=0.5*(cos(a(m,1))+1.d0)
                else
                 am1sq=a(m,1)*a(m,1)
                 sslbd(m)=1.d0-am1sq/(2.d0*am1sq-2*abs(a(m,1))+1.d0)
                endif
               endif
               if(j-i.ne.3 .or. m.ne.2) then
                   a(m,2)=alphacos(1)*(psi-psi0bb(1))
                   a(m,4)=alphacos(1)
                   rsig=rbb(1)
               else
                   a(m,2)=alphacos(2)*(psi-psi0bb(2))
                   a(m,4)=alphacos(2)
                   rsig=rbb(2)
               endif
               if(a(m,2).gt.vmp .or.a(m,2).lt.-vmp) then
                if(j-i.ne.3 .or. m.ne.1) then
                    a(m,2)=alphacos(2)*(psi-psi0bb(2))
                    a(m,4)=alphacos(2)
                    rsig=rbb(2)
                endif
               endif
               if(a(m,2).lt.vmp .and.a(m,2).gt.-vmp) then
                if(lcospid) then
                 bblbd(m)=0.5*(cos(a(m,2))+1.d0)
                else
                 am2sq=a(m,2)*a(m,2)
                 bblbd(m)=1.d0-am2sq/(2.d0*am2sq-2*abs(a(m,2))+1.d0)
                endif
               endif
               if(0.00005 .lt.sslbd(m).or.0.00005 .lt.bblbd(m)) then 
                   drkj=sqrt(d2rkj)
                   do k=1,3
                    fi(m,k)=rm(m,k)*drkj/d2m
                    fl(m,k)=-rn(m,k)*drkj/d2n
                   enddo
                   rijrkj=0.d0
                   rklrkj=0.d0
                   do k=1,3
                    rijrkj=rijrkj+rij(m,k)*rkj(m,k)
                    rklrkj=rklrkj+rkl(m,k)*rkj(m,k)
                   enddo
                   do k=1,3
                    df=(fi(m,k)*rijrkj-fl(m,k)*rklrkj)/d2rkj
                    fj(m,k)=-fi(m,k)+df
                    fk(m,k)=-fl(m,k)-df
                   enddo
               endif
            enddo
            sslmbd=sslbd(1)*sslbd(2)
            bblmbd=bblbd(1)*bblbd(2)
            if(bblmbd.lt.0.00005 .and. sslmbd.lt.0.00005) goto 465
            fce=0.d0
            r = sqrt(rsq)
            if(bblmbd.gt.0.00005) then
                if(r.lt.rsig*cntfct) then
                    ncord=ncord+2
                    if(kqist(4,kk,jq).gt.0) then    ! same chain
                        intrhc=intrhc+1
                        icnt(4)=icnt(4)+1                       
                    else
                        intehc=intehc+1
                        icnt(1)=icnt(1)+1
                    endif
                endif
                if(lsink.and.r.lt.rsig*c216) then
                    enelj=-1.d0*epsbb
                else
                    if(lbar) then
                      rsi=rsig*c216/r
                      r6=rsi**6
                      enelj=r6*(4.d0*r6-18.d0*rsi+13.d0)*epsbb
                fce=fce+6.d0*r6*(21.d0*rsi-8.d0*r6-13.d0)/r*bblmbd*epsbb
                    else
                      rsi=rsig/r
                      r6=rsi**6
                      enelj=4.d0*r6*(r6-1.d0)*epsbb
                      fce=fce+24.d0*r6*(1.d0-2.d0*r6)/r*bblmbd*epsbb
                    endif
                endif
                epot=epot+enelj*bblmbd
                do m=1,2
                 if(lcospid) then
                 dvdp(m)=dvdp(m)-0.5*a(m,4)*sin(a(m,2))*bblbd(3-m)*enelj
                 else
                  if(a(m,2).gt.0.d0) then
                    dgdx=2.d0*a(m,2)*(a(m,2)-1.d0)
        dvdp(m)=dvdp(m)+a(m,4)*(dgdx/(dgdx+1)**2)*bblbd(3-m)*enelj
c        write(*,*) "bb",i,j,m,psi,a(m,4)*(dgdx/(dgdx+1)**2),bblbd(m)
                  else
                    dgdx=-2.d0*a(m,2)*(a(m,2)+1.d0)
        dvdp(m)=dvdp(m)+a(m,4)*(dgdx/(dgdx-1)**2)*bblbd(3-m)*enelj
c        write(*,*) "bb",i,j,m,psi,a(m,4)*(dgdx/(dgdx-1)**2),bblbd(m)
                  endif
                 endif
                enddo
            endif
            if(epsmj.gt.0.00005 .and. sslmbd.gt.0.00005) then
                rsig=sigma1(abs(kqist(4,kk,jq)))
                if(r.lt.rsig*cntfct) then
                    ncord=ncord+2
                    if(kqist(4,kk,jq).gt.0) then    ! same chain
                        intrhc=intrhc+1
                        icnt(5)=icnt(5)+1                       
                    else
                        intehc=intehc+1
                        icnt(2)=icnt(2)+1
                    endif
                endif
                
                if(lepid.and.ksdchns(i).gt.2 .and.ksdchns(j).gt.2) then
                    screenr=-1.d0*r/screend
                    expc=exp(screenr)
                    if(ksdchns(i).eq.ksdchns(j)) then
                        coulpotcoeff=expc*coul
                    else ! opposite signs
                        coulpotcoeff=-expc*coul
                    endif
                    if(lecperm) then
                        enelj=coulpotcoeff/r
                        fce=fce+enelj*(screenr-1.d0)/r
                    else
                        enelj=coulpotcoeff/rsq
                        fce=fce+enelj*(screenr-2.d0)/r
                    endif
                elseif(lsink.and.r.lt.rsig*c216) then
                    enelj=-1.d0*epsmj
                elseif(lbar) then
                    rsi=rsig*c216/r
                    r6=rsi**6
                    enelj=r6*(4.d0*r6-18.d0*rsi+13.d0)*epsmj
                fce=fce+6.d0*r6*(21.d0*rsi-8.d0*r6-13.d0)/r*sslmbd*epsmj
                else
                    rsi=rsig/r
                    r6=rsi**6
                    enelj=4.d0*r6*(r6-1.d0)*epsmj
                    fce=fce+24.d0*r6*(1.d0-2.d0*r6)/r*sslmbd*epsmj
                endif
                epot=epot+enelj*sslmbd
                
                do m=1,2
                 if(lcospid) then
                 dvdp(m)=dvdp(m)-0.5*a(m,3)*sin(a(m,1))*sslbd(3-m)*enelj
                 else
                  if(a(m,1).gt.0.d0) then
                    dgdx=2.d0*a(m,1)*(a(m,1)-1.d0)
        dvdp(m)=dvdp(m)+a(m,3)*(dgdx/(dgdx+1)**2)*sslbd(3-m)*enelj
c        write(*,*) "ss",i,j,m,psi,a(m,3)*(dgdx/(dgdx+1)**2),sslbd(m)
                  else
                    dgdx=-2.d0*a(m,1)*(a(m,1)+1.d0)
        dvdp(m)=dvdp(m)+a(m,3)*(dgdx/(dgdx-1)**2)*sslbd(3-m)*enelj
c        write(*,*) "ss",i,j,m,psi,a(m,3)*(dgdx/(dgdx-1)**2),sslbd(m)
                  endif
                 endif
                enddo
            endif
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
            do m=1,2
                i1=kqist(m,kk,jq)
                i2=kqist(m,kk,jq)+1
                i3=kqist(m,kk,jq)-1
                i4=kqist(3-m,kk,jq)
                fx(i1)=fx(i1)-dvdp(m)*fi(m,1)
                fy(i1)=fy(i1)-dvdp(m)*fi(m,2)
                fz(i1)=fz(i1)-dvdp(m)*fi(m,3)
                fx(i2)=fx(i2)-dvdp(m)*fj(m,1)
                fy(i2)=fy(i2)-dvdp(m)*fj(m,2)
                fz(i2)=fz(i2)-dvdp(m)*fj(m,3)
                fx(i3)=fx(i3)-dvdp(m)*fk(m,1)
                fy(i3)=fy(i3)-dvdp(m)*fk(m,2)
                fz(i3)=fz(i3)-dvdp(m)*fk(m,3)
                fx(i4)=fx(i4)-dvdp(m)*fl(m,1)
                fy(i4)=fy(i4)-dvdp(m)*fl(m,2)
                fz(i4)=fz(i4)-dvdp(m)*fl(m,3)
            enddo
465     continue
        enddo
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE ENERGY AND FORCE OF THE GO MODEL

        subroutine evalgo(epot,chi)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lconftm,lwal,lwals,ldens,lcpot,lelectr,ljwal,lecperm,ltr
        logical lsimpang,lfrompdb(len),lconect(len),l3rdcn(len),lrepcoul
        logical lwritemap,lsink,lel,lmrs,ldet,lki,lkj,lpid,lsqpbc,lsim
        logical ldynss,lss,lpbcx,lpbcy,lpbcz,lssn,lsselj,lwalup,lwaldown
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/restr/bckbmin,bckb2min,sdchnmax,sfact,kmaxbckb,lecperm
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
c         chirality count does not include the gap
        
!$omp parallel default(shared), private(i,j,xi,yi,zi,dx,dy,dz,
!$omp& rsq,r,rb,rb2,ene,fce,repx,repy,repz,rsi,r6)
!$omp do reduction(+:epot)
        do 495 i=1,men-1
        if(lconect(i)) then ! lconect is true if i and i+1 are connected
            j=i+1   ! contacts between neighbour residues
            xi=x0(i)
            yi=y0(i)
            zi=z0(i)
            dx = xi-x0(j)
            dy = yi-y0(j)
            dz = zi-z0(j)
            rsq=dx*dx+dy*dy+dz*dz
            nei(1,i)=nei(1,i)+1
            nei(1,j)=nei(1,j)+1
            r = sqrt(rsq)
            rb = r - b(i)
            rb2 = rb*rb
            ene = H1*rb2 + H2*rb2*rb2
            fce = (2*H1+4*H2*rb2)*rb ! NEGATIVE FCE MEANS REPULSION
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
            epot=epot+ene
            if(lconect(j)) then ! i,i+2 contacts purely repulsive
                j=i+2
                dx = xi-x0(j)
                dy = yi-y0(j)
                dz = zi-z0(j)
                rsq=dx*dx+dy*dy+dz*dz
                if(rsq.lt.rneisq) then
                    nei(1,i)=nei(1,i)+1
                    nei(1,j)=nei(1,j)+1
                endif
                if(rsq.lt.cutsq) then
                    r = sqrt(rsq)
                    rsi=sigma0/r
                    r6=rsi**6
                    ene=4.d0*r6*(r6-1.d0)+1.d0
                    fce=24.d0*r6*(1.d0-2.d0*r6)/r
                    !r10=rsi**10.
                    !r12=r10*rsi*rsi
                    !ene=5.d0*r12-6.d0*r10
                    !fce=60.d0*(r10-r12)/r
                    if(fce.gt.1.d+3) fce=1.d+3
                    if(fce.lt.-1.d+3) fce=-1.d+3
                    fce=-fce/r
                    repx=fce*dx
                    repy=fce*dy
                    repz=fce*dz
                    fx(i) = fx(i) + repx
                    fx(j) = fx(j) - repx
                    fy(i) = fy(i) + repy
                    fy(j) = fy(j) - repy
                    fz(i) = fz(i) + repz
                    fz(j) = fz(j) - repz
                    epot=epot+ene
                endif
            endif
        endif
495     continue
!$omp enddo nowait
!$omp end parallel
!$omp parallel default(shared),
!$omp& private(i,j,k,dx,dy,dz,r,rsq,icm,rsi,r6,ene,fce,repx,repy,repz)
!$omp do
!$omp& reduction(+:epot)
        do 446 k=1,kront ! repulsive contacts
            i=krist(1,k)
            j=krist(2,k)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.lt.rneisq) then
                nei(1,i)=nei(1,i)+1
                nei(1,j)=nei(1,j)+1
            endif
            if(rsq.gt.cutsq) goto 446
            r = sqrt(rsq)
            rsi=sigma0/r
            r6=rsi**6
            ene=4.d0*r6*(r6-1.d0)+1.d0
            fce= 24.d0*r6*(1.d0-2.d0*r6)/r
            epot=epot+ene
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
446     continue
!$omp enddo nowait
!$omp end parallel
        do 496 k=1,kcont ! all other contacts (including native)
            i=kcist(1,k)
            j=kcist(2,k)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.lt.rneisq) then
                nei(1,i)=nei(1,i)+1
                nei(1,j)=nei(1,j)+1
            endif
            if(rsq.gt.rcutsq) goto 496
            r = sqrt(rsq)
            icm=kcist(3,k)
            if(icm.gt.0) then ! NATIVE CONTACTS
                if(abs(klist(3,icm)).lt.631) then ! L-J potential
                    rsig=sont(icm)
                    if(r.le.rsig*cntfct) then
                        icn=icn+1
                        ncord=ncord+2
                        !if(lwritemap) then
                        !    ncnt=ncnt+1
                        !    nli(1,ncnt)=i
                        !    nli(2,ncnt)=j
                        !    nli(3,ncnt)=abs(klist(3,icm))
                        !    nli(4,ncnt)=0
                        !    nli(5,ncnt)=sign(icm,klist(3,icm))
                        !endif
                        if(imap(icm).eq.0) then
                         imap(icm)=1
                         if(lcpot) then
                          if(khbful(2,i).lt.ksdchn(inameseq(i),2))then
                            khbful(2,i)=khbful(2,i)+1
                          endif
                          if(khbful(2,j).lt.ksdchn(inameseq(j),2))then
                            khbful(2,j)=khbful(2,j)+1
                          endif
                         endif
                        endif
                    else
                        if(imap(icm).eq.1) then
                         imap(icm)=0
                         if(lcpot) then
                          if(khbful(2,i).gt.0)then
                           khbful(2,i)=khbful(2,i)-1
                          endif
                          if(khbful(2,j).gt.0)then
                           khbful(2,j)=khbful(2,j)-1
                          endif
                         endif
                        endif
                    endif
                    rsi=rsig/r
                    r6=rsi**6
                    ene=4.d0*r6*(r6-1.d0)
                    fce= 24.d0*r6*(1.d0-2.d0*r6)/r
                else if(abs(klist(3,icm)).eq.631) then ! for SSbonds
c                   ene=ene*disul
c                   fce=fce*disul
                    icn=icn+1
                    rb = r - 6.d0/unit
                    rb2 = rb*rb
                    ene = H1*rb2 + H2*rb2*rb2
                    fce = (2*H1+4*H2*rb2)*rb
                endif
            else
                kremainder=mod(icm,2)
                icm=icm/2
                if(icm.eq.-1) then
                    if(r.gt.rcut) goto 496   ! rcut=20 A
                    if(lcpot.or.ldynss) then !EXCLUSIVE POT FOR SS BONDS
                        !if(lwritemap) then
                        !    ncnt=ncnt+1
                        !    nli(1,ncnt)=i
                        !    nli(2,ncnt)=j
                        !    nli(3,ncnt)=-1
                        !    nli(4,ncnt)=-1
                        !    if(kremainder.eq.0) then
                        !        nli(5,ncnt)=88
                        !    else
                        !        nli(5,ncnt)=-88
                        !    endif
                        !endif
                        if(l3rdcn(i).or.l3rdcn(j)) then
                            if(l3rdcn(i)) then
                                if(l3rdcn(j)) then
                                    adiab=min(adia(i),adia(j))
                                else
                                    adiab=adia(i)
                                endif
                            else
                                adiab=adia(j)
                            endif
                            if(knct3rd(i).eq.j.and.knct3rd(j).eq.i)then
                                !if(lwritemap) nli(4,ncnt)=3
                                totcoeff=ad
                                if(adiab.gt.0) then
                                    adia(i)=adia(i)-1
                                    adia(j)=adia(j)-1
                                endif
                                if(r.gt.sigss*cntfct) then !BREAK SSBOND
                                    l3rdcn(i)=.false.
                                    l3rdcn(j)=.false.
                                    adia(i)=0
                                    adia(j)=0
                                    if(kremainder.eq.0) then
                                        intesc=intesc-1
                                    else
                                        intrsc=intrsc-1
                                    endif
                                    icheck=0
                                    do kss=1,nssb
                          if(ksb(1,kss).eq.i.and.ksb(2,kss).eq.j) then
                                        icheck=1
                                        icnss=icnss-1
                          endif
                          if(ksb(2,kss).eq.i.and.ksb(1,kss).eq.j) then
                                        icheck=1
                                        icnss=icnss-1
                          endif
                                    enddo
                                    if(icheck.eq.0) icdss=icdss-1
                                endif
                            else ! OTHER CYS ARE NOT ATTRACTED TO SSBOND
                                !if(lwritemap) nli(4,ncnt)=2
                                if(adiab.gt.0) then
                                    totcoeff=adiab
                                else
                                    totcoeff=0.d0
                                endif
                            endif
                        else ! IF NO SSBOND IS FORMED, MANY CYS COMPETE
                            totcoeff=ad
                            if(r.le.sigss) then
                                l3rdcn(i)=.true.   ! ESTABLISH FULL
                                l3rdcn(j)=.true.   ! SSBOND FOR THE
                                knct3rd(i)=j         ! REMAINING TIME
                                knct3rd(j)=i         ! OF SIMULATION
                                adia(i)=ad
                                adia(j)=ad
                                if(kremainder.eq.0) then
                                    intesc=intesc+1
                                else ! INTESC/INTRSC ARE CHANGED ONCE
                                    intrsc=intrsc+1
                                endif
                                icheck=0
                                do kss=1,nssb
                          if(ksb(1,kss).eq.i.and.ksb(2,kss).eq.j) then
                                    icheck=1
                                    icnss=icnss+1
                          endif
                          if(ksb(2,kss).eq.i.and.ksb(1,kss).eq.j) then
                                    icheck=1
                                    icnss=icnss+1
                          endif
                                enddo
                                if(icheck.eq.0) icdss=icdss+1
                            endif
                        endif
                        if(lsselj) then
                         rsi=sigss/r
                         r6=rsi**6
                         if(r.lt.rmrs) then
                          ene=dislj*(4*r6*(r6-1.d0)+1)
                          fce=dislj*24*r6*(1.d0-2.d0*r6)/r
                         else
                          ene=dislj*totcoeff/ad*(4*r6*(r6-1.d0)+1)
                          fce=dislj*totcoeff/ad*24*r6*(1.d0-2.d0*r6)/r
                         endif
                        else
                         expm=exp(amrs*(rmrs-r))
                         ene=totcoeff*dmrs/ad*(1.d0-expm)**2
                         fce=2*totcoeff*amrs*dmrs*expm*(1.0-expm)/ad
                         if(r.lt.cut) then
                          rsi=sigma0/r
                          r6=rsi**6
                          ene=ene+(ad-totcoeff)*(4*r6*(r6-1.d0)+1)/ad
                         fce=fce+(ad-totcoeff)*24*r6/ad*(1.d0-2.d0*r6)/r
                         endif
                        endif
                    else
                        if(r.gt.cut) goto 496
                        rsi=sigma0/r
                        r6=rsi**6
                        ene=4.d0*r6*(r6-1.d0)+1.d0
                        fce= 24.d0*r6*(1.d0-2.d0*r6)/r
                    endif
                else if(icm.eq.-2) then
                    if(r.gt.rcut) goto 496      ! rcut=20 A
                    if(lcpot) then ! ELECTROSTATIC POTENTIAL
!                        if(lwritemap) then
!                            ncnt=ncnt+1
!                            nli(1,ncnt)=i
!                            nli(2,ncnt)=j
!                            nli(3,ncnt)=i*men+j!each cont. has uniq nr
!                            nli(4,ncnt)=1
!                     if(ldet) call compute_details(r,i,j,dx,dy,dz,0.d0)
!                        endif
                        coulpotcoeff=0.d0
                        adiab=max(adia(i),adia(j))
                        if(lrepcoul) adiab=0.d0
                        if(adiab.lt.ad) then
                            screenr=-1.d0*r/screend
                            expc=exp(screenr)
                            if(ksdchns(i).eq.ksdchns(j)) then
                                coulpotcoeff=expc*coul*(ad-adiab)/ad
                            else ! opposite signs
                                coulpotcoeff=expc*coul*(adiab-ad)/ad
                                if(lrepcoul) coulpotcoeff=0.d0
                            endif
                        endif
                        if(lecperm) then
                            ene=coulpotcoeff/r
                            fce=ene*(screenr-1.d0)/r
                        else
                            ene=coulpotcoeff/rsq
                            fce=ene*(screenr-2.d0)/r
                        endif
                        !fce=-2.d0*ene/r
                    endif
                    if(r.lt.cut) then
                        rsi=sigma0/r
                        r6=rsi**6
                        ene=ene+4.d0*r6*(r6-1.d0)+1.d0
                        fce=fce+24.d0*r6*(1.d0-2.d0*r6)/r
                    endif
                endif
            endif
            epot=epot+ene
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
496     continue

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  ASSIGN INITIAL VELOCITIES TO ATOMS

        subroutine intvel3d(aheat,part,nen1)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        sumx=0.
        sumy=0.
        sumz=0.
        do 200 i=nen1,men
        xx=2.d0*(ran2(0)-0.5d0)
        yy=2.d0*(ran2(0)-0.5d0)
        zz=2.d0*(ran2(0)-0.5d0)
        xyz=1.d0/ dsqrt(xx*xx+yy*yy+zz*zz)
        x1(i)=xx*xyz
        y1(i)=yy*xyz
        z1(i)=zz*xyz
        x2(i)=0.d0
        y2(i)=0.d0
        z2(i)=0.d0
        x3(i)=0.d0
        y3(i)=0.d0
        z3(i)=0.d0
        x4(i)=0.d0
        y4(i)=0.d0
        z4(i)=0.d0
        x5(i)=0.d0
        y5(i)=0.d0
        z5(i)=0.d0
        sumx=sumx+x1(i)
        sumy=sumy+y1(i)
        sumz=sumz+z1(i)
200     continue

        ! SCALE VELOCITIES SO THAT TOTAL MOMENTUM = ZERO
        x=0.d0
        do 210 i=nen1,men
        x1(i)=x1(i)-sumx/part
        y1(i)=y1(i)-sumy/part
        z1(i)=z1(i)-sumz/part
        x=x+x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i)
210     continue

        ! SCALE VELOCITIES TO DESIRED TEMPERATURE
        heat= dsqrt(aheat/x)
        do 220 i=nen1,men
        x1(i)=x1(i)*heat
        y1(i)=y1(i)*heat
        z1(i)=z1(i)*heat
220     continue
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  THE LANGEVIN NOISE

        subroutine lang(twopi,gamma2,const2,nen1)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/vel/x1(len),y1(len),z1(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        ! X-COMPONENT
        do 10 i=nen1,men
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         x1(i)=x1(i)+const2*gam
         fx(i)=fx(i)-gamma2*x1(i)
10      continue

        ! Y-COMPONENT
        do 20 i=nen1,men
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         y1(i)=y1(i)+const2*gam
         fy(i)=fy(i)-gamma2*y1(i)
20      continue

        ! Z-COMPONENT
        do 30 i=nen1,men
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         z1(i)=z1(i)+const2*gam
         fz(i)=fz(i)-gamma2*z1(i)
30      continue

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  THE LANGEVIN NOISE WHICH TAKES INTO ACCOUNT THE MASSES OF A.A.
        subroutine lang_mass(twopi,gamma2,const2,nen1)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/vel/x1(len),y1(len),z1(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/mass/rmas(len),rsqmas(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        do 10 i=nen1,men
         con2=const2/rsqmas(i)
         gam2=gamma2/rmas(i)
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         x1(i)=x1(i)+con2*gam
         fx(i)=fx(i)-gam2*x1(i)
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         y1(i)=y1(i)+con2*gam
         fy(i)=fy(i)-gam2*y1(i)
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         z1(i)=z1(i)+con2*gam
         fz(i)=fz(i)-gam2*z1(i)
10        continue

        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE MAKES FCC BEADS
        subroutine make_fcc()
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
c                            xup=xup+0.5*walmindst
c                            xdown=xdown-0.5*walmindst
c                            yup=yup+0.5*walmindst
c                            ydown=ydown-0.5*walmindst
            xaf=sqrt(3.0)/2.0*af        ! af is lattice constant
            ncxt=int((xup-xdown)/xaf) ! xaf is triangle height
            ncyt=int((yup-ydown)/af)
            do kwx=1,ncxt
                do kwy=1,ncyt
                    menw=menw+1
                    x0(men+menw)=xdown+(kwx-1)*xaf
                    y0(men+menw)=ydown+(kwy-1)*af+mod(kwx,2)*0.5*af
                    z0(men+menw)=zdown
                enddo
            enddo
            do kwx=1,ncxt
                do kwy=1,ncyt
                    menw=menw+1
                    x0(men+menw)=xdown+(kwx-1)*xaf
                    y0(men+menw)=ydown+(kwy-1)*af+mod(kwx-1,2)*0.5*af
                    z0(men+menw)=zup
                enddo
            enddo
            
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE CONNECTS BEADS TO THE WALL INSTANTENOUSLY
        subroutine connect_to_wal()
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
c                            xup=xup+0.5*walmindst
c                            xdown=xdown-0.5*walmindst
c                            yup=yup+0.5*walmindst
c                            ydown=ydown-0.5*walmindst
            
            do ib=1,men
                z0temp(ib)=z0(ib)
                ksorted(ib)=ib
            enddo
            call sort2(men,z0temp,ksorted) ! sort by z0
            j=1
            ip1=ipwn/2
            ip2=ipwn/2
            do k=1,ipwn/2
                i=ksorted(k)
                xpul(j)=x0(i)
                ypul(j)=y0(i)
                zpul(j)=z0(i)-zdown
c                aseq(i)='TYR'
                ipw(1,j)=-1
                ipw(2,j)=i
                j=j+1
                i=ksorted(men+1-k)
                xpul(j)=x0(i)
                ypul(j)=y0(i)
                zpul(j)=z0(i)-zup
c                aseq(i)='TYR'
                ipw(1,j)=1
                ipw(2,j)=i
                j=j+1
            enddo
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE CONNECTS BEADS TO THE WALL ONE BY ONE
        subroutine connect_to_wal_one_by_one()
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
            
            kstartconnected=ip1+1
            kendconnected=men-ip2
            do j=kstartconnected,kendconnected
                k=ksorted(j)
                if(z0(k).lt.zdown+walmindst) then
                    ip1=ip1+1
                    ksorted(j)=ksorted(ip1)
                    ksorted(ip1)=k
                    ipw(1,ip1)=-1
                    ipw(2,ip1)=k
                endif
                if(z0(k).gt.zup-walmindst) then
                    ip2=ip2+1
                    ksorted(j)=ksorted(ip2)
                    ksorted(ip2)=k
                    ipw(1,ip2)=1
                    ipw(2,ip2)=k
                endif
                if(z0(k).lt.zdown+walmindst
     +          .or.z0(k).gt.zup-walmindst) then
                    xpul(k)=x0(k)
                    ypul(k)=y0(k)
                    if(z0(k).lt.zdown+walmindst) zpul(k)=z0(k)-zdown
                    if(z0(k).gt.zup-walmindst) zpul(k)=z0(k)-zup
                    z0temp(k)=0
                    goto 9876
                endif
            enddo
9876    continue
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine displace(iprota,jprotb,away)
        implicit double precision(a-h,o-z)
        parameter(len=10000) !maximum number of all residues together
        logical lconect(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        dimension rcena(3),rcenb(3),qcena(3),qcenb(3),tow(3)

        do k=1,3
         rcena(k)=0
         rcenb(k)=0
         qcena(k)=0
         qcenb(k)=0
        enddo

        if(iprota.gt.nchains.or.jprotb.gt.nchains) then
            write(*,*) "INCORRECT NUMBERS OF CHAINS TO DISPLACE"
            stop
        endif
        
        mena = menchain(iprota+1)-menchain(iprota)
        menb = menchain(jprotb+1)-menchain(jprotb)
        
        do i=menchain(iprota)+1,menchain(iprota+1)
            rcena(1)=rcena(1)+x0(i)
            rcena(2)=rcena(2)+y0(i)
            rcena(3)=rcena(3)+z0(i)
            qcena(1)=qcena(1)+x0(i)*x0(i)
            qcena(2)=qcena(2)+y0(i)*y0(i)
            qcena(3)=qcena(3)+z0(i)*z0(i)
        enddo

        do i=menchain(jprotb)+1,menchain(jprotb+1)
            rcenb(1)=rcenb(1)+x0(i)
            rcenb(2)=rcenb(2)+y0(i)
            rcenb(3)=rcenb(3)+z0(i)
            qcenb(1)=qcenb(1)+x0(i)*x0(i)
            qcenb(2)=qcenb(2)+y0(i)*y0(i)
            qcenb(3)=qcenb(3)+z0(i)*z0(i)
        enddo
        
        do k=1,3
            rcena(k)=rcena(k)/mena
            qcena(k)=qcena(k)/mena
            qcena(k)=qcena(k)-rcena(k)*rcena(k)
            qcena(k)=sqrt(qcena(k))
            rcenb(k)=rcenb(k)/menb
            qcenb(k)=qcenb(k)/menb
            qcenb(k)=qcenb(k)-rcenb(k)*rcenb(k)
            qcenb(k)=sqrt(qcenb(k))
        enddo
        
        write(1,*)'#initial centers of mass and gyrations'
        write(1,*)'#A      xcen    ycen    zcen         xg    yg    zg'
        write(1,801) '#',(rcena(k)*unit,k=1,3),(qcena(l)*unit,l=1,3)
        write(1,*)'#B      xcen    ycen    zcen         xg    yg    zg'
        write(1,801) '#',(rcenb(k)*unit,k=1,3),(qcenb(l)*unit,l=1,3)
801     format(a,6x,3(f6.2,2x),4x,3(f6.2))
        
        toward=0
        do k=1,3
            tow(k)=rcenb(k)-rcena(k) ! direction from CM of A to CM of B
            toward=toward + (tow(k))**2
        enddo
        toward=dsqrt(toward)
        write(1,*)'CM separation between A and B ',sngl(toward*unit)

        do k=1,3
            tow(k)=tow(k)/toward
        enddo

        rxb=0.d0
        ryb=0.d0
        rzb=0.d0
        do i=menchain(jprotb)+1,menchain(jprotb+1)
c           bring the center of mass of B to the location of CM of A
            x0(i)=x0(i)+ rcena(1)-rcenb(1)
            y0(i)=y0(i)+ rcena(2)-rcenb(2)
            z0(i)=z0(i)+ rcena(3)-rcenb(3)
c           then move away
            x0(i)=x0(i)+tow(1)*away
            y0(i)=y0(i)+tow(2)*away
            z0(i)=z0(i)+tow(3)*away
            rxb=rxb+x0(i)
            ryb=ryb+y0(i)
            rzb=rzb+z0(i)
        enddo
        rxb=rxb/menb
        ryb=ryb/menb
        rzb=rzb/menb
        write(1,*)'shifted CM of B ',
     +  sngl(rxb*unit),sngl(ryb*unit),sngl(rzb*unit)
     
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE UPDATES VERLET LIST

        subroutine update_verlet_list(verlcut,nen1)
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        logical loverlap,lconect(len),lsamechain,lxov,lyov,lzov,lendchn
        logical lconftm,lcpot,lsim,lmrs,l3rdcn(len),lii4,lwal
        logical lpbcx,lpbcy,lpbcz,lpid,lfcc,lcintr,lepid
        character aseq*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        dimension pxmin(len),pymin(len),pzmin(len)
        dimension pxmax(len),pymax(len),pzmax(len)
        kront = 0
        kcont = 0
        kqont2 = kqont
        kqont = 0
        kfccw2=kfccw
        kfccw=0
        kactual = 1
        kqactual = 1
        kfactual = 1
        jq2=jq
        jq=3-jq ! if jq was 2, it is 1, if it was 1, it is 2
        cutoff=verlcut+rcut
        cutoff2=cutoff**2
        smallcutoff=verlcut+cut
        smcut2=smallcutoff**2
        
        do ib=1,men
            oxv(1,ib)=x0(ib)
            oxv(2,ib)=y0(ib)
            oxv(3,ib)=z0(ib)
            if(lfcc) then
              if(z0(ib)-zdown.lt.cutoff) then
                xi=x0(ib)
                yi=y0(ib)
                zi=z0(ib)
                do i=men+1,men+menw/2
                    dx = xi-x0(i)-shear
                    dy = yi-y0(i)
                    dz = zi-z0(i)
                    if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                    if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                    rsq=dx*dx+dy*dy+dz*dz
                    if(rsq.lt.cutoff2) then
                        kfccw=kfccw+1
                        kfcc(1,kfccw,jq)=ib
                        kfcc(2,kfccw,jq)=i
                        kfcc(3,kfccw,jq)=0
9129                    continue
                        if(kfactual.le.kfccw2) then
                         kf1=kfcc(1,kfactual,jq2)
                         kf2=kfcc(2,kfactual,jq2)
                         if(kf1.lt.ib .or.(kf1.eq.ib.and.kf2.lt.i))then
                            kfactual=kfactual+1
                            goto 9129
                         else if(kf1.eq.ib .and. kf2.eq.i) then
                            kfcc(3,kfccw,jq)=kfcc(3,kfactual,jq2)
                            kfactual=kfactual+1
                         endif
                        endif
                    endif
                enddo
              endif
              if(zup-z0(ib).lt.cutoff) then
                xi=x0(ib)
                yi=y0(ib)
                zi=z0(ib)
                do i=men+1+menw/2,men+menw
                    dx = xi-x0(i)+shear
                    dy = yi-y0(i)
                    dz = zi-z0(i)
                    if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                    if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                    rsq=dx*dx+dy*dy+dz*dz
                    if(rsq.lt.cutoff2) then
                        kfccw=kfccw+1
                        kfcc(1,kfccw,jq)=ib
                        kfcc(2,kfccw,jq)=i
                        kfcc(3,kfccw,jq)=0
9219                    continue
                        if(kfactual.le.kfccw2) then
                         kf1=kfcc(1,kfactual,jq2)
                         kf2=kfcc(2,kfactual,jq2)
                         if(kf1.lt.ib .or.(kf1.eq.ib.and.kf2.lt.i))then
                            kfactual=kfactual+1
                            goto 9219
                         else if(kf1.eq.ib .and. kf2.eq.i) then
                            kfcc(3,kfccw,jq)=kfcc(3,kfactual,jq2)
                            kfactual=kfactual+1
                         endif
                        endif
                    endif
                enddo
              endif
            endif
        enddo
        if((lpbcx.and.xsep.lt.0.001).or.(lpbcy.and.ysep.lt.0.001)
     +  .or.(lpbcz.and.zsep.lt.0.001)) then
            write(1,*) 'BROKEN PBC. PLS CHECK!' 
            stop
        endif
        
        ic=2
        do i=1,men
            if(i.gt.menchain(ic)) ic=ic+1
                
            if(lconect(i).and.lconect(i+1)) then
                kdist=i+3
            else if(lconect(i)) then
                kdist=i+2
            else
                kdist=i+1
            endif
            do j=kdist,men
                lsamechain=j.le.menchain(ic)
                lendchn=j.eq.i+1
                if(i.lt.nen1 .and. j.lt.nen1) goto 1129
                xi=x0(i)
                yi=y0(i)
                zi=z0(i)
                dx = xi-x0(j)
                dy = yi-y0(j)
                dz = zi-z0(j)
                if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                if(lpbcz) dz = dz-zsep*nint(dz*zinv)
                rsq=dx*dx+dy*dy+dz*dz
                if(rsq.lt.cutoff2) then
9797                continue
                    if(kactual.le.klont) then ! find native contacts
                        k1=klist(1,kactual)
                        k2=klist(2,kactual)
                        if(k1.lt.i .or. (k1.eq.i .and. k2.lt.j)) then
                            kactual=kactual+1
                            goto 9797
                        else if(k1.eq.i .and. k2.eq.j) then
                            kcont=kcont+1
                            kcist(1,kcont)=i
                            kcist(2,kcont)=j
                            kcist(3,kcont)=kactual
                            kactual=kactual+1
                            goto 1129
                        endif
                    endif
                    ! 3579 is the label for "the rest" of contacts
                    if(lendchn) goto 3579
                    iname1=inameseq(i)
                    iname2=inameseq(j)
                    if(lsamechain) then ! 3580 skips electrostatics
                        if(.not.lcintr) goto 3580
                        if(abs(j-i).eq.4 .and. .not. lii4) goto 3579
                    endif
                    if(iname1.eq.4 .and.iname2.eq.4 .and.lmrs) goto 3579
                    if(lpid) goto 9696 !no need 2 keep kqist for pid
9898                continue
                    if(kqactual.le.kqont2) then ! find non-native c.
                        kq1=kqist(1,kqactual,jq2)
                        kq2=kqist(2,kqactual,jq2)
                        if(kq1.lt.i .or.(kq1.eq.i.and.kq2.lt.j)) then
                            kqactual=kqactual+1
                            goto 9898
                        else if(kq1.eq.i .and. kq2.eq.j) then
                            kqont=kqont+1
                            kqist(1,kqont,jq)=i
                            kqist(2,kqont,jq)=j
                            kqist(3,kqont,jq)=kqist(3,kqactual,jq2)
                            kqist(4,kqont,jq)=kqist(4,kqactual,jq2)
                            kqactual=kqactual+1
                            !if(lpid) goto 3579
                            goto 1129
                        endif
                    endif
9696                continue
                    ! make a new non-native contact
                    kqont=kqont+1
                    kqist(1,kqont,jq)=i
                    kqist(2,kqont,jq)=j
                    kqist(3,kqont,jq)=0
                    if(lsamechain) then
                        kqist(4,kqont,jq)=1
                    else
                        kqist(4,kqont,jq)=-1
                    endif
                    if(lpid) then
                  kqist(4,kqont,jq)=kqist(4,kqont,jq)*(iname1*21+iname2)
                        if(lepid) goto 3580
                    else
                        goto 1129 ! TODO check if it works for lpid
                    endif
                        ! electrostatic, disulfide or repulsive contacts
3579                continue
                    if(.not.lcpot) goto 3581 !skip ssbond and electr
                    if(ksdchns(i).gt.2 .and.ksdchns(j).gt.2) then
                        kcont=kcont+1
                        kcist(1,kcont)=i
                        kcist(2,kcont)=j
                        kcist(3,kcont)=-4 ! electrostatic interaction
                        if(lsamechain) kcist(3,kcont)=-5
                        goto 1129
                    endif
3580                continue
                    iname1=inameseq(i)
                    iname2=inameseq(j)
                    if(iname1.eq.4 .and.iname2.eq.4 .and.lmrs) then
                        kcont=kcont+1
                        kcist(1,kcont)=i
                        kcist(2,kcont)=j
                        kcist(3,kcont)=-2 ! nonnative SSbond
                        if(lsamechain) kcist(3,kcont)=-3
                        goto 1129
                    endif
3581                continue
                    if(rsq.lt.smcut2) then ! REPULSIVE CONTACT
                        kront=kront+1
                        krist(1,kront)=i
                        krist(2,kront)=j
                        krist(3,kront)=0
                    endif
                endif
1129            continue
            enddo
        enddo
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE GENERATE THE STARTING CONFIGURATION

        subroutine confstart(sdens,confcut)
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        logical lwal,lwals,ldens,lsim,lconect(len),lcpb
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
        dimension phi(len),theta(len)
        dimension T(len,3,3),R(len,3)
        
        ! this is "infinite" box, for close-packed one insert zeroes
        if(ldens) then ! starting density is 0.0001 residues/A^3
            startvolume=(men)/(sdens) !*unit**3)
            startboxsize=0.5*startvolume**(1.0/3.0)
            xmin=-startboxsize
            ymin=-startboxsize
            zmin=-startboxsize
            xmax=startboxsize
            ymax=startboxsize
            zmax=startboxsize
            if(lcpb) then
                xdown=xmin
                xup=xmax
                ydown=ymin
                yup=ymax
                zdown=zmin
                zup=zmax
                xsep=xup-xdown
                ysep=yup-ydown
                zsep=zup-zdown
                xinv=1.d0/xsep
                yinv=1.d0/ysep
                zinv=1.d0/zsep
            endif
        else
            if(lcpb) then
                write(2,*)'lcpb must be used together with ldens'
                stop
            endif
            xmin=0.d0
            ymin=0.d0
            zmin=0.d0
            xmax=0.d0
            ymax=0.d0
            zmax=0.d0
        endif
        
        do 1000 ic=1,nchains
        kter=0
100     continue
        ! ASSIGN RANDOM VALUE FOR PHI AND THETA
        phi(menchain(ic)+1)=Pi/2
        theta(menchain(ic)+1)=0.d0
        phi(menchain(ic)+2)=0.d0
        theta(menchain(ic)+2)=ran2(0)*Pi/3
        do i=menchain(ic)+3,menchain(ic+1)-1
        phi(i)=(2.d0*ran2(0)-1.d0)*Pi
        theta(i)=ran2(0)*Pi/3        ! theta <= Pi/3
        enddo

        ! COMPUTE THE TRANSFER MATRICES
        do ib=menchain(ic)+1,menchain(ic+1)-1
        T(ib,1,1)=cos(theta(ib))
        T(ib,1,2)=sin(theta(ib))
        T(ib,1,3)=0.d0
        T(ib,2,1)=sin(theta(ib))*cos(phi(ib))
        T(ib,2,2)=-cos(theta(ib))*cos(phi(ib))
        T(ib,2,3)=sin(phi(ib))
        T(ib,3,1)=sin(theta(ib))*sin(phi(ib))
        T(ib,3,2)=-cos(theta(ib))*sin(phi(ib))
        T(ib,3,3)=-cos(phi(ib))
        enddo
        ! COMPUTE RANDOM DIRECTION OF PROTEIN CHAIN
        theta0=dacos(1.d0-2.d0*ran2(0))
        phi0=2.d0*Pi*ran2(0)
        ! COMPUTE THE BACK-BONE VECTORS
        do ib=menchain(ic)+1,menchain(ic+1)-1
        r1=bond*sin(theta0)*cos(phi0)
        r2=bond*sin(theta0)*sin(phi0)
        r3=bond*cos(theta0)
        do i=ib,menchain(ic)+1,-1
        R(ib,1)=T(i,1,1)*r1+T(i,1,2)*r2+T(i,1,3)*r3
        R(ib,2)=T(i,2,1)*r1+T(i,2,2)*r2+T(i,2,3)*r3
        R(ib,3)=T(i,3,1)*r1+T(i,3,2)*r2+T(i,3,3)*r3
        r1=R(ib,1)
        r2=R(ib,2)
        r3=R(ib,3)
        enddo
        enddo
        ranx=(xmax-xmin)*ran2(0)+xmin
        rany=(ymax-ymin)*ran2(0)+ymin
        ranz=(zmax-zmin)*ran2(0)+zmin
        ! COMPUTE THE POSITIONS OF MONOMERS
        do ib=menchain(ic)+1,menchain(ic+1)
        x0(ib)=ranx
        y0(ib)=rany
        z0(ib)=ranz
        do i=menchain(ic)+1,ib-1
        x0(ib)=x0(ib)+R(i,1)
        y0(ib)=y0(ib)+R(i,2)
        z0(ib)=z0(ib)+R(i,3)
        enddo
        enddo
        ! CHECK THE CONFORMATION
        do i=1,menchain(ic+1)-3
        do j=i+3,menchain(ic+1)
        dx = x0(j) - x0(i)
        dy = y0(j) - y0(i)
        dz = z0(j) - z0(i)
        if(lcpb) then
            dx = dx-xsep*nint(dx*xinv)
            dy = dy-ysep*nint(dy*yinv)
            dz = dz-zsep*nint(dz*zinv)
        endif
        rs2 = dx*dx+dy*dy+dz*dz
        if(rs2.lt.confcut*confcut) then
          kter = kter+1
          if(kter.lt.9000) then
            goto 100
          else
            write(2,*)'Confstart: FAILURE'
          stop
          endif
        endif
        enddo
        enddo
        if(.not.lcpb) then
            do ib=1,menchain(ic+1)
            if(x0(ib).lt.xmin) xmin=x0(ib)
            if(y0(ib).lt.ymin) ymin=y0(ib)
            if(z0(ib).lt.zmin) zmin=z0(ib)
            if(x0(ib).gt.xmax) xmax=x0(ib)
            if(y0(ib).gt.ymax) ymax=y0(ib)
            if(z0(ib).gt.zmax) zmax=z0(ib)
            enddo
        endif
1000    continue
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE COMPUTES THE RADIUS OF GYRATION

        subroutine gyration(rg)
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/masscenter/xcm,ycm,zcm,xmcm,ymcm,zmcm
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/mass/rmas(len),rsqmas(len)

        xcm=0.d0
        ycm=0.d0
        zcm=0.d0
        xmcm=0.d0
        ymcm=0.d0
        zmcm=0.d0
        do i=1,men
        xcm=xcm+x0(i)
        ycm=ycm+y0(i)
        zcm=zcm+z0(i)
        xmcm=xmcm+x0(i)*rmas(i)
        ymcm=ymcm+y0(i)*rmas(i)
        zmcm=zmcm+z0(i)*rmas(i)
        enddo
        xcm=xcm/men
        ycm=ycm/men
        zcm=zcm/men
        xmcm=xmcm/men-xdown
        ymcm=ymcm/men-ydown
        zmcm=zmcm/men-zdown
        rg=0.d0
        do i=1,men
        dx=x0(i)-xcm
        dy=y0(i)-ycm
        dz=z0(i)-zcm
        rg=rg+dx*dx+dy*dy+dz*dz
        enddo
        rg=rg/men
        rg=sqrt(rg)
        return
        end

C THIS SUBROUTINE COMPUTES W AND RADIUS OF GYRATION OF INDIVIDUAL CHAINS        
        subroutine cgyration()
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        
        logical lconect(len),lkmt
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/mass/rmas(len),rsqmas(len)
        common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        double precision mit(3,3),vmit(3,3),ra(3)

        do ic=1,nchains
        klewy=0
        kprawy=0
        if(lkmt) call kmt(klewy,kprawy,menchain(ic)+1,menchain(ic+1))
        chainlength=menchain(ic+1)-menchain(ic)
        if(kprawy-klewy.eq.chainlength) then
            knts(1,ic) = 0 ! no knots
            knts(2,ic) = 0 ! no knots
        else
            knts(1,ic)=klewy  ! left end of the knot
            knts(2,ic)=kprawy ! right end of the knot
        endif
        cxcm=0.d0
        cycm=0.d0
        czcm=0.d0
        do i=1,3
          do j=1,3
            mit(i,j)=0.d0
          enddo
        enddo
        do i=menchain(ic)+1,menchain(ic+1)
        cxcm=cxcm+x0(i)
        cycm=cycm+y0(i)
        czcm=czcm+z0(i)
        enddo
        cxcm=cxcm/chainlength
        cycm=cycm/chainlength
        czcm=czcm/chainlength
        xyzcm(1,ic)=cxcm
        xyzcm(2,ic)=cycm
        xyzcm(3,ic)=czcm
        crg=0.d0
        do i=menchain(ic)+1,menchain(ic+1)
        dx=x0(i)-cxcm
        dy=y0(i)-cycm
        dz=z0(i)-czcm
        crg=crg+dx*dx+dy*dy+dz*dz
        mit(1,1)=mit(1,1)+dy*dy+dz*dz
        mit(2,2)=mit(2,2)+dz*dz+dx*dx
        mit(3,3)=mit(3,3)+dx*dx+dy*dy
        mit(1,2)=mit(1,2)-dx*dy
        mit(2,3)=mit(2,3)-dy*dz
        mit(3,1)=mit(3,1)-dz*dx
        enddo
        mit(2,1)=mit(1,2)
        mit(3,2)=mit(2,3)
        mit(1,3)=mit(3,1)
        crg=crg/chainlength
        cgyr(ic)=sqrt(crg)
        call jacobi(mit,3,3,ra,vmit,nrot)
        do i=1,3
         ra(i)=sqrt(ra(i)/chainlength)
        enddo
        call sort(3,ra)
        rb=(ra(1)+ra(3))*0.5d0
        dr=ra(2)-rb
        w(ic)=dr/rb
        dx=x0(menchain(ic)+1)-x0(menchain(ic+1))
        dy=y0(menchain(ic)+1)-y0(menchain(ic+1))
        dz=z0(menchain(ic)+1)-z0(menchain(ic+1))
        cend(ic)=sqrt(dx*dx+dy*dy+dz*dz) !/chainlength
        enddo
        
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES DETAILS OF INTER-RESIDUE ANGLES

        subroutine compute_details(i,j)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lpbcx,lpbcy,lpbcz
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        
        xi=x0(i)
        yi=y0(i)
        zi=z0(i)
        dx = xi-x0(j)
        dy = yi-y0(j)
        dz = zi-z0(j)
        if(lpbcx) dx = dx-xsep*nint(dx*xinv)
        if(lpbcy) dy = dy-ysep*nint(dy*yinv)
        if(lpbcz) dz = dz-zsep*nint(dz*zinv)
        rsq=dx*dx+dy*dy+dz*dz
        r = sqrt(rsq)
        rangli(1)=r*unit
        bbir=(vxv(1,i)*dx+vxv(2,i)*dy+vxv(3,i)*dz)/r
        bbjr=(vxv(1,j)*dx+vxv(2,j)*dy+vxv(3,j)*dz)/r
        bbij=vxv(1,i)*vxv(1,j)+vxv(2,i)*vxv(2,j)+vxv(3,i)*vxv(3,j)

        i0=i-1
        sdir=(v(1,i0)-v(1,i))*dx+(v(2,i0)-v(2,i))*dy+(v(3,i0)-v(3,i))*dz
        sdirnrm=0.0
        if(sdir.lt.0.0) then ! if <0, no need to compute norm
            sdchni=0.0
        else
            do l=1,3
                sdirnrm=sdirnrm+(v(l,i)-v(l,i0))**2
            enddo
            sdir=sdir**2/(sdirnrm*r*r) ! cos^2(n,r)
            sdchni=sdir!min(1.0,max(sdchnmax-sdir,0.0)/sdchnmin)
        endif

        j0=j-1
        sdjr=(v(1,j)-v(1,j0))*dx+(v(2,j)-v(2,j0))*dy+(v(3,j)-v(3,j0))*dz
        sdjrnrm=0.0
        if(sdjr.lt.0.0) then
            sdchnj=0.0
        else
            do l=1,3
                sdjrnrm=sdjrnrm+(v(l,j)-v(l,j0))**2
            enddo
            sdjr=sdjr**2/(sdjrnrm*r*r)
            sdchnj=sdjr!min(1.0,max(sdchnmax-sdjr,0.0)/sdchnmin)
        endif

        rangli(2)=abs(bbir)
        rangli(3)=abs(bbjr)
        rangli(4)=abs(bbij)
        rangli(5)=sdchni
        rangli(6)=sdchnj
        phi=0.d0
        call countpid(i,j,phi)
        rangli(7)=phi
        phi=0.d0
        call countpid(j,i,phi)
        rangli(8)=phi
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES SS-BONDS WITHIN GO-MODEL
        subroutine compute_ssbonds()
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical l3rdcn(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        dimension rsortssbond(len),ksortssbond(len)

        do i=1,men            ! ZERO THE 3RD ORDER CONNECTION TABLE
        l3rdcn(i)=.false. ! true if residue i forms a disulfide bond
        knct3rd(i)=0        ! index of residue that is bonded to i
        enddo
        l=0
        do 496 k=1,kront ! all other contacts
            i=krist(1,k)
            j=krist(2,k)
            icm=krist(3,k)
            xi=x0(i)
            yi=y0(i)
            zi=z0(i)
            dx = xi-x0(j)
            dy = yi-y0(j)
            dz = zi-z0(j)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.gt.rcutsq) goto 496
            r = sqrt(rsq)
            if((icm.eq.-1).and.(r.le.rmrs*1.5d0)) then
                l=l+1
                ksortssbond(l)=k
                rsortssbond(l)=r
            endif
496     continue
        if(l.gt.1) then
            call sort2(l,rsortssbond,ksortssbond)
        endif    ! ALL NON-NATIVE POSSIBLE SSBONDS ARE NOW SORTED BY R
        do m=1,l ! CLOSEST CYSTEINS ARE CONNECTED FIRST
            k=ksortssbond(m)
            i=krist(1,k)
            j=krist(2,k)
            if((.not.l3rdcn(i)).and.(.not.l3rdcn(j))) then
                icm=krist(3,k)
                kremainder=mod(icm,2)
                l3rdcn(i)=.true.
                l3rdcn(j)=.true.
                knct3rd(i)=j
                knct3rd(j)=i 
                if(kremainder.eq.0) then
                    intesc=intesc+1
                else
                    intrsc=intrsc+1
                endif
            endif
        enddo
        
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE NATIVE COUPLINGS WITHIN GO-MODEL
c   based on cut-off length dnat
        subroutine compute_cmap(dnat,lcdnat)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lpbcx,lpbcy,lpbcz,lconect(len),lcdnat
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        character aseq*3,ares*3,filn*32,bb*2,buffer*128,ch1*1,ch2*1,ch*1
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
        klont=0

        do 10000 ic=1,nchains
            do 10000 icc=ic,nchains
                do 2000 ib1=menchain(ic)+1,menchain(ic+1)
                    if(ic.eq.icc) then
                        kdist=ib1+3
                    else
                        kdist=menchain(icc)+1
                    endif
                    do 2000 ib2=kdist,menchain(icc+1)
                        dx=xn(ib1)-xn(ib2)
                        dy=yn(ib1)-yn(ib2)
                        dz=zn(ib1)-zn(ib2)
                        if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                        if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                        if(lpbcz) dz = dz-zsep*nint(dz*zinv)
                        dal=dx*dx+dy*dy+dz*dz
                        dal=dsqrt(dal)
                        dcut=dnat
             if(lcdnat) dcut=c216*sigma1(inameseq(ib1)*21+inameseq(ib2))
                        if(dal.le.dcut) then
                          klont=klont+1
                          klist(1,klont)=ib1
                          klist(2,klont)=ib2
                          if(ic.eq.icc) then
                            klist(3,klont)=1
                          else
                            klist(3,klont)=-1
                          endif
                        endif
2000            continue
10000   continue

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE READS THE NATIVE COUPLINGS WITHIN GO-MODEL FROM FILE
        subroutine load_cmap(cmapf)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character cmapf*32
        logical lfrompdb(len)
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        
        open(98,file=cmapf,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,2x,a)')'ERROR OPENING FILENAME',cmapf
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        read(98,*) icn
        read(98,*) men
        do i=1,icn
          read(98,*)(klist(j,i),j=1,2),sont(i)
          klist(3,i)=1 ! TODO make it -1 for different chains
        enddo

        klont=icn
        do i=1,men
            read(98,*) the0(i),phi0(i)
            lfrompdb(i)=.true.
        enddo
        close(98)
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE NATIVE COUPLINGS WITHIN GO-MODEL
c   based on all-atom VdW spheres covering
        subroutine compute_contact_map(filename)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,filename*32,aname*3
        logical lconect(len)
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/rall/rx_a(len,14),ry_a(len,14),rz_a(len,14)
        common/nall/na(len),aname(len,14)
        common/radi/vrad(len,14)
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        alpha=(26/7.)**(1./6.)

        call load_allatom(filename)
        call assign_VdW_radius

        klont=0

        ! BUILD CONTACT MAP
        do 10000 ic=1,nchains
            do 10000 icc=ic,nchains
c               do 2000 ib1=1,men-2
                do 2000 ib1=menchain(ic)+1,menchain(ic+1)
                    na1=na(ib1)
                    if(ic.eq.icc) then
                        kdist=ib1+3
                    else
                        kdist=menchain(icc)+1
                    endif
c                   do 2000 ib2=ib1+2,men
                    do 2000 ib2=kdist,menchain(icc+1)
                        na2=na(ib2)
                        rmin=1.E+6
                        kcc=0
                        kccbb=0
                        kccbs=0
                        kccss=0
                        kccbsbs=0
                        kccbssb=0
                        do 200 ja1=1,na1
                        do 200 ja2=1,na2
                            rx=rx_a(ib1,ja1)-rx_a(ib2,ja2)
                            ry=ry_a(ib1,ja1)-ry_a(ib2,ja2)
                            rz=rz_a(ib1,ja1)-rz_a(ib2,ja2)
                            rr=sqrt(rx*rx+ry*ry+rz*rz)
                            vrsum=vrad(ib1,ja1)+vrad(ib2,ja2)
                            if(rr.le.vrsum*alpha) then
                                if(ja1.lt.5.and.ja2.lt.5) then
                                    kccbb=kccbb+1
                                else if(ja1.gt.4.and.ja2.gt.4) then
                                    kccss=kccss+1
                                else
                                    kccbs=kccbs+1
                                    if(ja1.lt.5.and.ja2.gt.4) then
                                        kccbsbs=kccbsbs+1
                                    endif
                                    if(ja1.gt.4.and.ja2.lt.5) then
                                        kccbssb=kccbssb+1
                                    endif
                                endif
                                kcc=kcc+1
                            endif
                            if(rr.lt.rmin) rmin=rr
200                     continue
                        if(kcc.gt.0) then
                        klont=klont+1
                        klist(1,klont)=ib1
                        klist(2,klont)=ib2
                        if(ic.eq.icc) then
                            klist(3,klont)=1
                        else
                            klist(3,klont)=-1
                        endif
                        if(kccbb.gt.0) klist(3,klont)=klist(3,klont)*2
                        if(kccbs.gt.0) klist(3,klont)=klist(3,klont)*3
                        if(kccss.gt.0) klist(3,klont)=klist(3,klont)*5
                        if(kccbsbs.gt.0) klist(3,klont)=klist(3,klont)*3
                        if(kccbssb.gt.0) klist(3,klont)=klist(3,klont)*7
                        endif
2000            continue
10000   continue
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ran2(iseed)
      implicit double precision(a-h,o-z)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
c      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      common/rans/idum
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if(iseed.ne.0) idum=-iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,)B.n.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE LOADS THE NATIVE BACKBONE CONFORMATION OF A PROTEIN
C FROM ITS PDB FILE
        subroutine load_protein(filn,lunwrap)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lfrompdb(len),lsimpang,lpdb,lcoilang,ldynss
        logical l3rdcn(len),lunwrap,lpbcx,lpbcy,lpbcz,lcheck
        character aseq*3,ares*3,filn*32,bb*2,buffer*128,ch1*1,ch2*1,ch*1
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        dimension ch(len)
        character bb2*2,bb4*4
        open(8,file=filn,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,2x,a)')'ERROR OPENING FILENAME',filn
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        
        do ib=1,len
            lconect(ib)=.false.
        enddo
        nchains=0
        menchain(1)=0
        ib=0
        xdown=0.d0
        ydown=0.d0
        zdown=0.d0
15      read(8,'(a)',end=20) buffer
        if(buffer(1:4).ne.'ATOM') then
            if(buffer(1:3).eq.'END') goto 20
            if(buffer(1:3).eq.'TER') then
                if(ib.gt.menchain(nchains+1)) then !filter out DNA chain
                    lconect(ib)=.false.
                    nchains=nchains+1
                    menchain(nchains+1)=ib
                endif
            elseif(buffer(1:6).eq.'CRYST1') then
                read(buffer(1:33),'(6x,3f9.3)') xup,yup,zup
                xup=xup/unit
                yup=yup/unit
                zup=zup/unit
                xsep=xup-xdown
                ysep=yup-ydown
                zsep=zup-zdown
                xinv=1.d0/xsep
                yinv=1.d0/ysep
                zinv=1.d0/zsep
            endif
            goto 15
        endif
        if((buffer(17:17).ne.'A').and.(buffer(17:17).ne.' ')) goto 15
        read(buffer,'(13x,a4,a3,x,a1,i4,4x,3f8.3)')
     +  bb4,ares,ch1,ival,xval,yval,zval
        read(bb4,'(a2,2x)') bb
        read(bb4,'(x,a2,x)') bb2
        if(bb.eq.'CA'.or.bb2.eq.'CA') then
        ib=ib+1
        lconect(ib)=.true.
        if(lpdb.and..not.lcoilang) lfrompdb(ib)=.true.
        aseq(ib)=ares
        iseq(ib)=ival
        ch(ib)=ch1
        if(ib.gt.2) then
            if(lconect(ib-1).and.ch(ib-1).ne.ch(ib)) then
                lconect(ib-1)=.false.
                nchains=nchains+1
                menchain(nchains+1)=ib-1
            endif
        endif
        xn(ib)=xval/unit
        yn(ib)=yval/unit
        zn(ib)=zval/unit
        if(aseq(ib).eq.'GLY') then
        inameseq(ib)=1
        else if(aseq(ib).eq.'PRO') then
        inameseq(ib)=2
        else if(aseq(ib).eq.'GLN') then
        inameseq(ib)=3
        else if(aseq(ib).eq.'CYS') then
        inameseq(ib)=4
        else if(aseq(ib).eq.'ALA') then
        inameseq(ib)=5
        else if(aseq(ib).eq.'SER') then
        inameseq(ib)=6
        else if(aseq(ib).eq.'VAL') then
        inameseq(ib)=7
        else if(aseq(ib).eq.'THR') then
        inameseq(ib)=8
        else if(aseq(ib).eq.'ILE') then
        inameseq(ib)=9
        else if(aseq(ib).eq.'LEU') then
        inameseq(ib)=10
        else if(aseq(ib).eq.'ASN') then
        inameseq(ib)=11
        else if(aseq(ib).eq.'ASP') then
        inameseq(ib)=12
        else if(aseq(ib).eq.'LYS') then
        inameseq(ib)=13
        else if(aseq(ib).eq.'GLU') then
        inameseq(ib)=14
        else if(aseq(ib).eq.'MET') then
        inameseq(ib)=15
        else if(aseq(ib).eq.'HIS') then
        inameseq(ib)=16
        else if(aseq(ib).eq.'PHE') then
        inameseq(ib)=17
        else if(aseq(ib).eq.'ARG') then
        inameseq(ib)=18
        else if(aseq(ib).eq.'TYR') then
        inameseq(ib)=19
        else if(aseq(ib).eq.'TRP') then
        inameseq(ib)=20
        endif
        endif
        goto 15
20      continue
        close(8)
        men=ib
        
c        if(.not.ldynss) then
        ! SSBONDS
        !nssb=0
        open(8,file=filn,status='old',iostat=ierr)
16      read(8,'(a)',end=21) buffer
        if(buffer(1:6).eq.'SSBOND') then
            read (buffer(16:21), '(a1,i5)' ) ch1,icys1
            read (buffer(30:35), '(a1,i5)' ) ch2,icys2
            nssb=nssb+1
            do j=1,men
             if (ch(j).eq.ch1 .and. iseq(j).eq.icys1) then
              ksb(1,nssb)=j
             endif
             if (ch(j).eq.ch2 .and. iseq(j).eq.icys2) then
              ksb(2,nssb)=j
             endif
            enddo
        endif
        if(buffer(1:3).eq.'END') goto 21
        goto 16
21      continue
        close(8)
c        endif
        if(men.eq.0) then
            write(*,'(a,a)')'NO ATOMS IN THE FILE ',filn
            write(*,*)'PROGRAM STOPPED.'
            stop
        endif
        if(lconect(men)) then ! if pdb file does not have TER or CH tags
            lconect(men)=.false.
            nchains=nchains+1
            menchain(nchains+1)=men
        endif
        if(lunwrap) then
            xcentrall=0.d0
            ycentrall=0.d0
            zcentrall=0.d0
            xicent=0.d0
            yicent=0.d0
            zicent=0.d0
            xjcent=0.d0
            yjcent=0.d0
            zjcent=0.d0
            Pi=dacos(-1.d0)! PI
            pi2=2.d0*Pi    ! 2PI
            do i=1,men
                xicent=xicent+xsep/pi2*cos(xn(i)/xsep*pi2)
                yicent=yicent+ysep/pi2*cos(yn(i)/ysep*pi2)
                zicent=zicent+zsep/pi2*cos(zn(i)/zsep*pi2)
                xjcent=xjcent+xsep/pi2*sin(xn(i)/xsep*pi2)
                yjcent=yjcent+ysep/pi2*sin(yn(i)/ysep*pi2)
                zjcent=zjcent+zsep/pi2*sin(zn(i)/zsep*pi2)
            enddo
            xicent=xicent/men
            yicent=yicent/men
            zicent=zicent/men
            xjcent=xjcent/men
            yjcent=yjcent/men
            zjcent=zjcent/men
            if(lpbcx) xcentrall=xsep/pi2*(atan2(-xjcent,-xicent)+Pi)
            if(lpbcy) ycentrall=ysep/pi2*(atan2(-yjcent,-yicent)+Pi)
            if(lpbcz) zcentrall=zsep/pi2*(atan2(-zjcent,-zicent)+Pi)
            
            do ic=1,nchains
                xcentr=0.d0
                ycentr=0.d0
                zcentr=0.d0
                xicent=0.d0
                yicent=0.d0
                zicent=0.d0
                xjcent=0.d0
                yjcent=0.d0
                zjcent=0.d0
                do i=menchain(ic)+1,menchain(ic+1)
                    xicent=xicent+xsep/pi2*cos(xn(i)/xsep*pi2)
                    yicent=yicent+ysep/pi2*cos(yn(i)/ysep*pi2)
                    zicent=zicent+zsep/pi2*cos(zn(i)/zsep*pi2)
                    xjcent=xjcent+xsep/pi2*sin(xn(i)/xsep*pi2)
                    yjcent=yjcent+ysep/pi2*sin(yn(i)/ysep*pi2)
                    zjcent=zjcent+zsep/pi2*sin(zn(i)/zsep*pi2)
                enddo
                xicent=xicent/(menchain(ic+1)-menchain(ic))
                yicent=yicent/(menchain(ic+1)-menchain(ic))
                zicent=zicent/(menchain(ic+1)-menchain(ic))
                xjcent=xjcent/(menchain(ic+1)-menchain(ic))
                yjcent=yjcent/(menchain(ic+1)-menchain(ic))
                zjcent=zjcent/(menchain(ic+1)-menchain(ic))
                if(lpbcx) xcentr=xsep/pi2*(atan2(-xjcent,-xicent)+Pi)
                if(lpbcy) ycentr=ysep/pi2*(atan2(-yjcent,-yicent)+Pi)
                if(lpbcz) zcentr=zsep/pi2*(atan2(-zjcent,-zicent)+Pi)
                do i=menchain(ic)+1,menchain(ic+1)
                    if(lpbcx.and.abs(xn(i)-xcentr).gt.0.5*xsep) then
                        if(xn(i).gt.xcentr) then
                            xn(i) = xn(i)-xsep
                        else
                            xn(i) = xn(i)+xsep
                        endif
                    endif
                    if(lpbcy.and.abs(yn(i)-ycentr).gt.0.5*ysep) then
                        if(yn(i).gt.ycentr) then
                            yn(i) = yn(i)-ysep
                        else
                            yn(i) = yn(i)+ysep
                        endif
                    endif
                    if(lpbcz.and.abs(zn(i)-zcentr).gt.0.5*zsep) then
                        if(zn(i).gt.zcentr) then
                            zn(i) = zn(i)-zsep
                        else
                            zn(i) = zn(i)+zsep
                        endif
                    endif
                enddo
                do i=menchain(ic)+2,menchain(ic+1)
                    j=i-1
                    xdif=abs(xn(i)-xn(j))
                    ydif=abs(yn(i)-yn(j))
                    zdif=abs(zn(i)-zn(j))
                    lcheck=.false.
                    lcheck=lcheck.or.xdif.gt.0.5*xsep
                    lcheck=lcheck.or.ydif.gt.0.5*ysep
                    lcheck=lcheck.or.zdif.gt.0.5*zsep
                    if(lcheck) then
             if(i-menchain(ic).gt.(menchain(ic+1)-menchain(ic))/2) then
                            i3=j
                            kleftend=i
                            krightend=menchain(ic+1)
                        else
                            i3=i
                            kleftend=menchain(ic)+1
                            krightend=j
                        endif
                        do i2=kleftend,krightend
                            if(xdif.gt.0.5*xsep) then
                                if(xn(i2).gt.xn(i3)) then
                                    xn(i2)=xn(i2)-xsep
                                else
                                    xn(i2)=xn(i2)+xsep
                                endif
                            endif
                            if(ydif.gt.0.5*ysep) then
                                if(yn(i2).gt.yn(i3)) then
                                    yn(i2)=yn(i2)-ysep
                                else
                                    yn(i2)=yn(i2)+ysep
                                endif
                            endif
                            if(zdif.gt.0.5*zsep) then
                                if(zn(i2).gt.zn(i3)) then
                                    zn(i2)=zn(i2)-zsep
                                else
                                    zn(i2)=zn(i2)+zsep
                                endif
                            endif
                        enddo
                    endif
                enddo
                xcentr=0.0
                ycentr=0.0
                zcentr=0.0
                do i=menchain(ic)+1,menchain(ic+1)
                    xcentr=xcentr+xn(i)
                    ycentr=ycentr+yn(i)
                    zcentr=zcentr+zn(i)
                enddo
                xcentr=xcentr/(menchain(ic+1)-menchain(ic))
                ycentr=ycentr/(menchain(ic+1)-menchain(ic))
                zcentr=zcentr/(menchain(ic+1)-menchain(ic))
                if(lpbcx.and.abs(xcentr-xcentrall).gt.0.5*xsep) then
                    do i=menchain(ic)+1,menchain(ic+1)
                        if(xcentr.gt.xcentrall) then
                            xn(i) = xn(i)-xsep
                        else
                            xn(i) = xn(i)+xsep
                        endif
                    enddo
                endif
                if(lpbcy.and.abs(ycentr-ycentrall).gt.0.5*ysep) then
                    do i=menchain(ic)+1,menchain(ic+1)
                        if(ycentr.gt.ycentrall) then
                            yn(i) = yn(i)-ysep
                        else
                            yn(i) = yn(i)+ysep
                        endif
                    enddo
                endif
                if(lpbcz.and.abs(zcentr-zcentrall).gt.0.5*zsep) then
                    do i=menchain(ic)+1,menchain(ic+1)
                        if(zcentr.gt.zcentrall) then
                            zn(i) = zn(i)-zsep
                        else
                            zn(i) = zn(i)+zsep
                        endif
                    enddo
                endif
            enddo
        endif
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOAD_SEQUENCE loads amino acid sequences to make random configurations
        subroutine load_sequence(filename)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lfrompdb(len),lmap,lsimpang
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        character aseq*3,ares*3,filename*32,bb*2,buffer*1024,mapfile*32
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        open(8,file=filename,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,a)')'ERROR OPENING FILENAME ',filename
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        ib=0
        
        read(8,'(a)') buffer
        if(buffer(1:7).eq.'screend') then
            read(buffer(8:),*) screend
            read(8,*) nchains ! number of protein chains to read
        else
            read(buffer,*) nchains
        endif
        
        menchain(1)=0
        do ic=2,nchains+1
        read(8,*) menchain(ic) ! number of amino acids in protein chain
        read(8,*) buffer ! sequence in FASTA format, smaller than 1024aa
        do ibb=1,menchain(ic)
        ib=ib+1
        lmap=.false.
        lconect(ib)=.true.
        if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'G') then
        aseq(ib) = 'GLY'
        inameseq(ib)=1
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'P') then
        aseq(ib) = 'PRO'
        inameseq(ib)=2
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'Q') then
        aseq(ib) = 'GLN'
        inameseq(ib)=3
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'C') then
        aseq(ib) = 'CYS'
        inameseq(ib)=4
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'A') then
        aseq(ib) = 'ALA'
        inameseq(ib)=5
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'S') then
        aseq(ib) = 'SER'
        inameseq(ib)=6
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'V') then
        aseq(ib) = 'VAL'
        inameseq(ib)=7
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'T') then
        aseq(ib) = 'THR'
        inameseq(ib)=8
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'I') then
        aseq(ib) = 'ILE'
        inameseq(ib)=9
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'L') then
        aseq(ib) = 'LEU'
        inameseq(ib)=10
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'N') then
        aseq(ib) = 'ASN'
        inameseq(ib)=11
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'D') then
        aseq(ib) = 'ASP'
        inameseq(ib)=12
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'K') then
        aseq(ib) = 'LYS'
        inameseq(ib)=13
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'E') then
        aseq(ib) = 'GLU'
        inameseq(ib)=14
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'M') then
        aseq(ib) = 'MET'
        inameseq(ib)=15
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'H') then
        aseq(ib) = 'HIS'
        inameseq(ib)=16
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'F') then
        aseq(ib) = 'PHE'
        inameseq(ib)=17
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'R') then
        aseq(ib) = 'ARG'
        inameseq(ib)=18
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'Y') then
        aseq(ib) = 'TYR'
        inameseq(ib)=19
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'W') then
        aseq(ib) = 'TRP'
        inameseq(ib)=20
        else if(buffer(ib-menchain(ic-1):ib-menchain(ic-1)).eq.'X') then
        inameseq(ib)=0
        endif
        iseq(ib)=ib-menchain(ic-1)
        xn(ib)=0.0
        yn(ib)=0.0
        zn(ib)=0.0
        enddo
        do ibb=1,menchain(ic)
            the0(ibb)=sigma1(min(inameseq(ibb),3))
        enddo
        menchain(ic)=menchain(ic)+menchain(ic-1)
        lconect(ib)=.false.
        read(8,*) nmaps
        if(nmaps.gt.0) then
            do imaps=1,nmaps
            read(8,*) mapfile ! map file for given sequence
            open(28,file=mapfile,status='old',iostat=ierr)
            if(ierr.ne.0) then
            write(*,'(a,a)')'ERROR OPENING FILENAME ',mapfile
            write(*,*)'PROGRAM STOPPED.'
            stop
            endif
            read(28,*) iframeshift ! how many residues from the start
            read(28,*) npdblen ! how long is the chain fragment with map
            read(28,*) ncntact ! how many contacts there are
            do k=1,ncntact
                klont=klont+1
                read(28,*) i,j,sont(klont)
                klist(1,klont)=i+menchain(ic-1)+iframeshift
                klist(2,klont)=j+menchain(ic-1)+iframeshift
                klist(3,klont)=1 ! L-J TODO minus 1 for different chains
            enddo
            do i=1,npdblen
                ibtrue=i+menchain(ic-1)+iframeshift
                read(28,*) the0(ibtrue),phi0(ibtrue)
                lfrompdb(ibtrue)=.true.
            enddo
            close(28)
            enddo
        endif
        enddo
        close(8)
        men=ib
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE LOADS PARAMETERS FOR NON-SPECIFIC ANGLE POTENTIALS
        subroutine load_paramfile(paramfile)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lfrompdb(len),lsimpang,ldi,lconftm,lcpot,lsim,lradii
        logical lpid,lenetab,lmj
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/mr/tene(18001,3,3,2),pp,lenetab
        dimension sig1r(24)
        character paramfile*32,buffer*128,aa1sig*3,aa2sig*3
        
        open(88,file=paramfile,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,2x,a)')'ERROR OPENING FILENAME',paramfile
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        ! 1=GLY, 2=PRO, 3=X
c        if(lsimpang) then
            read(88,*)
            read(88,*)(sigma1(i),i=1,7)
            sigma1(8)=sigma1(7)
            sigma1(7)=sigma1(6)
c        else
            read(88,*)
            read(88,*)(angpot(i,1,1),i=1,7)
            read(88,*)(angpot(i,1,2),i=1,7)
            read(88,*)(angpot(i,1,3),i=1,7)
            read(88,*)(angpot(i,2,1),i=1,7)
            read(88,*)(angpot(i,2,2),i=1,7)
            read(88,*)(angpot(i,2,3),i=1,7)
            read(88,*)(angpot(i,3,1),i=1,7)
            read(88,*)(angpot(i,3,2),i=1,7)
            read(88,*)(angpot(i,3,3),i=1,7)
c        endif
c        if(ldi) then
            read(88,*)
            read(88,*)(dihpot(i,1,1),i=1,6)
            read(88,*)(dihpot(i,1,2),i=1,6)
            read(88,*)(dihpot(i,1,3),i=1,6)
            read(88,*)(dihpot(i,2,1),i=1,6)
            read(88,*)(dihpot(i,2,2),i=1,6)
            read(88,*)(dihpot(i,2,3),i=1,6)
            read(88,*)(dihpot(i,3,1),i=1,6)
            read(88,*)(dihpot(i,3,2),i=1,6)
            read(88,*)(dihpot(i,3,3),i=1,6)
c        if(lpid) then
c            read(88,*)(dihpot(i,4,4),i=1,6)
c            read(88,*)(dihpot(i,5,5),i=1,6)
c            do i=1,6
c                dihpot(i,6,6)=dihpot(i,4,4)
c                dihpot(i,7,7)=dihpot(i,5,5)
c            enddo
c        endif
            do i=1,3
                do j=1,3
c                    do k=1,7 ! derivative
c                        angpot(k+7,j,i)=(k-1.0)*angpot(k,j,i)
                        angpot(8,j,i)=0.d0
                        angpot(9,j,i)=angpot(2,j,i)
                        angpot(10,j,i)=2.d0*angpot(3,j,i)
                        angpot(11,j,i)=3.d0*angpot(4,j,i)
                        angpot(12,j,i)=4.d0*angpot(5,j,i)
                        angpot(13,j,i)=5.d0*angpot(6,j,i)
                        angpot(14,j,i)=6.d0*angpot(7,j,i)
c                    enddo
                enddo
            enddo
c        endif
        read(88,*)
        read(88,*)(ksdchn(i,1),i=1,20) ! aatype
        read(88,*)(ksdchn(i,2),i=1,20) ! aatype s limit
        read(88,*)(ksdchn(i,3),i=1,20) ! aatype s limit for nonpolar aa
        read(88,*)(ksdchn(i,4),i=1,20) ! aatype s limit for polar aa
c        if(lcpot.and.(.not.lsim)) then
        if(.not.lsim) then
                read(88,*)
                read(88,*)(sig1r(i),i=1,20)
            if(lmj.or..not.lradii) then
                read(88,*)
                do i=1,20
                    do j=i,20
                        if(lmj) then
                            read(88,*) aa1sig,aa2sig,sigmaij,epsilonij
                            emj(i*21+j)=epsilonij
                            emj(j*21+i)=epsilonij
                        else
                            read(88,*) aa1sig,aa2sig,sigmaij
                        endif
                        sigma1(j*21+i)=sigmaij/unit*0.5d0**(1.d0/6.d0)
                        sigma1(i*21+j)=sigmaij/unit*0.5d0**(1.d0/6.d0)
                    enddo
                enddo
            endif
            if(lradii) then
                do i=1,20
                    do j=i,20
            sigma1(j*21+i)=(sig1r(i)+sig1r(j))/unit*0.5d0**(1.d0/6.d0)
            sigma1(i*21+j)=(sig1r(i)+sig1r(j))/unit*0.5d0**(1.d0/6.d0)
                    enddo
                enddo
            endif
        endif
        close(88)
        
        if (lenetab) then ! MR vvvvvvvvvvvvv
            open(88,file='backbone.txt',status='old',iostat=ierr)
            if(ierr.ne.0) then
                write(*,'(a)')'ERROR OPENING FILENAME backbone.txt'
                write(*,*)'PROGRAM STOPPED.'
                stop
            endif
            do it=1,18001
              read(88,*)(((tene(it,j,k,l),k=1,3),j=1,3),l=1,2)
            enddo
            close(88)
            pp=atan(1.0)/4500
        endif ! MR ^^^^^^^^^^^^^^^^
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine kmt(klewy,kprawy,kbegin,kend)
       implicit double precision(a-h,o-z)
       parameter(len=10000)
       common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
       dimension xk(kend-kbegin+1)
       dimension yk(kend-kbegin+1)
       dimension zk(kend-kbegin+1)
       dimension na(kend-kbegin+1)
       dimension nb(kend-kbegin+1)
       knot=0
       knot1=0
       kbinl=0
       kbinr=0
       kendl=0
       kendr=0
       koba=0

       nres=kend-kbegin+1

       do i=1,nres!kbegin,kend
         xk(i)=x0(i+kbegin-1)
         yk(i)=y0(i+kbegin-1)
         zk(i)=z0(i+kbegin-1)
         na(i)=i+kbegin-1
         nb(i)=i+kbegin-1
       enddo

       ii=0
       j=0
       i=2
       n=nres
       modify=0

566    continue   
       ifun3=0
         det0=0d0
         det1=0d0
         det2=0d0
         ifun=0
         ro=0

!#######################################################################

         do j=1,n-1
       if (((j.lt.(i-2)).or.(j.gt.(i+1))).and.(ifun.ne.1)) then
    
          ro=0.d0
          xro=0.d0
          yro=0.d0
          zro=0.d0
        
        det1=(xk(j)-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(zk(j)-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(yk(j)-yk(i-1))*(zk(i)-zk(i-1))-
     & (zk(j)-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(xk(j)-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(yk(j)-yk(i-1))*(xk(i)-xk(i-1))
        
        det2=(xk(j+1)-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(zk(j+1)-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(yk(j+1)-yk(i-1))*(zk(i)-zk(i-1))-
     & (zk(j+1)-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(xk(j+1)-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(yk(j+1)-yk(i-1))*(xk(i)-xk(i-1))

        if (det1*det2.gt.0) ifun=0

        if (det1*det2.le.0) then 
         ! przeciecie sie plaszczyzny trojkata i-1,i,i+1
         ! z punktem wyzanaczonym przez prosta j,j+1
         ! dane wzorem na ro, F0_(i,1,i,i+1)(0)

       xji=xk(j+1)-xk(j)     ! wyznaczenie wektora od j do j+1
       yji=yk(j+1)-yk(j)
       zji=zk(j+1)-zk(j)
       
       det3=(xji-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(zji-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(yji-yk(i-1))*(zk(i)-zk(i-1))-
     & (zji-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(xji-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(yji-yk(i-1))*(xk(i)-xk(i-1))
       
       det0=(-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(-yk(i-1))*(zk(i)-zk(i-1))-
     & (-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(-yk(i-1))*(xk(i)-xk(i-1))
         
       ro=det1/(det3-det0)
       
       xro=xk(j)-ro*xji 
       yro=yk(j)-ro*yji
       zro=zk(j)-ro*zji
       
       f12xy=(xro-xk(i-1))*(yk(i)-yk(i-1))-(xk(i)-xk(i-1))*(yro-yk(i-1))
       f23xy=(xro-xk(i))*(yk(i+1)-yk(i))-(xk(i+1)-xk(i))*(yro-yk(i))
       f31xy=(xro-xk(i+1))*(yk(i-1)-yk(i+1))-(xk(i-1)-xk(i+1))*(yro-
     & yk(i+1))

        if ((f12xy*f23xy.le.0).or.(f23xy*f31xy.le.0)) then 
         ifun=0

        elseif ((f12xy*f23xy.gt.0).and.(f23xy*f31xy.gt.0)) then
          ifun=1

        endif
       endif 
       endif          
         enddo 


!#######################################################################

       if (ifun.eq.0) then 

       do ii=i,n-1
        xk(ii)=xk(ii+1)
        yk(ii)=yk(ii+1)
        zk(ii)=zk(ii+1)
        na(ii)=na(ii+1)
       enddo
          
        xk(n)=0d0
        yk(n)=0d0
        zk(n)=0d0
        modify=modify+1
        n=n-1        
       elseif (ifun.eq.1) then
        i=i+1
       endif
!#######################################################################

       if (i.lt.n) then
        goto 566
       elseif ((i.eq.n).and.(modify.ge.1)) then
        i=2
        modify=0
        if(n.ge.3) goto 566
       endif

!#######################################################################

       if (i.eq.2) then
!-----------------------------------------------------------------------         
         if (knot.eq.0) knot1=1
         knot=1
!-----------------------------------------------------------------------       
         if ((kbinl.eq.0).and.(kbinr.eq.0)) then
            kbinl=0
            kbinr=0
         elseif (kendl.eq.0) then
            kendl=1     
            n=nres-kbinr      
            do i=1,n-1
             xk(i)=x0(i+kbegin-1)
             yk(i)=y0(i+kbegin-1)
             zk(i)=z0(i+kbegin-1)
            enddo
            n=n-1  
            kbinr=kbinr+1  
            i=2
            modify=0
            goto 566
         elseif (kendr.eq.0) then
            kendr=1      
            n=nres-koba  
            n=n-2     
            koba=koba+1 
            do i=1,n
             xk(i)=x0(i+kbegin-1+1)
             yk(i)=y0(i+kbegin-1+1)
             zk(i)=z0(i+kbegin-1+1)
            enddo
            i=2
            modify=0
            goto 566
         elseif (koba.gt.0) then
            koba=koba
         endif
       
!-----------------------------------------------------------------------       

       elseif (i.ne.2) then
!-----------------------------------------------------------------------             
         if (knot.eq.1) knot1=1
         knot=0
!-----------------------------------------------------------------------               
         if ((kbinl.lt.200).and.(kendl.eq.0)) then
            n=nres-kbinl
            n=n-1
            kbinl=kbinl+1
            do i=1,n
             xk(i)=x0(i+kbegin-1+kbinl)
             yk(i)=y0(i+kbegin-1+kbinl)
             zk(i)=z0(i+kbegin-1+kbinl)
            enddo
            i=2
            modify=0
            goto 566
         elseif ((kbinl.ge.200).and.(kendl.eq.0)) then
            kendl=1  
            n=nres-kbinr
            n=n-1
            kbinr=kbinr+1
            do i=1,n
             xk(i)=x0(i+kbegin-1)
             yk(i)=y0(i+kbegin-1)
             zk(i)=z0(i+kbegin-1)
            enddo
            i=2
            modify=0
            goto 566
         elseif ((kbinr.lt.200).and.(kendr.eq.0)) then
            n=nres-kbinr
            n=n-1        
            kbinr=kbinr+1  
            do i=1,n
             xk(i)=x0(i+kbegin-1)
             yk(i)=y0(i+kbegin-1)
             zk(i)=z0(i+kbegin-1)
            enddo
            i=2
            modify=0
            goto 566
         elseif ((kbinr.ge.200).and.(kendr.eq.0)) then
            kendr=1  
            n=nres-koba
            koba=koba+1
            n=n-2
            do i=1,n
             xk(i)=x0(i+kbegin-1+1)
             yk(i)=y0(i+kbegin-1+1)
             zk(i)=z0(i+kbegin-1+1)
            enddo
            i=2
            modify=0
            goto 566
         elseif (koba.le.1) then
            n=nres-(koba*2)
            koba=koba+1
            n=n-2
            do i=1,n-1
             xk(i)=x0(i+kbegin-1+koba+1)
             yk(i)=y0(i+kbegin-1+koba+1)
             zk(i)=z0(i+kbegin-1+koba+1)
            enddo
            i=2
            modify=0
            goto 566
         endif
!-----------------------------------------------------------------------       
       endif
!#######################################################################     
      
       kprawy=nres-kbinr-1
       klewy=kbinl-1
       return
       end
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINES COMPUTE THE RMSD OF THE C ALPHA BACKBONE TO THE
C NATIVE BACKBONE TAKEN FROM PDB
        subroutine compute_rmsd(rms)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        dimension r(len,3),rn(len,3)
        do ib=1,men
        r(ib,1)=x0(ib)
        r(ib,2)=y0(ib)
        r(ib,3)=z0(ib)
        rn(ib,1)=xn(ib)
        rn(ib,2)=yn(ib)
        rn(ib,3)=zn(ib)
        enddo
        call kabsch(men,r,rn,rms)
        return
        end

      subroutine kabsch(ns,q,q_trg,d)
      implicit double precision(a-h,o-z)
      parameter(len=10000)
      parameter(nsx=len)
      dimension q(nsx,3),q_trg(nsx,3),qk(nsx,3),qkk(nsx,3)
      dimension r(3,3),s(3,3),u(3,3),a(3,3),b(3,3),c(3,3)
      dimension q_trg_cm(3),q_cm(3),aumu(3),v1(3),v2(3),v3(3)
c*********************************************************************
c                    no translations
c*********************************************************************
      e_0=0.E0
      do k=1,3
        q_trg_cm(k) =  0.E0
        q_cm(k) = 0.E0
        do i=1,ns
          q_trg_cm(k) = q_trg_cm(k) + q_trg(i,k)
          q_cm(k)     = q_cm(k) + q(i,k)
        end do
        q_trg_cm(k) = q_trg_cm(k)/ns
        q_cm(k)     = q_cm(k)/ns
        do i=1,ns
          qk(i,k)    = q(i,k) - q_cm(k)
          q_trg(i,k) = q_trg(i,k) - q_trg_cm(k)
c           now they have the same center of mass: (0,0,0)
          e_0=e_0+qk(i,k)*qk(i,k)+q_trg(i,k)*q_trg(i,k)
        end do
      end do
      e_0=e_0/(1.0*ns)
c   no more need for q, just qk!!!  calcolo matrice r
      do k=1,3
         do l=1,3
            r(k,l)=0.E0
            do i=1,ns
               r(k,l) = r(k,l) + q_trg(i,k)*qk(i,l)
            end do
            r(k,l)=2.E0*r(k,l)/(1.0*ns)      
         end do    
      end do
      
c*****************************************************
c    calculate matrix s = r_t*r
c*****************************************************
      do k=1,3
         do l=1,3
            s(k,l)=0.E0      
            do m=1,3
               s(k,l)=s(k,l)+r(m,k)*r(m,l)
            end do
         end do
      end do
      call jacobi(s,3,3,aumu,a,nrot)
      call eigsrt(aumu,a,3,3)
c calcolo a3 come a1xa2
      do i=1,3
         v1(i)=a(i,1)
         v2(i)=a(i,2)
      end do
      call pvector (v1,v2,v3,vmod)
      do k=1,3
         a(k,3)=v3(k)
      end do
c     calcolo i vettori b
      do k=1,3
         do l=1,3
            b(k,l)=0.E0
            do m=1,3
               b(k,l)=b(k,l)+r(k,m)*a(m,l)
            end do
            c(k,l)=b(k,l)
         end do
      end do
      
c********************************************************************      
c           normalization
c********************************************************************
      call norma(b)
      do i=1,3
        v1(i)=b(i,1)
        v2(i)=b(i,2)
      end do
      call pvector (v1,v2,v3,vmod)
      do i=1,3
         b(i,3)=v3(i)
      end do
c controllo i prodotti scalari degli a_k con i b_k 
      bdotc=0.E0
      do l=1,3
         bdotc=bdotc+b(l,3)*c(l,3)
      end do
c se il prodotto scalare del terzo e'<0, allora cambia il segno
      sigma3 = 1.0
      if (bdotc.lt.0.E0) sigma3=-1.0

c***********************************************************************
c                matrice u
c***********************************************************************
      do l=1,3
         do m=1,3
            u(l,m)=0.E0
            do k=1,3
               u(l,m)=u(l,m)+b(l,k)*a(m,k)
            end do
         end do
      end do

c******************** rotazione di kabsch   ****************************
      do i=1,ns
        do k=1,3
          qkk(i,k)=0.E0
          do l=1,3
             qkk(i,k)=qkk(i,k)+u(k,l)*qk(i,l)
          end do
        end do
      end do
c***********************************************************************
c             distanza di kabsch
c***********************************************************************
      d=0.E0
      do i=1,ns
c         print *,(qkk(i,k),k=1,3)  ! serve per scrivere le coordinate
        do k =1,3
          d=d+(q_trg(i,k)-qkk(i,k))**2
        end do
      end do
      d=sqrt(d/ns)
      e = sqrt(e_0-sqrt(aumu(1))-sqrt(aumu(2))-sigma3*sqrt(aumu(3)))
      return
      end

      subroutine norma (a)
      implicit double precision(a-h,o-z)
      dimension a(3,3)
      do i=1,3
         sum = 0.0 
         do j=1,3
            sum = sum + a(j,i)*a(j,i)
         end do
         do j=1,3
            a(j,i) = a(j,i)/sqrt(sum)
         end do 
      end do              
      return
      end 

        subroutine pvector(u,v,z,zmod)
        implicit double precision(a-h,o-z)
        dimension u(3),v(3),z(3)
        z(1)=u(2)*v(3)-u(3)*v(2)
        z(2)=u(3)*v(1)-u(1)*v(3)
        z(3)=u(1)*v(2)-u(2)*v(1)
        zmod=sqrt(z(1)*z(1)+z(2)*z(2)+z(3)*z(3))
        return
        end

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      real*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      real*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      write(*,*) 'too many iterations in jacobi'
      stop
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software W"..

      SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np
      real*8 d(np),v(np,np)
      INTEGER i,j,k
      real*8 p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software W"..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        SUBROUTINE sort(n,arr)
        implicit double precision (a-h,o-z)
        INTEGER n,M,NSTACK         !    from "Numerical Recipes"
        dimension arr(n)
        PARAMETER (M=7,NSTACK=500)
        INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
C       REAL a,temp
        jstack=0
        l=1
        ir=n
1       if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
          if(arr(i).le.a)goto 2
          arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
        else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          write(*,*)'NSTACK too small in sort'
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
        endif
        goto 1
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE PRINTS THE ORIGINAL CONTACT MAP FILE
        subroutine print_cmap(iun,ktime)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lfrompdb(len)
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        
c        write(iun,*) "# total chain length",men
        write(iun,'(a,i9)') '# number of contacts',klont
        write(iun,'(a,i7)') ' time i j nr type Ri Rj r, length = ',men
        do i=1,klont
            sigrang=sont(i)*c216
            ib1=klist(1,i)
            ib2=klist(2,i)
            write(iun,'(a,i9,2i5,i9,i3,x,a,x,a,f7.3)')"K",ktime,ib1,ib2,
     +          ib1*men+ib2,klist(3,i),aseq(ib1),aseq(ib2),sigrang*unit
        enddo
        
        do i=1,men
               write(iun,'(a,i9,2f8.4,x,a)') 
     +         "A",ktime,the0(i),phi0(i),aseq(i)
        enddo
        
        write(iun,*)
        write(iun,*)
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C THIS SUBROUTINE PRINTS THE CONTACT MAP FILE
        subroutine print_map(iun,ktime)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lpbcx,lpbcy,lpbcz,l3rdcn(len)
        logical lfrompdb(len),lconect(len),lchiral,langle,ldet,lwrtang
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/restr/bckbmin,bckb2min,sdchnmax,sfact,kmaxbckb,lecperm
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
c        write(iun,*) "# total chain length",men

        ncnt=0
        do k=1,kqont ! non-native attractive contacts
            if(abs(kqist(4,k,jq)).gt.8) then
                iconttype=5
            else
                iconttype=abs(kqist(4,k,jq)) ! 4-8 for backbone contact
            endif
            rsig=sigma1(abs(kqist(4,k,jq)))
            i=kqist(1,k,jq) ! kqist 1 or 2 = residue numbers
            j=kqist(2,k,jq) ! jq = number of verlet list (1 or 2)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.le.(rsig*cntfct)**2 .and.iconttype.gt.1) then
                ncnt=ncnt+1
                nli(1,ncnt)=i
                nli(2,ncnt)=j
                nli(3,ncnt)=0
                nli(4,ncnt)=iconttype
                nli(5,ncnt)=kqist(4,k,jq)
            endif
        enddo

        
        do k=1,kcont ! all other contacts (including native)
            i=kcist(1,k)
            j=kcist(2,k)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            icm=kcist(3,k)
            if(icm.gt.0) then ! NATIVE CONTACTS
                if(abs(klist(3,icm)).lt.631) then ! L-J potential
                    rsig=sont(icm)
                    if(rsq.le.(rsig*cntfct)**2) then
                        ncnt=ncnt+1
                        nli(1,ncnt)=i
                        nli(2,ncnt)=j
                        nli(3,ncnt)=abs(klist(3,icm))
                        nli(4,ncnt)=0
                        nli(5,ncnt)=sign(icm,klist(3,icm))
                    endif
                endif
            else
                kremainder=mod(icm,2)
                icm=icm/2
                if(icm.eq.-1) then
                    if(rsq.le.(rmrs*cntfct)**2) then ! FOR SS BONDS
                        ncnt=ncnt+1
                        nli(1,ncnt)=i
                        nli(2,ncnt)=j
                        nli(3,ncnt)=-1
                        nli(4,ncnt)=-1
                        if(kremainder.eq.0) then
                            nli(5,ncnt)=88
                        else
                            nli(5,ncnt)=-88
                        endif
                        if(l3rdcn(i).or.l3rdcn(j)) then
                            if(knct3rd(i).eq.j.and.knct3rd(j).eq.i)then
                                nli(4,ncnt)=3
                            else ! OTHER CYS ARE NOT ATTRACTED TO SSBOND
                                nli(4,ncnt)=2
                            endif
                        endif
                    endif
                endif
            endif
        enddo
        
        write(iun,'(a,i9)') '# number of contacts',ncnt
        write(iun,'(5(a,f9.4))')'bbij ',bckbmin,' sdchnmax^2 ',sdchnmax,
     +  ' ssPID ',psi0ss,' bbPID1 ',psi0bb(1),' bbPID2 ',psi0bb(2)
        write(iun,'(2a)')'         time i j nr type Ri Rj r bbir bbjr ',
     +  'bbij sdchni^2 sdchnj^2 conttype'
        do i=1,ncnt
            if(nli(4,i).eq.0) then
                sigrang=sont(abs(nli(5,i)))*c216
            else
                sigrang=sigma1(abs(nli(5,i)))*c216
            endif
            ib1=nli(1,i)
            ib2=nli(2,i)
            if(ldet) then
                call compute_details(ib1,ib2)
        !sigrang=sigma1(inameseq(ib1)*21+inameseq(ib2))*c216
         write(iun,'(a,i9,2i6,i9,i4,x,a,x,a,7f7.3,i4,2f7.3,2i3,i7)')"K",
     +  ktime,ib1,ib2,ib1*men+ib2,nli(4,i),aseq(ib1),aseq(ib2),rangli(1)
     +  ,rangli(2),rangli(3),rangli(4),rangli(5),rangli(6),sigrang*unit,
     +  nli(3,i),rangli(7),rangli(8),nei(1,ib1),nei(1,ib2),nli(5,i)
            else
                dx=x0(ib1)-x0(ib2)
                dy=y0(ib1)-y0(ib2)
                dz=z0(ib1)-z0(ib2)
                if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                if(lpbcz) dz = dz-zsep*nint(dz*zinv)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                write(iun,'(a,i9,2i5,i9,i3,x,a,x,a,f7.3)')"K",ktime,ib1,
     +          ib2,ib1*men+ib2,nli(4,i),aseq(ib1),aseq(ib2),r*unit
            endif
        enddo
        !ind=1
        if(lwrtang) then
          if(lchiral) then
            do i=2,men-3
            if(lconect(i-1).and.lconect(i).and.lconect(i+1)) then
            write(iun,*) "A",ktime,thetemp(i+1),phitemp(i+1),chir(i)
            !ind=ind+1
            endif
            enddo
          else
c            write(iun,*) 0.0,0.0
c            write(iun,*) the0(2),0.0
            do i=3,men-2
               if(lconect(i-2).and.lconect(i-1).and.lconect(i)) then
               write(iun,'(a,i9,2f8.4,x,a)') 
     +         "A",ktime,thetemp(i),phitemp(i),aseq(i)
               !ind=ind+1
               endif
            enddo
c            write(iun,*) 0.0,0.0
          endif
        endif
        write(iun,*)
        write(iun,*)
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE PRINTS THE CHAIN COORDINATES IN PDB FORMAT
        subroutine print_conformation(iun,time,energy,rms)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character*3 aseq,ares
        character*2 chainid(len)
        logical lwal,lwals,ldens,lsim,lconect(len),l3rdcn(len),lposcrd
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
        dimension nserial(len)
        
        if(lwal) then
            xbox=(xup-xdown)*unit
            ybox=(yup-ydown)*unit
            zbox=(zup-zdown)*unit
            xmin=xdown
            ymin=ydown
            zmin=zdown
        else
            xmin=x0(1)
            ymin=y0(1)
            zmin=z0(1)
            do ib=2,men
                if(x0(ib).lt.xmin) xmin=x0(ib)
                if(y0(ib).lt.ymin) ymin=y0(ib)
                if(z0(ib).lt.zmin) zmin=z0(ib)
            enddo
        endif
        
        if(lposcrd) then
          write(iun,'(/,2a,f10.1,2x,a,f12.2,2x,a,3f9.3)') 'REMARK   0 ', 
     +    'TIME =',time,'ENERGY =',energy,'SHIFT =',xmin,ymin,zmin
        else
          write(iun,'(/,a,f10.1,2x,a,f12.2,2x,a,f9.4)')
     +   'REMARK   0 TIME =',time,'ENERGY =',energy,'RMSD =',rms*unit
          xmin=0.d0
          ymin=0.d0
          zmin=0.d0
        endif
        if(lwal) write(iun,'(a,3f9.3,3f7.2,x,a,11x,a)')
     +      'CRYST1',xbox,ybox,zbox,90.0,90.0,90.0,'P 1','1'
        
        ic=1
        iatom=0
        cmin=-999/unit
        do 1000 ib=1,men
        ares=aseq(ib)
        iatom=iatom+1
        ic2=ic
1532    continue
        if(ic2 .gt. 153) then ! program supports up to 153 unique chains
            ic2=ic2-153
            goto 1532
        endif
        if(ic2 .le. 26) then ! A IS ASCII 65, Z IS 90 (26 letters)
            icchar=ic2+64 ! A,B,C,D...
        else
            if(ic2 .le. 52) then
                icchar=ic2+70 ! a,b,c,d...
            else
                icchar=ic2-5 ! 0,1,2,3...
            endif
        endif
        if(ic2 .le. 62) then
            write(chainid(ib),'(a,a)') ' ',char(icchar)
        else
            write(chainid(ib),'(i2)') ic2-53 ! 10,11,12,13...
        endif
        nserial(ib)=iatom
        if(x0(ib)-xmin.lt.cmin.or.y0(ib)-ymin.lt.cmin.or.
     +  z0(ib)-zmin.lt.cmin) then
        write(iun,'(a,i7,2x,a,a,a,i4,f12.2,2f8.2,f6.2)')
     +  'ATOM',iatom,'CA  ',ares,chainid(ib),iseq(ib),(x0(ib)-xmin)*unit
     +  ,(y0(ib)-ymin)*unit,(z0(ib)-zmin)*unit,float(nei(1,ib))
        else
        write(iun,'(a,i7,2x,a,a,a,i4,f12.3,2f8.3,f6.2)')
     +  'ATOM',iatom,'CA  ',ares,chainid(ib),iseq(ib),(x0(ib)-xmin)*unit
     +  ,(y0(ib)-ymin)*unit,(z0(ib)-zmin)*unit,float(nei(1,ib))
        endif
        if(.not. lconect(ib)) then
        iatom=iatom+1
        write(iun,'(a,i7,6x,a,a,i6,28x)') 'TER ',iatom,ares,
     +  chainid(ib),iseq(ib)
        iatom=iatom+1
        if(xyzcm(1,ic)-xmin.lt.cmin.or.xyzcm(2,ic)-ymin.lt.cmin.or.
     +  xyzcm(3,ic)-zmin.lt.cmin) then
        write(iun,'(a,i5,2x,3a,i4,f12.2,2f8.2)') 'HETATM',iatom,'C   ',
     +  'COG',chainid(ib),iseq(ib)+1,(xyzcm(1,ic)-xmin)*unit,
     +  (xyzcm(2,ic)-ymin)*unit,(xyzcm(3,ic)-zmin)*unit
        else
        write(iun,'(a,i5,2x,3a,i4,f12.3,2f8.3)') 'HETATM',iatom,'C   ',
     +  'COG',chainid(ib),iseq(ib)+1,(xyzcm(1,ic)-xmin)*unit,
     +  (xyzcm(2,ic)-ymin)*unit,(xyzcm(3,ic)-zmin)*unit
        endif
        write(iun,'(a,i4,a,f10.1,a,f8.2,a,f7.2,a,i4,a,i3,a,i3,a,f6.3)') 
     +  'REMARK',ic,' T= ',time,' RG= ',cgyr(ic)*unit,
     +  ' R_end_to_end= ',cend(ic)*unit,' N= ',menchain(ic+1)
     +  -menchain(ic),' K1= ',knts(1,ic),' K2= ',knts(2,ic),' W= ',w(ic)
c       write(iun,'(a,i4,a,f10.1,i6,3f8.3)')'COG   ',ic,' T= ',time,
c    +  menchain(ic+1)-menchain(ic),(xyzcm(1,ic)-xmin)*unit,
c    +  (xyzcm(2,ic)-ymin)*unit,(xyzcm(3,ic)-zmin)*unit
        ic=ic+1
        endif
1000    continue
        nssb2=0
        do 1100 ib=1,men
        if(l3rdcn(ib).and.(ib.lt.knct3rd(ib))) then
            nssb2=nssb2+1
            j=knct3rd(ib)
            write(iun,'(a,i4,2a,i5,2a,i5)') 'SSBOND',nssb2,' CYS',
     +      chainid(ib),iseq(ib),'    CYS',chainid(j),iseq(j)
            write(iun,'(a,2i5)') 'CONECT',nserial(ib),nserial(j)
        endif
1100    continue        
        write(iun,'(a)')'END'
        return
        end
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE PRINTS THE CHAIN COORDINATES IN XYZ FORMAT
        subroutine print_conf_xyz(iun,time,energy,rms)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character*3 aseq,ares
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        write(iun,'(a,f10.1,4x,a,f10.4,4x,a,f8.4)')
     +  'TIME =',time,'ENERGY =',energy,'RMSD =',rms*unit
        iatom=0
        do 1000 ib=1,men
c        write(iun,'(3f12.5)')x0(ib)*unit,y0(ib)*unit,z0(ib)*unit
        write(iun,'(3f9.3)')x0(ib)*unit,y0(ib)*unit,z0(ib)*unit
1000    continue
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE PRINTS THE RESTART FILE
        subroutine print_restart(kb,itraj)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character rstfile*64,stafile*32,filname*32
        logical l3rdcn(len),lwal
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/restart/delta,work,sep0,rstfile,stafile,filname,klenstr
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
        common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        
        time=kb*delta
        write(stafile,*) '(a',klenstr,',i2.2,i10.10,a)'
        write(rstfile,stafile) filname,itraj,nint(time),'.rst'
        open(21,file=rstfile,status='unknown')
        write(21,'(2f18.3)')time,work
        write(21,'(8f14.6)')zup,zdown,yup,ydown,xup,xdown,sep0,shear
        write(21,'(7i10)')(kbwal(i),i=1,8)
        write(21,'(2i8)')intrsc,intesc
        do i=1,men
            if(l3rdcn(i)) write(21,*) i,knct3rd(i)
        enddo
        write(21,*)(ksorted(i),i=1,men)
        write(21,*)(xpul(i),ypul(i),zpul(i),i=1,men)
        write(21,*)(x0(i),y0(i),z0(i),i=1,men)
        write(21,*)(x1(i),y1(i),z1(i),i=1,men)
        write(21,*)icnss,icdss,ip1,ip2
        if(lwal) then
            write(21,*)(ipw(1,i),i=1,men)
            write(21,*)(ipw(2,i),i=1,men)
        endif
        write(21,*)kbperiod
        close(21)
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE APPLIES THE CONSTANT AFM FORCE TO THE PROTEIN
        subroutine afm(aforce,epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)

        ! STICK THE N-TERMINUS WITH A HARMONIC POTENTIAL
        dx = x0(ip1)-xn(ip1)
        dy = y0(ip1)-yn(ip1)
        dz = z0(ip1)-zn(ip1)
        rsq=dx*dx+dy*dy+dz*dz
        r = sqrt(rsq)
        ene1 = HH1*rsq + HH2*rsq*rsq
        epot = epot + ene1
        fce = (2*HH1+4*HH2*rsq)*r
        if(r.gt.1.d-10) then
        fce=-fce/r
        repx=fce*dx
        repy=fce*dy
        repz=fce*dz
        fx(1) = fx(1) + repx
        fy(1) = fy(1) + repy
        fz(1) = fz(1) + repz
        endif

        ! APPLY A FORCE TO THE DIRECTION PARRALLEL TO THE LINE
        ! CONNECTING FIRST AND LAST MONOMERS IN THE NATIVE STATE
        fax=aforce*afx
        fay=aforce*afy
        faz=aforce*afz
        fx(ip2)=fx(ip2)+fax
        fy(ip2)=fy(ip2)+fay
        fz(ip2)=fz(ip2)+faz

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C IMPLEMENTATION OF THE CONSTANT VELOCITY OF ONE END
        subroutine vafm(fresist,fresistperp,epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lwal,ldens,lshear,lobo,ljwal,lfcc,lpbcx,lpbcy,lpbcz
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/vel/x1(len),y1(len),z1(len)
        common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
        common/equil/kteql,ktrest,nratvel,nratveld,kconnecttime
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
        common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
        common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
        common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
        
        fresist=0.0
        fresistperp=0.0
        if(lwal) then ! MANY CHAINS
          if(lfcc) then
!$omp parallel default(shared), private(i,j,dx,dy,dz,
!$omp& rsq,r,rb,rb2,ene,fce,repx,repy,repz,rsi,rsig,r6,potpref)
!$omp do reduction(+:epot,fresist,fresistperp)
            do k=1,kfccw
                i=kfcc(1,k,jq)
                j=kfcc(2,k,jq)
                if(j.lt.men+menw/2) then
                    shearpull=shear
                else
                    shearpull=-shear
                endif
                dx = x0(i)-x0(j)-shearpull
                dy = y0(i)-y0(j)
                dz = z0(i)-z0(j)
                if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                if(lpbcy) dy = dy-ysep*nint(dy*yinv)
c               if(lpbcz) dz = dz-zsep*nint(dz*zinv)
                rsq=dx*dx+dy*dy+dz*dz
                r = sqrt(rsq)
                rsig=sigma1(9)
                if(r.le.rsig*cntfct) then
                    kfcc(3,k,jq)=min(abs(kfcc(3,k,jq))+1,int(ad))
                    if(j.lt.men+menw/2) then
                        icw(1)=icw(1)+1
                    else
                        icw(2)=icw(2)+1
                    endif
                else
                    if(kfcc(3,k,jq).lt.0) then
                        kfcc(3,k,jq)=kfcc(3,k,jq)+1
                    else
                        kfcc(3,k,jq)=min(0,1-kfcc(3,k,jq))
                    endif
                endif
                rsi=rsig/r
                r6=rsi**6
                if(r.le.walmindst) then
                    potpref=cofp*4.d0*r6
                else
                    potpref=cofp*abs(kfcc(3,k,jq))/ad*4.d0*r6
                endif
                ene=potpref*(r6-1.d0)
                fce=potpref*(6.d0-12.d0*r6)/r
                epot=epot+ene
                if(fce.gt.1.d+3) fce=1.d+3
                if(fce.lt.-1.d+3) fce=-1.d+3
                fce=-fce/r
                repx=fce*dx
                repy=fce*dy
                repz=fce*dz
                fx(i) = fx(i) + repx
c                fx(j) = fx(j) - repx
                fy(i) = fy(i) + repy
c                fy(j) = fy(j) - repy
                fz(i) = fz(i) + repz
c                fz(j) = fz(j) - repz
                    if(j.lt.men+menw/2) then
                        fresist=fresist+repz ! z-component of the force
                fresistperp=fresistperp+repx ! x-component
                    else
                        fresist=fresist-repz ! z-component of the force
                fresistperp=fresistperp-repx ! x-component
                    endif
            enddo
!$omp enddo nowait
!$omp end parallel
          else
           do ip=1,ip1+ip2
            if(ipw(1,ip).eq.0) goto 7865
            if(ipw(1,ip).lt.0) then
                zpull=zdown+zpul(ip)
                !zpull=zdown+walmindst
                shearpull=shear
                icw(1)=icw(1)+1
            endif
            if(ipw(1,ip).gt.0) then
                zpull=zup+zpul(ip)
                !zpull=zup-walmindst
                shearpull=-shear
                icw(2)=icw(2)+1
            endif
            i=ipw(2,ip)
            dx = x0(i)-xpul(ip)-shearpull
            dy = y0(i)-ypul(ip)
            dz = z0(i)-zpull
            rsq=dx*dx+dy*dy+dz*dz
            r = sqrt(rsq)
            if(ljwal) then
                rsi=sigma1(9)/r
                r6=rsi**6
                coefpul = cofp*z0temp(i)/ad
                ene1 = coefpul*4.d0*r6*(r6-1.d0)
                fce= coefpul*24.d0*r6*(1.d0-2.d0*r6)/r
                if(r.gt.sigma1(9)*1.5) then
                    if(z0temp(i).gt.0) then
                        z0temp(i)=z0temp(i)-1
                    else
                        z0temp(i)=0
                        ipw(1,ip)=0
                    endif
                else
                    if(z0temp(i).lt.ad) then
                        z0temp(i)=z0temp(i)+1
                    endif
                endif
            else
                ene1 = HH1*rsq + HH2*rsq*rsq
                fce = (2*HH1+4*HH2*rsq)*r
                if(lobo) then
                  if(z0temp(ip).lt.ad) z0temp(ip)=z0temp(ip)+1
                  ene1=ene1*cofp*z0temp(ip)/ad
                  fce=fce*cofp*z0temp(ip)/ad
                else
                  ene1=ene1*cofp*float(kbwal(kconnecttime))/ktrest
                  fce=fce*cofp*float(kbwal(kconnecttime))/ktrest
                endif
            endif
            if(r.gt.1.d-10) then
              fce=-fce/r
c              if(abs(fce).gt.1.d+3) write(*,*) 'ERR',r,z0(i),zdown,zup
              if(fce.gt.1.d+3) fce=1.d+3
              if(fce.lt.-1.d+3) fce=-1.d+3
              repx=fce*dx
              repy=fce*dy
              repz=fce*dz
              fx(i) = fx(i) + repx
              fy(i) = fy(i) + repy
              fz(i) = fz(i) + repz
              if(ipw(1,ip).gt.0) then
c              zuforce=zuforce-repz
                fresist=fresist-repz !z-component of the resisting force
                fresistperp=fresistperp-repx !x-component
              endif
              if(ipw(1,ip).lt.0) then
c              zdforce=zdforce+repz
                fresist=fresist+repz !z-component of the resisting force
                fresistperp=fresistperp+repx !x-component
              endif
            endif
            epot=epot+ene1
7865        continue
           enddo
          endif
        else ! ONLY ONE CHAIN
        ! STICK THE N-TERMINUS WITH A HARMONIC POTENTIAL
          dx = x0(ip1)-xn(ip1)
          dy = y0(ip1)-yn(ip1)
          dz = z0(ip1)-zn(ip1)
          rsq=dx*dx+dy*dy+dz*dz
          r = sqrt(rsq)
          ene1 = HH1*rsq + HH2*rsq*rsq
          fce = (2*HH1+4*HH2*rsq)*r
          ene1=ene1*cofp
          epot=epot+ene1
          fce=fce*cofp
          if(r.gt.1.d-10) then
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(ip1) = fx(ip1) + repx
            fy(ip1) = fy(ip1) + repy
            fz(ip1) = fz(ip1) + repz
          endif
        
          ! CONSTANT V IN THE DIRECTION PARRALLEL TO THE LINE
          ! CONNECTING FIRST AND LAST MONOMERS IN THE NATIVE STATE
          xpul(1)=xpul(1)+vpulx
          ypul(1)=ypul(1)+vpuly
          zpul(1)=zpul(1)+vpulz
          dx = x0(ip2)-xpul(1)
          dy = y0(ip2)-ypul(1)
          dz = z0(ip2)-zpul(1)
          rsq=dx*dx+dy*dy+dz*dz
          r = sqrt(rsq)
          ene1 = HH1*rsq + HH2*rsq*rsq
          fce = (2*HH1+4*HH2*rsq)*r
          ene1=ene1*cofp
          epot=epot+ene1
          fce=fce*cofp
          if(r.gt.1.d-10) then ! else fresist and fresistperp are 0
           fce=-fce/r
           repx=fce*dx
           repy=fce*dy
           repz=fce*dz
           fx(ip2) = fx(ip2) + repx
           fy(ip2) = fy(ip2) + repy
           fz(ip2) = fz(ip2) + repz
           fresist=repx*afx + repy*afy + repz*afz
c          perpendicular component of the resisting force, in xz plane
           fresistperp=repx*afz/(afx+afz**2/afx)-repz/(1.0+(afz/afx)**2)
          endif
        endif
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE CONNECTS THE DOMAINS FOR THE TITIN
        subroutine build_titin(ndomain)
        implicit double precision(a-h,o-z)

        parameter(len=10000)
        logical lconect(len)
        character aseq*3,pdbfile*32
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr

        if(men*ndomain.gt.len) then
        write(*,'(a)')'NUMBER OF DOMAIN EXCEEDS ALLOWED DIMENSION.'
        write(*,'(a)')'PROGRAM STOPPED.'
        stop
        endif

        rx=xn(men)-xn(1)
        ry=yn(men)-yn(1)
        rz=zn(men)-zn(1)
        rr=sqrt(rx*rx+ry*ry+rz*rz)
        rx=rx*(rr+bond)/rr
        ry=ry*(rr+bond)/rr
        rz=rz*(rr+bond)/rr
        men1=men
        
        ib=men
        do 1000 ido=2,ndomain
        do i=1,men1
        ib=ib+1
        xn(ib)=xn(i)+rx*(ido-1)
        yn(ib)=yn(i)+ry*(ido-1)
        zn(ib)=zn(i)+rz*(ido-1)
        iseq(ib)=iseq(i)+men1*(ido-1)
        aseq(ib)=aseq(i)
        enddo

        kk=(ido-1)*klont
        do k=1,klont
        i=klist(1,k)
        j=klist(2,k)
        l=klist(3,k)
        i1=(ido-1)*men1+i
        j1=(ido-1)*men1+j
        kk=kk+1
        klist(1,kk)=i1
        klist(2,kk)=j1
        klist(3,kk)=l
        enddo

1000        continue

        men=ib
        klont=kk

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE DISALLOWS CONTACTS BETWEEN DOMAINS
        subroutine interdomain(ndomain)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        dimension klst(3,len*20)
        
        idsize=men/ndomain
        kn=0
        do 1000 k=1,klont
        i=klist(1,k)
        j=klist(2,k)
        l=klist(3,k)
        if((i-1)/idsize.ne.(j-1)/idsize) then
        klist(3,k)=-1 ! inter-chain contacts should have -1 here
        else
        kn=kn+1
        klst(1,kn)=i
        klst(2,kn)=j
        klst(3,kn)=l
        endif
1000        continue
        klont=kn
        do k=1,klont
        klist(1,k)=klst(1,k)
        klist(2,k)=klst(2,k)
        klist(3,k)=klst(3,k)
        enddo
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE READS THE MASSES OF AMINO ACIDS
        subroutine amino_acid_mass
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,ares*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/mass/rmas(len),rsqmas(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        do 1000 ib=1,men
        ares=aseq(ib)
        if(ares.eq.'GLY') then
        rmas(ib)= 57.052
        else if(ares.eq.'ALA') then
        rmas(ib)= 71.079
        else if(ares.eq.'SER') then
        rmas(ib)= 87.078
        else if(ares.eq.'PRO') then
        rmas(ib)= 97.117
        else if(ares.eq.'VAL') then
        rmas(ib)= 99.133
        else if(ares.eq.'THR') then
        rmas(ib)=101.105
        else if(ares.eq.'CYS') then
        rmas(ib)=103.144
        else if(ares.eq.'ILE') then
        rmas(ib)=113.160
        else if(ares.eq.'LEU') then
        rmas(ib)=113.160
        else if(ares.eq.'ASN') then
        rmas(ib)=114.104
        else if(ares.eq.'ASP') then
        rmas(ib)=115.089
        else if(ares.eq.'GLN') then
        rmas(ib)=128.131
        else if(ares.eq.'LYS') then
        rmas(ib)=128.174
        else if(ares.eq.'GLU') then
        rmas(ib)=129.116
        else if(ares.eq.'MET') then
        rmas(ib)=131.198
        else if(ares.eq.'HIS') then
        rmas(ib)=137.142
        else if(ares.eq.'PHE') then
        rmas(ib)=147.177
        else if(ares.eq.'ARG') then
        rmas(ib)=156.188
        else if(ares.eq.'TYR') then
        rmas(ib)=163.170
        else if(ares.eq.'TRP') then
        rmas(ib)=186.213
        endif
1000        continue

        ! AVERAGE MASS
        avmas=0.d0
        do ib=1,men
        avmas=avmas+rmas(ib)
        enddo
        avmas=avmas/men

        ! REDUCED MASS AND ITS SQUARE ROOT
        do ib=1,men
        rma=rmas(ib)/avmas
        rmas(ib)=rma
        rsqmas(ib)=sqrt(rma)
        enddo

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine model_chirality
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lfrompdb(len),lconect(len),lchiral,langle,lsimpang
        character aseq*3
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
        dimension xs(len),ys(len),zs(len)

        do 1000 ib=1,men-1
        ib1=ib+1
        xs(ib)=xn(ib1)-xn(ib)
        ys(ib)=yn(ib1)-yn(ib)
        zs(ib)=zn(ib1)-zn(ib)
1000    continue

        do 2000 ib=1,men-3
        if(lconect(ib) .and. lconect(ib+1)) then
        ib1=ib+1
        ib2=ib+2
        xa=ys(ib)*zs(ib1)-zs(ib)*ys(ib1)
        ya=zs(ib)*xs(ib1)-xs(ib)*zs(ib1)
        za=xs(ib)*ys(ib1)-ys(ib)*xs(ib1)
        as=xa*xs(ib2)+ya*ys(ib2)+za*zs(ib2)
        aa=b(ib)*b(ib1)*b(ib2)
c        aa=bond**3
        chirn(ib)=as/aa
        if(lsimpang) chirn(ib)=sigma1(min(inameseq(ib),3))
        endif
2000    continue
c        do ib=1,men-3
c        write(1,'(i6,f12.4)')ib,chirn(ib)
c        enddo
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE COMPUTES THE ENERGY AND FORCE GIVEN BY THE
C CHIRALITY POTENTIALS
        subroutine eval_chirality(enechi)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lchiral,langle
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
        common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        dimension xs(len),ys(len),zs(len)

        do 1000 ib=1,men-1
        ib1=ib+1
        xs(ib) = v(1,ib) !=x0(ib1)-x0(ib)
        ys(ib) = v(2,ib) !=y0(ib1)-y0(ib)
        zs(ib) = v(3,ib) !=z0(ib1)-z0(ib)
1000    continue

        do 2000 ib=1,men-3
        if(lconect(ib) .and. lconect(ib+1)) then
        ib1=ib+1
        ib2=ib+2
        xa = -vxv(1,ib1) !=ys(ib)*zs(ib1)-zs(ib)*ys(ib1)
        ya = -vxv(2,ib1) !=zs(ib)*xs(ib1)-xs(ib)*zs(ib1)
        za = -vxv(3,ib1) !=xs(ib)*ys(ib1)-ys(ib)*xs(ib1)
        aa=xa*xs(ib2)+ya*ys(ib2)+za*zs(ib2)
        chir(ib)=aa/b(ib2) !**3
        endif
2000    continue

        enechi=0.d0
        do 3000 ib=1,men-3
        if(lconect(ib) .and. lconect(ib+1)) then
        chira=chir(ib)
        if(chira*chirn(ib).gt.0.d0) goto 3000
        ib1=ib+1
        ib2=ib+2
        ib3=ib+3
        ! COMPUTE THE ENERGY
        enechi=enechi+echi*chira*chira/2.d0
        ! COMPUTE THE FORCE
        fchi=-echi*chira/bond**3
        ! FIRST RESIDUE
        repx=ys(ib2)*zs(ib1)-zs(ib2)*ys(ib1)
        repy=zs(ib2)*xs(ib1)-xs(ib2)*zs(ib1)
        repz=xs(ib2)*ys(ib1)-ys(ib2)*xs(ib1)
        fx(ib)=fx(ib)+fchi*repx
        fy(ib)=fy(ib)+fchi*repy
        fz(ib)=fz(ib)+fchi*repz
        ! SECOND RESIDUE
        xa=xs(ib)+xs(ib1)
        ya=ys(ib)+ys(ib1)
        za=zs(ib)+zs(ib1)
        repx=ya*zs(ib2)-za*ys(ib2)
        repy=za*xs(ib2)-xa*zs(ib2)
        repz=xa*ys(ib2)-ya*xs(ib2)
        fx(ib1)=fx(ib1)+fchi*repx
        fy(ib1)=fy(ib1)+fchi*repy
        fz(ib1)=fz(ib1)+fchi*repz
        ! THIRD RESIDUE
        xa=xs(ib1)+xs(ib2)
        ya=ys(ib1)+ys(ib2)
        za=zs(ib1)+zs(ib2)
        repx=ya*zs(ib)-za*ys(ib)
        repy=za*xs(ib)-xa*zs(ib)
        repz=xa*ys(ib)-ya*xs(ib)
        fx(ib2)=fx(ib2)+fchi*repx
        fy(ib2)=fy(ib2)+fchi*repy
        fz(ib2)=fz(ib2)+fchi*repz
        ! FOURTH RESIDUE
        repx=ys(ib)*zs(ib1)-zs(ib)*ys(ib1)
        repy=zs(ib)*xs(ib1)-xs(ib)*zs(ib1)
        repz=xs(ib)*ys(ib1)-ys(ib)*xs(ib1)
        fx(ib3)=fx(ib3)+fchi*repx
        fy(ib3)=fy(ib3)+fchi*repy
        fz(ib3)=fz(ib3)+fchi*repz
        endif
3000    continue

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE LOADS THE NATIVE BACKBONE CONFORMATION OF A PROTEIN 
C FROM ITS PDB FILE AND THEN CALCULATES THE TORSION ANGLES OF THE
C NATIVE STATE
        subroutine load_allatom(filename)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,ares*3,filename*32,bb*3,buffer*128,bbo*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        character aname*3
        common/rall/rx_a(len,14),ry_a(len,14),rz_a(len,14)
        common/nall/na(len),aname(len,14)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        open(15,file=filename,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,2x,a)')'ERROR OPENING FILENAME',filename
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        ivalold=-1
        ib=0
        do i=1,len
        na(i)=0
        enddo
        bbo='   '
15        read(15,'(a)',end=20) buffer
        if(buffer(1:3).eq.'END') goto 20
        if(buffer(1:4).ne.'ATOM') goto 15
        read(buffer,'(13x,a3,1x,a3,2x,i4,4x,3f8.3)')
     +  bb,ares,ival,xval,yval,zval
        if(bb.eq.bbo) goto 15
        if((buffer(17:17).ne.'A').and.(buffer(17:17).ne.' ')) goto 15
        bbo=bb
        if(ival.ne.ivalold) then !if(bb.eq.'N  ')then
        ivalold=ival
        ib=ib+1
        ja=0
        aseq(ib)=ares
        iseq(ib)=ival
        endif
        if(bb(1:1).eq.'N'.or.bb(1:1).eq.'C'.or.bb(1:1).eq.'O'.or.
     +     bb(1:1).eq.'S')then
        ja=ja+1
        rx_a(ib,ja)=xval
        ry_a(ib,ja)=yval
        rz_a(ib,ja)=zval
        na(ib)=ja
        aname(ib,ja)=bb
        endif
        goto 15
20        continue
        close(15)
        men=ib
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine assign_VdW_radius
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,ares*3,ana*3,filename*32,bb*3,buffer*128
        common/sequence/iseq(len),inameseq(len),aseq(len)
        character aname*3
        common/rall/rx_a(len,14),ry_a(len,14),rz_a(len,14)
        common/nall/na(len),aname(len,14)
        common/radi/vrad(len,14)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        dimension nb(len)

        do 1000 ib=1,men
        vrad(ib,1)=1.64                ! BACKBONE NITROGEN
        vrad(ib,2)=1.88                ! BACKBONE C-ALPHA
        vrad(ib,3)=1.61                ! BACKBONE C'
        vrad(ib,4)=1.42                ! BACKBONE OXYGEN
        ares=aseq(ib)
        if(ares.eq.'GLY') then
          nb(ib)=4
          goto 1000
        else if(ares.eq.'PRO') then
          nb(ib)=7
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
        else if(ares.eq.'GLN') then
          nb(ib)=9
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.61
          vrad(ib,8)=1.42
          vrad(ib,9)=1.64
        else if(ares.eq.'CYS') then
          nb(ib)=6
          vrad(ib,5)=1.88
          vrad(ib,6)=1.77
        else if(ares.eq.'VAL') then
          nb(ib)=7
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
        else if(ares.eq.'PHE') then
          nb(ib)=11
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.61
          vrad(ib,8)=1.76
          vrad(ib,9)=1.76
          vrad(ib,10)=1.76
          vrad(ib,11)=1.76
        else if(ares.eq.'MET') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.77
          vrad(ib,8)=1.88
        else if(ares.eq.'ILE') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.88
        else if(ares.eq.'ASP') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.46
          vrad(ib,8)=1.42
        else if(ares.eq.'GLU') then
          nb(ib)=9
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.61
          vrad(ib,8)=1.46
          vrad(ib,9)=1.42
        else if(ares.eq.'LYS') then
          nb(ib)=9
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.88
          vrad(ib,9)=1.64
        else if(ares.eq.'ARG') then
          nb(ib)=11
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.64
          vrad(ib,9)=1.61
          vrad(ib,10)=1.64
          vrad(ib,11)=1.64
        else if(ares.eq.'SER') then
          nb(ib)=6
          vrad(ib,5)=1.88
          vrad(ib,6)=1.46
        else if(ares.eq.'THR') then
          nb(ib)=7
          vrad(ib,5)=1.88
          vrad(ib,6)=1.46
          vrad(ib,7)=1.88
        else if(ares.eq.'TYR') then
          nb(ib)=12
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.76
          vrad(ib,8)=1.76
          vrad(ib,9)=1.76
          vrad(ib,10)=1.76
          vrad(ib,11)=1.61
          vrad(ib,12)=1.46
        else if(ares.eq.'HIS') then
          nb(ib)=10
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.64
          vrad(ib,8)=1.76
          vrad(ib,9)=1.76
          vrad(ib,10)=1.64
        else if(ares.eq.'ASN') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.42
          vrad(ib,8)=1.64
        else if(ares.eq.'TRP') then
          nb(ib)=14
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.76
          vrad(ib,8)=1.61
          vrad(ib,9)=1.64
          vrad(ib,10)=1.61
          vrad(ib,11)=1.76
          vrad(ib,12)=1.76
          vrad(ib,13)=1.76
          vrad(ib,14)=1.76
        else if(ares.eq.'ALA') then
          nb(ib)=5
          vrad(ib,5)=1.88
        else if(ares.eq.'LEU') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.88
        endif
1000        continue

        ! CHECK AND CORRECTION
        do ib=1,men
        nat=na(ib)
        if(nat.lt.nb(ib)) then
        write(1,*)'AMINO ACID HAS FEWER ATOMS THAN SHOULD BE'
        write(1,'(i5,2x,a3,2x,2i6)')iseq(ib),aseq(ib),nat,nb(ib)
        else if(nat.gt.nb(ib)) then
c       write(1,*)'AMINO ACID HAS MORE ATOMS THAN SHOULD BE'
c       write(1,'(i5,2x,a3,2x,2i6)')iseq(ib),aseq(ib),nat,nb(ib)
c       do j=1,nat
c       write(1,'(a3,2x,i3,2x,a3)')aseq(ib),iseq(ib),aname(ib,j)
c       enddo
        nat=nb(ib)
        na(ib)=nb(ib)
        endif
        do j=1,nat
        ares=aseq(ib)
        ana=aname(ib,j)
        rad=vrad(ib,j)
        if(ana(1:1).eq.'N'.and.rad.ne.1.64) then
        vrad(ib,j)=1.64
        else if(ana(1:1).eq.'S'.and.rad.ne.1.77) then
        vrad(ib,j)=1.77
        else if(ana(1:1).eq.'O'.and.rad.ne.1.42.and.rad.ne.1.46) then
        vrad(ib,j)=1.46
        else if(ana(1:1).eq.'C'.and.rad.ne.1.88.and.rad.ne.1.76.
     +          and.rad.ne.1.61) then
        vrad(ib,j)=1.88
        endif
        if(rad.eq.0.) then
         write(1,'(a)')'ATOM ERROR:'
         write(1,'(a3,2x,i3,2x,a3)')aseq(ib),iseq(ib),aname(ib,j)
         stop
        endif
        enddo
        enddo
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine gopotential(asigma)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lpbcx,lpbcy,lpbcz
        common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/vel/x1(len),y1(len),z1(len)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        dimension histo(len)

          icn=0
          nhist=0
          delh=0.5
          mah=0
          do i=1,200
          histo(i)=0
          enddo
          bmin=1000
          bmax=-1000

        dsum1=0.d0
        do 2000 k=1,klont
        i=klist(1,k)
        j=klist(2,k)
        dx=xn(i)-xn(j)
        dy=yn(i)-yn(j)
        dz=zn(i)-zn(j)
        if(lpbcx) dx = dx-xsep*nint(dx*xinv)
        if(lpbcy) dy = dy-ysep*nint(dy*yinv)
        if(lpbcz) dz = dz-zsep*nint(dz*zinv)
        dal=dx*dx+dy*dy+dz*dz
        dal=dsqrt(dal)
        dsum1=dsum1+dal
        sont(k)=dal*0.5d0**(1.d0/6.d0) ! dal for 1012
c        write(1,'(3i6,f12.3)')k,i,j,dal*unit
           histl=dal*unit
           if(histl.gt.bmax) bmax=histl
           if(histl.lt.bmin) bmin=histl
           mhi=nint(histl/delh)
           if(mhi.eq.0) mhi=1
           if(mhi.gt.mah) mah=mhi
           nhist=nhist+1
           histo(mhi)=histo(mhi)+1
2000    continue

c           write(1,*)'contacts between ',bmin,'  and ',bmax
c           write(1,*)'histogram'
c           do 144 i=1,mah
c           hh=histo(i)
c           if(hh.eq.0) go to 144
c           sep=i*delh
c           hh=hh/nhist
c           write(1,*)sep,hh
c144        continue
        if(klont.gt.0) then
            asigma=dsum1/klont
        else
            asigma=0.d0
        endif

c       attractive non-native contacts
c       sigma corresponds to 5 A

c        arut=5.d0/unit
c        do i=1,kront
c        if(krist(3,i).eq.2) then
c        sigma1(2)=arut*0.5d0**(1.d0/6.d0)
c        endif
c        enddo
c        enddo
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine compute_native_angles
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lfrompdb(len),lsimpang
        common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
        common/four/rx(4),ry(4),rz(4)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang

        do 100 ib=2,men-1
        do ibb=1,3
        ib1=ib+ibb-2
        rx(ibb)=xn(ib1)
        ry(ibb)=yn(ib1)
        rz(ibb)=zn(ib1)
        enddo
        call bondangle(theta)
        the0(ib)=theta
100        continue

        do 200 ib=3,men-1
        do ibb=1,4
        ib1=ib+ibb-3
        rx(ibb)=xn(ib1)
        ry(ibb)=yn(ib1)
        rz(ibb)=zn(ib1)
        enddo
        call dihedral(phi)
        phi0(ib)=phi
200        continue

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine angeval(epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lfrompdb(len),lsimpang,langle,ldi,ldisimp
        character*3 aseq,ares
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
        common/four/rx(4),ry(4),rz(4)
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        ! CALCULATE ENERGY AND FORCE RELATED TO THE BOND ANGLE
        do 100 ib=2,men-1
        if(lconect(ib-1).and.lconect(ib)) then
        do ibb=1,3
        ib1=ib+ibb-2
        rx(ibb)=x0(ib1)
        ry(ibb)=y0(ib1)
        rz(ibb)=z0(ib1)
        enddo
        call bondangle(theta)
        thetemp(ib)=theta
        
        if(langle) then
        if(lfrompdb(ib-1).and.lfrompdb(ib+1)) then
            theta=theta-the0(ib)
            ene=CBA*theta*theta
        else
            j=min(inameseq(ib),3)
            k=min(inameseq(ib+1),3)
c            if(j.eq.2) theta=min(theta,3.03)
            theta2=theta*theta
            theta3=theta2*theta
            ene=angpot(2,j,k)*theta+angpot(3,j,k)*theta2+
     +      angpot(4,j,k)*theta3+angpot(5,j,k)*theta2**2
     +      +angpot(6,j,k)*theta2*theta3+angpot(7,j,k)*theta3**2
        endif
        epot=epot+ene

        do 101 ibb=1,3

        ib1=ib+ibb-2

        rxo=rx(ibb)
        ryo=ry(ibb)
        rzo=rz(ibb)

        rx(ibb)=rxo+dar
        call bondangle(theta)
        if(lfrompdb(ib-1).and.lfrompdb(ib+1)) then
            theta=theta-the0(ib)
            ene1=CBA*theta*theta
        else
c            if(j.eq.2) theta=min(theta,3.03)
            theta2=theta*theta
            theta3=theta2*theta
            ene1=angpot(2,j,k)*theta+angpot(3,j,k)*theta2+
     +      angpot(4,j,k)*theta3+angpot(5,j,k)*theta2**2
     +      +angpot(6,j,k)*theta2*theta3+angpot(7,j,k)*theta3**2
        endif
        fce=(ene1-ene)/dar
        if(abs(fce).gt.1.d+3) fce=sign(1.d+3,fce)
        fx(ib1)=fx(ib1)-fce
        rx(ibb)=rxo

        ry(ibb)=ryo+dar
        call bondangle(theta)
        if(lfrompdb(ib-1).and.lfrompdb(ib+1)) then
            theta=theta-the0(ib)
            ene1=CBA*theta*theta
        else
c            if(j.eq.2) theta=min(theta,3.03)
            theta2=theta*theta
            theta3=theta2*theta
            ene1=angpot(2,j,k)*theta+angpot(3,j,k)*theta2+
     +      angpot(4,j,k)*theta3+angpot(5,j,k)*theta2**2
     +      +angpot(6,j,k)*theta2*theta3+angpot(7,j,k)*theta3**2
        endif
        fce=(ene1-ene)/dar
        if(abs(fce).gt.1.d+3) fce=sign(1.d+3,fce)
        fy(ib1)=fy(ib1)-fce
        ry(ibb)=ryo

        rz(ibb)=rzo+dar
        call bondangle(theta)
        if(lfrompdb(ib-1).and.lfrompdb(ib+1)) then
            theta=theta-the0(ib)
            ene1=CBA*theta*theta
        else
c            if(j.eq.2) theta=min(theta,3.03)
            theta2=theta*theta
            theta3=theta2*theta
            ene1=angpot(2,j,k)*theta+angpot(3,j,k)*theta2+
     +      angpot(4,j,k)*theta3+angpot(5,j,k)*theta2**2
     +      +angpot(6,j,k)*theta2*theta3+angpot(7,j,k)*theta3**2
        endif
        fce=(ene1-ene)/dar
        if(abs(fce).gt.1.d+3) fce=sign(1.d+3,fce)
        fz(ib1)=fz(ib1)-fce
        rz(ibb)=rzo

101        continue
        endif
        endif
100        continue


        ! CALCULATE ENERGY AND FORCE RELATED TO THE DIHEDRAL ANGLE
        do 200 ib=3,men-1
        if(lconect(ib-2) .and. lconect(ib-1) .and. lconect(ib)) then
        do ibb=1,4
        ib1=ib+ibb-3
        rx(ibb)=x0(ib1)
        ry(ibb)=y0(ib1)
        rz(ibb)=z0(ib1)
        enddo

        call dihedral(phi)
        phitemp(ib)=phi
        if(langle.and.ldi) then
        if(lfrompdb(ib-2).and.lfrompdb(ib+1)) then
            phi=phi-phi0(ib)
            ene=CDA*(1.d0-dcos(phi))+CDB*(1.d0-dcos(3.d0*phi))
        else
            j=min(inameseq(ib-1),3)
            k=min(inameseq(ib),3)
            sinfi=dsin(phi)
            cosfi=dcos(phi)
            ene=dihpot(2,j,k)*sinfi+dihpot(3,j,k)*cosfi
     +      +dihpot(4,j,k)*sinfi**2+dihpot(5,j,k)*cosfi**2
     +      +dihpot(6,j,k)*sinfi*cosfi
        endif
        epot=epot+ene    

        do 201 ibb=1,4

        ib1=ib+ibb-3

        rxo=rx(ibb)
        ryo=ry(ibb)
        rzo=rz(ibb)

        rx(ibb)=rxo+dar
        call dihedral(phi)
        if(lfrompdb(ib-2).and.lfrompdb(ib+1)) then
            phi=phi-phi0(ib)
            ene1=CDA*(1.d0-dcos(phi))+CDB*(1.d0-dcos(3.d0*phi))
        else
            sinfi=dsin(phi)
            cosfi=dcos(phi)
            ene1=dihpot(2,j,k)*sinfi+dihpot(3,j,k)*cosfi
     +      +dihpot(4,j,k)*sinfi**2+dihpot(5,j,k)*cosfi**2
     +      +dihpot(6,j,k)*sinfi*cosfi
        endif
        fx(ib1)=fx(ib1)-(ene1-ene)/dar
        rx(ibb)=rxo

        ry(ibb)=ryo+dar
        call dihedral(phi)
        if(lfrompdb(ib-2).and.lfrompdb(ib+1)) then
            phi=phi-phi0(ib)
            ene1=CDA*(1.d0-dcos(phi))+CDB*(1.d0-dcos(3.d0*phi))
        else
            sinfi=dsin(phi)
            cosfi=dcos(phi)
            ene1=dihpot(2,j,k)*sinfi+dihpot(3,j,k)*cosfi
     +      +dihpot(4,j,k)*sinfi**2+dihpot(5,j,k)*cosfi**2
     +      +dihpot(6,j,k)*sinfi*cosfi
        endif
        fy(ib1)=fy(ib1)-(ene1-ene)/dar
        ry(ibb)=ryo

        rz(ibb)=rzo+dar
        call dihedral(phi)
        if(lfrompdb(ib-2).and.lfrompdb(ib+1)) then
            phi=phi-phi0(ib)
            ene1=CDA*(1.d0-dcos(phi))+CDB*(1.d0-dcos(3.d0*phi))
        else
            sinfi=dsin(phi)
            cosfi=dcos(phi)
            ene1=dihpot(2,j,k)*sinfi+dihpot(3,j,k)*cosfi
     +      +dihpot(4,j,k)*sinfi**2+dihpot(5,j,k)*cosfi**2
     +      +dihpot(6,j,k)*sinfi*cosfi
        endif
        fz(ib1)=fz(ib1)-(ene1-ene)/dar
        rz(ibb)=rzo

201     continue
        endif
        endif
200     continue
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine countpid(i,j,phi)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lene
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        dimension fi(3),fj(3),fk(3),fl(3)
        dimension rm(3),rn(3),rij(3),rkj(3),rkl(3)

        ! CALCULATE ENERGY AND FORCE RELATED TO THE DIHEDRAL ANGLE
        phi=99.d0 ! IF NO ERROR, THIS SHOULD BE REWRITTEN
        if(lconect(i-1) .and. lconect(i)) then
        ! COMPUTING FORCE (THE SAME FOR EVERY TYPE OF POTENTIAL)
C     B EKKER, H., BERENDSEN, H. J. C. AND VAN GUNSTEREN, W. F. (1995)
C     FORCE AND VIRIAL OF TORSIONAL-ANGLE-DEPENDENT POTENTIALS
C     J. COMPUT. CHEM., 16: 527-533. DOI: 10.1002/JCC.540160502
        ! v(i1)=-rij, v(i2)=rkj, v(i3)=-rkl
        i1=i
        i2=i+1
        i3=i-1
        i4=j
        rij(1)=x0(i1)-x0(i2)
        rij(2)=y0(i1)-y0(i2)
        rij(3)=z0(i1)-z0(i2)
        rkj(1)=x0(i3)-x0(i2)
        rkj(2)=y0(i3)-y0(i2)
        rkj(3)=z0(i3)-z0(i2)
        rkl(1)=x0(i3)-x0(i4)
        rkl(2)=y0(i3)-y0(i4)
        rkl(3)=z0(i3)-z0(i4)
        rm(1)=rij(2)*rkj(3)-rij(3)*rkj(2)
        rm(2)=rij(3)*rkj(1)-rij(1)*rkj(3)
        rm(3)=rij(1)*rkj(2)-rij(2)*rkj(1)
        rn(1)=rkj(2)*rkl(3)-rkj(3)*rkl(2)
        rn(2)=rkj(3)*rkl(1)-rkj(1)*rkl(3)
        rn(3)=rkj(1)*rkl(2)-rkj(2)*rkl(1)
        d2n=0.d0
        d2m=0.d0
        d2rkj=0.d0
        do k=1,3
           d2n=d2n+rn(k)*rn(k)
           d2m=d2m+rm(k)*rm(k)
           d2rkj=d2rkj+rkj(k)*rkj(k)
        enddo
        if (abs(d2n).lt.0.01d0 .or. abs(d2m).lt.0.01d0) then
           do k=1,3
            fi(k)=0.d0
            fj(k)=0.d0
            fk(k)=0.d0
            fl(k)=0.d0
           enddo
        else
           dmn=0.d0
           do k=1,3
             dmn=dmn+rn(k)*rm(k)
           enddo
           cosphi=dmn/sqrt(d2n*d2m)
           cosphi=min(cosphi,1.d0)
           cosphi=max(cosphi,-1.d0)
           rijn=0.d0
           do k=1,3
            rijn=rijn+rij(k)*rn(k)
           enddo
           if (rijn.lt.0.d0) then
            phi=-1.d0*acos(cosphi)
           else
            phi=acos(cosphi)
           endif
        endif
        endif
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine evalangles(epot,lsldh,cofdih)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lfrompdb(len),lsimpang,langle,ldi,ldisimp
        logical lenetab,lcoildih,lsldh
        character*3 aseq,ares
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
        common/four/rx(4),ry(4),rz(4)
        common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
        common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
        common/mr/tene(18001,3,3,2),pp,lenetab 
        dimension fi(3),fj(3),fk(3),fl(3),rp(3),qij(3),qkj(3),qp(3)
        dimension ra(3),rb(3),rm(3),rn(3)
        ! CALCULATE ENERGY AND FORCE RELATED TO THE BOND ANGLE
!$omp parallel default(shared), private(ib,i1,i2,i3,rijrkj,
!$omp& d2ij,d2kj,k,dij,dkj,costh,theta,fi,fj,fk,rp,drp,qij,
!$omp& qkj,qp,ra,rb,ene,dvdp,j,theta2,theta3,theta4,theta5)
!$omp do reduction(+:epot)
        do 100 ib=2,men-1
        if(lconect(ib-1).and.lconect(ib)) then
        ! COMPUTING FORCE (THE SAME FOR EVERY TYPE OF POTENTIAL)
C      SWOPE, W. C. AND FERGUSON, D. M. (1992) ALTERNATIVE EXPRESSIONS 
C      FOR ENERGIES AND FORCES DUE TO ANGLE BENDING AND TORSIONAL ENERGY
C      J. COMPUT. CHEM. 13: 585-594. DOI: 10.1002/JCC.540130508
        ! v(i1)=-rij, v(i2)=rkj
        i1=ib-1
        i2=ib
        i3=ib+1
        rijrkj=0.d0
        d2ij=0.d0
        d2kj=0.d0
        do k=1,3
         rijrkj=rijrkj-v(k,i1)*v(k,i2)
         d2ij=d2ij+v(k,i1)*v(k,i1)
         d2kj=d2kj+v(k,i2)*v(k,i2)
        enddo
        dij=sqrt(d2ij)
        dkj=sqrt(d2kj)
        costh=rijrkj/(dij*dkj)
        costh=min(costh,1.d0)
        costh=max(costh,-1.d0)
        theta=acos(costh)
        if (costh.eq.0.d0) then
         do k=1,3
          fi(k)=0.d0
          fj(k)=0.d0
          fk(k)=0.d0
         enddo
        else
         do k=1,3
          qij(k)=-v(k,i1)/dij
          qkj(k)=v(k,i2)/dkj
          qp(k)=vxv(k,i2)
         enddo
         ra(1)=qij(2)*qp(3)-qij(3)*qp(2)
         ra(2)=qij(3)*qp(1)-qij(1)*qp(3)
         ra(3)=qij(1)*qp(2)-qij(2)*qp(1)
         rb(1)=qkj(2)*qp(3)-qkj(3)*qp(2)
         rb(2)=qkj(3)*qp(1)-qkj(1)*qp(3)
         rb(3)=qkj(1)*qp(2)-qkj(2)*qp(1)
         do k=1,3
          fi(k)=ra(k)/dij
          fk(k)=-rb(k)/dkj
          fj(k)=-fi(k)-fk(k)
         enddo
        endif
        ! COMPUTING POTENTIAL-DEPENDENT PART
        thetemp(ib)=theta
        if(langle) then
        if(lfrompdb(ib-1).and.lfrompdb(ib+1)) then
            theta=theta-the0(ib)
            ene=CBA*theta*theta
            dvdp=2.d0*CBA*theta
        else
            j=min(inameseq(ib),3)
            k=min(inameseq(ib+1),3)
            if (lenetab) then ! tabularized potential
              xx=theta/pp ! xx - theta in 0.01 degree units
              ixx=INT(xx)  ! ixx - number of the table row
              xi=xx-ixx   ! xi - fractional part to interpolate
              aa=tene(ixx,j,k,1) ! 1-energy, 2-force
              ene=aa+xi*(tene(ixx+1,j,k,1)-aa)
              aa=tene(ixx,j,k,2)
              dvdp=(aa+xi*(tene(ixx+1,j,k,2)-aa))/pp/(-660)
            else ! polynomial potential
              theta2=theta*theta
              theta3=theta2*theta
              theta4=theta3*theta
              theta5=theta4*theta
              ene=angpot(1,j,k)+angpot(2,j,k)*theta+angpot(3,j,k)*theta2
     +        +angpot(4,j,k)*theta3+angpot(5,j,k)*theta4
     +        +angpot(6,j,k)*theta5+angpot(7,j,k)*theta3**2
              dvdp=angpot(9,j,k)+angpot(10,j,k)*theta
     +        +angpot(11,j,k)*theta2+angpot(12,j,k)*theta3
     +        +angpot(13,j,k)*theta4+angpot(14,j,k)*theta5
            endif
        endif
        epot=epot+ene
        fx(i1)=fx(i1)-dvdp*fi(1)
        fy(i1)=fy(i1)-dvdp*fi(2)
        fz(i1)=fz(i1)-dvdp*fi(3)
        fx(i2)=fx(i2)-dvdp*fj(1)
        fy(i2)=fy(i2)-dvdp*fj(2)
        fz(i2)=fz(i2)-dvdp*fj(3)
        fx(i3)=fx(i3)-dvdp*fk(1)
        fy(i3)=fy(i3)-dvdp*fk(2)
        fz(i3)=fz(i3)-dvdp*fk(3)
        endif
        endif
100     continue
!$omp enddo nowait
!$omp end parallel 
        ! CALCULATE ENERGY AND FORCE RELATED TO THE DIHEDRAL ANGLE
!$omp parallel default(shared), private(ib,i1,i2,i3,i4,rm,rn,d2n,d2m,
!$omp& d2rkj,k,phi,fi,fj,fk,fl,dmn,cosphi,rijn,drkj,rijrkj,rklrkj,df,
!$omp& ene,dvdp,j,sinfi,cosfi,sin2fi,cos2fi,sincosfi)
!$omp do reduction(+:epot)
        do 200 ib=3,men-1
        if(lconect(ib-2) .and. lconect(ib-1) .and. lconect(ib)) then
        ! COMPUTING FORCE (THE SAME FOR EVERY TYPE OF POTENTIAL)
C     B EKKER, H., BERENDSEN, H. J. C. AND VAN GUNSTEREN, W. F. (1995)
C     FORCE AND VIRIAL OF TORSIONAL-ANGLE-DEPENDENT POTENTIALS
C     J. COMPUT. CHEM., 16: 527-533. DOI: 10.1002/JCC.540160502
        ! v(i1)=-rij, v(i2)=rkj, v(i3)=-rkl
        i1=ib-2
        i2=ib-1
        i3=ib
        i4=ib+1
        rm(1)=vxv(1,i2)
        rm(2)=vxv(2,i2)
        rm(3)=vxv(3,i2)
        rmnrm=vnrm(i2)
        rn(1)=vxv(1,i3)
        rn(2)=vxv(2,i3)
        rn(3)=vxv(3,i3)
        rnnrm=vnrm(i3)
        !d2n=0.d0
        !d2m=0.d0
        d2rkj=0.d0
        do k=1,3
        ! d2n=d2n+rn(k)*rn(k)
        ! d2m=d2m+rm(k)*rm(k)
         d2rkj=d2rkj+v(k,i2)*v(k,i2)
        enddo
!        if (d2n.eq.0.d0 .or. d2m.eq.0.d0) then
!         phi=0.d0
!         do k=1,3
!          fi(k)=0.d0
!          fj(k)=0.d0
!          fk(k)=0.d0
!          fl(k)=0.d0
!         enddo
!       else
         dmn=rn(1)*rm(1)+rn(2)*rm(2)+rn(3)*rm(3)
         cosphi=max(min(dmn,1.d0),-1.d0) !/sqrt(d2n*d2m)
         rijn=0.d0
         do k=1,3
          rijn=rijn-v(k,i1)*rn(k)
         enddo
         if (rijn.lt.0.d0) then
          phi=-1.d0*acos(cosphi)
         else
          phi=acos(cosphi)
         endif
         drkj=sqrt(d2rkj)
         do k=1,3
          fi(k)=rm(k)*drkj/rmnrm !/d2m
          fl(k)=-rn(k)*drkj/rnnrm !/d2n
         enddo
         rijrkj=0.d0
         rklrkj=0.d0
         do k=1,3
          rijrkj=rijrkj-v(k,i1)*v(k,i2)
          rklrkj=rklrkj-v(k,i3)*v(k,i2)
         enddo
         do k=1,3
          df=(fi(k)*rijrkj-fl(k)*rklrkj)/d2rkj
          fj(k)=-fi(k)+df
          fk(k)=-fl(k)-df
         enddo
!        endif
        ! CALCULATING PART DEPENDENT ON POTENTIAL
        phitemp(ib)=phi
        if(langle.and.ldi) then
        if(lfrompdb(ib-2).and.lfrompdb(ib+1).and..not.lcoildih) then
            phi=phi-phi0(ib)
            if(ldisimp) then
                ene=0.5*CDH*phi*phi
                dvdp=-CDH*phi
            else
                ene=CDA*(1.d0-dcos(phi))+CDB*(1.d0-dcos(3.d0*phi))
                dvdp=CDA*sin(phi)+3.d0*CDB*sin(3.d0*phi)
            endif
        else
            j=min(inameseq(ib-1),3)
            k=min(inameseq(ib),3)
            sinfi=dsin(phi)
            cosfi=cosphi
            sin2fi=sinfi**2
            cos2fi=cosfi**2
            sincosfi=sinfi*cosfi
            ene=dihpot(1,j,k)+dihpot(2,j,k)*sinfi+dihpot(3,j,k)*cosfi
     +      +dihpot(4,j,k)*sin2fi+dihpot(5,j,k)*cos2fi
     +      +dihpot(6,j,k)*sincosfi
            dvdp=dihpot(2,j,k)*cosfi-dihpot(3,j,k)*sinfi
     +      +2.d0*(dihpot(4,j,k)-dihpot(5,j,k))*sincosfi
     +      +dihpot(6,j,k)*(cos2fi-sin2fi)
        endif
        if(lsldh) then
            dvdp=dvdp*cofdih
            ene=ene*cofdih
        endif
        epot=epot+ene    
        fx(i1)=fx(i1)-dvdp*fi(1)
        fy(i1)=fy(i1)-dvdp*fi(2)
        fz(i1)=fz(i1)-dvdp*fi(3)
        fx(i2)=fx(i2)-dvdp*fj(1)
        fy(i2)=fy(i2)-dvdp*fj(2)
        fz(i2)=fz(i2)-dvdp*fj(3)
        fx(i3)=fx(i3)-dvdp*fk(1)
        fy(i3)=fy(i3)-dvdp*fk(2)
        fz(i3)=fz(i3)-dvdp*fk(3)
        fx(i4)=fx(i4)-dvdp*fl(1)
        fy(i4)=fy(i4)-dvdp*fl(2)
        fz(i4)=fz(i4)-dvdp*fl(3)
        endif
        endif
200     continue
!$omp enddo nowait
!$omp end parallel 
        return
        end

C THIS SUBROUTINE RETURN THE BOND ANGLE AT THE SECOND SITE
c ux1=x0(ib)-x0(ib-1)
c ux2=x0(ib)-x0(ib+1)
        subroutine bondangle(theta)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/four/rx(4),ry(4),rz(4)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        ux1=rx(2)-rx(1)
        uy1=ry(2)-ry(1)
        uz1=rz(2)-rz(1)
        uu1=sqrt(ux1*ux1+uy1*uy1+uz1*uz1)

        ux2=rx(2)-rx(3)
        uy2=ry(2)-ry(3)
        uz2=rz(2)-rz(3)
        uu2=sqrt(ux2*ux2+uy2*uy2+uz2*uz2)

        u12=ux1*ux2+uy1*uy2+uz1*uz2
        u12=u12/(uu1*uu2)

        u12=min(u12,1.d0)
        u12=max(u12,-1.d0)
        theta=dacos(u12)

        return
        end

C THIS SUBROUTINE RETURN THE DIHEDRAL ANGLE AT THE THIRD SITE

        subroutine dihedral(phi)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/four/rx(4),ry(4),rz(4)
        common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc

        ux1=rx(2)-rx(1)
        uy1=ry(2)-ry(1)
        uz1=rz(2)-rz(1)

        ux2=rx(3)-rx(2)
        uy2=ry(3)-ry(2)
        uz2=rz(3)-rz(2)

        ux3=rx(4)-rx(3)
        uy3=ry(4)-ry(3)
        uz3=rz(4)-rz(3)

        vx1=uy1*uz2-uz1*uy2
        vy1=uz1*ux2-ux1*uz2
        vz1=ux1*uy2-uy1*ux2
        vv1=vx1*vx1+vy1*vy1+vz1*vz1

        vx2=uy2*uz3-uz2*uy3
        vy2=uz2*ux3-ux2*uz3
        vz2=ux2*uy3-uy2*ux3
        vv2=vx2*vx2+vy2*vy2+vz2*vz2

        v12=vx1*vx2+vy1*vy2+vz1*vz2
        v12=v12/sqrt(vv1*vv2)
        v12=min(v12,1.d0)
        v12=max(v12,-1.d0)
        phi=dacos(v12)

        di=vx1*ux3+vy1*uy3+vz1*uz3
        if(di.lt.0.d0) phi=-phi

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE sort2(n,arr,ibarr)
        implicit double precision(a-h,o-z)
      INTEGER n,M,NSTACK           !    from "Numerical Recipes"
c      REAL arr(n)
       dimension arr(n)
                integer ibarr(n)
      PARAMETER (M=7,NSTACK=500)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
c      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
       a=arr(j)
       ib=ibarr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
       ibarr(i+1)=ibarr(i)
11        continue
          i=0
2         arr(i+1)=a
      ibarr(i+1)=ib
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
      itempb=ibarr(k)
        arr(k)=arr(l+1)
      ibarr(k)=ibarr(l+1)
        arr(l+1)=temp
      ibarr(l+1)=itempb
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
      itempb=ibarr(l+1)
          arr(l+1)=arr(ir)
       ibarr(l+1)=ibarr(ir)
          arr(ir)=temp
       ibarr(ir)=itempb
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
       itempb=ibarr(l)
          arr(l)=arr(ir)
       ibarr(l)=ibarr(ir)
          arr(ir)=temp
       ibarr(ir)=itempb
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
      itempb=ibarr(l+1)
          arr(l+1)=arr(l)
       ibarr(l+1)=ibarr(l)
          arr(l)=temp
       ibarr(l)=itempb
        endif
        i=l+1
        j=ir
        a=arr(l)
      ib=ibarr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
      itempb=ibarr(i)
        arr(i)=arr(j)
       ibarr(i)=ibarr(j)
        arr(j)=temp
       ibarr(j)=itempb
        goto 3
5       arr(l)=arr(j)
       ibarr(l)=ibarr(j)
        arr(j)=a
       ibarr(j)=ib
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          write(*,*)'NSTACK too small in sort'
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

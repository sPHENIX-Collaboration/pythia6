! ======================================================================
! eRHIC implementation of PYTHIA event generator.
! This file contains routines to read in a steer file, initialise
! PYTHIA and generate events.
!
! The code is written in fixed format for compatibility with
! the included PYTHIA and Radgen files.
!
! There is no program defined - the routines are called
! from C++ code, which provides the 'main' function.
! We therefore wrap routines and variables in a module to ease
! passing data back and forth between Fortran and C++.
!
! A note on variable naming: any variables, routines or types that
! are visible to C++ code are written with underscores between words
! rather than using mixed case i.e.
!  names_like_this rather than NamesLikeThis etc.
! This is because Fortran is case-insensitive so all variable
! and routine names appear as all lower case in C++. Using underscores
! to separate words gives greater readability on the C++ side.
!
! \todo Port radiative photon output to C++ code.
! ======================================================================
      module pythia6
      
      private ! Make things private by default
      
      ! These variables are all public so they can be
      ! shared with C++ code.
      
      integer, public:: printascii ! Flag to print ASCII output

      integer, public:: genevent ! Number of generations to get the current event
      double precision, public:: trueX ! Bjorken-x
      double precision, public:: trueW2 ! Invariant mass of hadronic final state
      double precision, public:: trueNu ! Photon energy in proton rest frame
      double precision, public:: F2
      double precision, public:: F1
      double precision, public:: R
      double precision, public:: sigma_rad
      double precision, public:: SigRadCor
      double precision, public:: EBrems
      double precision, public:: radgamE ! Energy of radiated photon
      double precision, dimension(3), public:: radgamp ! Momentum of radiated photon

      ! These need to be accessed by the module
      ! subroutines, but aren't needed externally.

      ! Logical unit to which ASCII output is written,
      ! is ascii output is requested.
      integer, parameter:: asciiLun = 29

      integer:: NEV ! Number of events to generate
      integer:: IEV ! Counts the number of the current event
      double precision:: pbeam ! Proton beam momentum (GeV/c)
      double precision:: ebeam ! Lepton beam momentum (GeV/c)
      double precision:: massp ! Proton mass
      double precision:: masse ! Lepton mass
      double precision:: pgamma ! Gamma factor (E/m) of beam proton
      double precision:: pbeta ! Beta (p/E) of beam proton

      ! The following subroutines are defined:
      !    initialise()
      !    finish()
      !    reset()
      ! The following functions are defined:
      !    generate()

      public:: initialise, finish, generate
      private:: reset
      
      contains

      ! ================================================================
      ! Initialisation routine.
      ! Sets module variables to default values.
      ! Reads input from steer file and directs PYTHIA to process it.
      ! Opens output ASCII file and prints header if ASCII output
      ! has been requested.
      ! Returns an integer with the following meaning:
      ! 0: Success
      ! <0: An error was encountered
      ! >0: Success, user requested only initialisation steps in steer
      ! If a value >0 is returned then event generation should not
      ! be performed - this indicates e.g. that radgen table
      ! initialisation was requested.
      ! ================================================================
      function initialise()

         include 'pythia.inc'              ! All PYTHIA commons blocks
         include "mc_set.inc"
         include "py6strf.inc"
         include "mcRadCor.inc"
         include "radgen.inc"
         include "phiout.inc"

         ! Added by liang 1/6/12
         ! Switches for nuclear correction
         common/PYNUCL/INUMOD, CHANUM, ORDER
         save/PYNUCL/
         
         ! Function variable
         integer:: initialise
         
         double precision INUMOD,CHANUM
         ! -------------------------------------------------------------
         ! Arrays storing times used in setting seed.
         ! -------------------------------------------------------------
         integer, dimension(3):: today
         integer, dimension(3):: now
         integer:: ORDER
         integer:: NPRT, I, ltype 
         integer:: idum1, idum2, initseed
         double precision:: sqrts
         double precision:: pbeamE, ebeamE, epznucl
         logical:: UseLut, GenLut
         ! -------------------------------------------------------------
         ! ASCII output file
         ! -------------------------------------------------------------
         character(len=256):: outname
         character(len=100) :: PARAM ! Parameters passed to pygive()
         
         ! Set some default values before reading anything.
         initialise = -1
         genevent = 0
         pbeam = 100. 
         ebeam = 4. 
         masse = PYMASS(11)
         massp = PYMASS(2212)
         NEV = 0
         IEV = 0
         iModel=0
         ltype=11
         
         ! Read the steer file passed on the command line.
         
         ! Read output file name
         read(*, *) outname
         ! Read lepton beam type 
         read(*, *) ltype 
         ! Read parameters for PYINIT call (beam and target particle energy).
         read(*, *) pbeam, ebeam
         ! Read number of events to generate, and to print.
         read(*, *) NEV,NPRT
         ! Read min/max x of radgen lookup table
         read(*, *) mcSet_XMin, mcSet_XMax
         ! Read min/max y of generation range      
         read(*, *) mcSet_YMin, mcSet_YMax
         ! Read min/max Q2 of generation range      
         read(*, *) mcSet_Q2Min, mcSet_Q2Max
         ! Read information for cross section used in radgen
         read(*, *) genSet_FStruct, genSet_R
         ! Read parameters of radcorr:
         ! do radcorr (1), generate look-up table (2)
         read(*, *) qedrad
         ! Read parameters for PYTHIA-Model = which generation is done     
         read(*, *) iModel
         ! Read target type mass and charge
         read(*, *) mcSet_TarA, mcSet_TarZ
         ! Read nuclear pdf parameter mass number A, charge number Z
         read(*, *) INUMOD, CHANUM
         ! Read nuclear pdf correction order
         read(*, *) ORDER
         ! Read information for cross section used in radgen
100      read(*, '(A)', end=200) PARAM
            call PYGIVE(PARAM)
         goto 100
         
         ! -------------------------------------------------------------
         ! Initialize PYTHIA.      
         ! -------------------------------------------------------------
200      write(*, *) '*********************************************'
         write(*, *) 'NOW all parameters are read by PYTHIA'
         write(*, *) '*********************************************'

         ! Getting the date and time of the event generation
         call idate(today)   ! today(1)=day, (2)=month, (3)=year
         call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        
         ! Take date as the SEED for the random number generation
         ! Use constant seed for debug purposes
         ! \todo Need to change this back for release!!!
         initseed = 12345!today(1) + 10*today(2) + today(3) + now(1) + 5*now(3)
         write(6,*) 'SEED = ', initseed
         call rndmq(idum1,idum2,initseed,' ')
        
         sqrts = sqrt(4. * pbeam * ebeam)
         write(*, *) '*********************************************'
         write(*, *) 'proton beam energy:', pbeam, 'GeV'
         write(*, *) 'lepton beam energy:', ebeam, 'GeV'
         write(*, *) 'resulting sqrt(s):', sqrts, 'GeV'
         write(*, *) '*********************************************'
         
         ! proton is defined in positive z and as target
         P(2,1)=0.0  
         P(2,2)=0.0  
         P(2,3)=pbeam
         ! lepton is defined in negative z and as beam
         P(1,1)=0.0  
         P(1,2)=0.0  
         P(1,3)=-ebeam
         ! Find the masses of the beam particles
         if(mcSet_TarZ.eq.0) then
            massp=PYMASS(2112)
         else
            massp=PYMASS(2212)
         endif
         masse=PYMASS(ltype)
         ! Compute beam-dependent properties
         pbeamE=sqrt(pbeam**2+massp**2)
         pbeta=pbeam/pbeamE
         pgamma=pbeamE/massp
         ebeamE=sqrt(ebeam**2+masse**2)
         ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-ebeam)
         epznucl=-pgamma*pbeta*(ebeamE)+pgamma*(-ebeam)
         write(*, *) ebeamEnucl, ebeamE, epznucl, -ebeam
         mcSet_EneBeam=sngl(ebeamEnucl)
         
         ! Take actions dependent on radgen options
         if(iModel.eq.0) then
            UseLUT=.false.
            GenLUT=.false.
            qedrad=0
            MSTP(199)=0
            mcRadCor_EBrems=0.
         elseif(iModel.eq.1) then
            if(qedrad.eq.0) then
               mcRadCor_EBrems=0.
               UseLUT=.false.
               GenLUT=.false.
               MSTP(199)=1
            elseif(qedrad.eq.1) then
               mcRadCor_EBrems=0.
               UseLUT=.true.
               GenLUT=.false.
               MSTP(199)=1
               call radgen_init(UseLUT,GenLUT)
               write(*, *) 'I have initialized radgen'
            elseif(qedrad.eq.2) then
               write(*, *) 'radgen lookup table will be generated'
               mcRadCor_EBrems=0.
               UseLUT=.true.
               GenLUT=.true.
               MSTP(199)=1
               call radgen_init(UseLUT,GenLUT)
               write(*, *) 'lookup table is generated;'
               write(*, *) 'to run now pythia change parameter qedrad to 1'
               ! Set the return value to greater than zero to
               ! indiate that event generation should not proceed.
               initialise = 1
               return
            endif
         endif
         
         ! Initialise PYTHIA beam species and momenta
         if((mcSet_TarZ.eq.1).and.(ltype.eq.11)) then
            call pyinit ('3MOM','gamma/e-','p+',WIN)
         elseif((mcSet_TarZ.eq.1).and.(ltype.eq.-11)) then
            call pyinit ('3MOM','gamma/e+','p+',WIN)
         elseif((mcSet_TarZ.eq.0).and.(ltype.eq.-11)) then
            call pyinit ('3MOM','gamma/e+','n0',WIN)
         elseif((mcSet_TarZ.eq.0).and.(ltype.eq.11)) then
            call pyinit ('3MOM','gamma/e-','n0',WIN)
         endif

         ! If we ever want to simulate fixed target we need to change this
         ! win=ebeam
         ! call pyinit('fixt','gamma/e-','p+', WIN)

         ! -------------------------------------------------------------
         ! Open ascii output file
         ! -------------------------------------------------------------
         if(printascii .ne. 0) then
            open(asciiLun, file=outname)
            write(*, *) 'the outputfile will be named:', outname
            ! This is what we write in the ascii-file
            write(asciilun,*)' PYTHIA EVENT FILE '
            write(asciilun,*)'============================================'
            write(asciilun,*) 'I, ievent, genevent, subprocess, nucleon,'//
     +         ' targetparton, xtargparton, beamparton, xbeamparton,'//
     +         ' thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu, leptonphi,'//
     +         ' s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R, sigma_rad,'//
     +         ' SigRadCor, EBrems, photonflux, t-diff, nrTracks'
            write(asciilun,*)'============================================'

            write(asciilun,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)'//
     +         '  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)'
            write(asciilun,*)'============================================'
         endif ! printascii
         
         initialise = 0
         
      end function initialise
      ! ================================================================


      ! ================================================================
      ! Generates a single event.
      ! Returns 0 upon success.
      ! Returns 1 upon an error, or if the total number of events
      ! (counted by NEV) requested has been reached.
      ! ================================================================
      function generate()

         include 'pythia.inc'              ! All PYTHIA commons blocks
         include "mc_set.inc"
         include "py6strf.inc"
         include "mcRadCor.inc"
         include "radgen.inc"
         include "phiout.inc"

         integer:: generate ! Function return value
         integer:: nrtrack
         integer:: i ! loop variable
         ! Tracks the number of trials between successful event
         ! generations
         integer, save:: lastgenevent = 0
         
         ! Return if we have reached the maximum event count
         if(.not. IEV .lt. NEV) then
            generate = 1
            return
         endif
         
         ! Clear previous contents of publicly available variables
         call reset()
         
         ! call PYTHIA event generation
999      call PYEVNT
         if(MSTI(61).eq.1) then
             write(*, *) 'go back to PYEVNT call'
             goto 999
         endif

         genevent = NGEN(0, 3) - lastgenevent

         trueX =  VINT(307)/VINT(309)/(4*pbeam*ebeam)
         trueW2 = massp**2 + VINT(307)*(1/trueX-1)
         trueNu = (trueW2 + VINT(307) - massp**2)/(2.*massp)
         
         ! Calculate values for radiated photon
         if(mcRadCor_EBrems.gt.0.) then
            ! dplabg is from common/phiout/ in phiout.inc.
            EBrems=sqrt(dplabg(1)**2+dplabg(2)**2+dplabg(3)**2)
            ! Boost energy and momentum from lab frame.
            radgamE=pgamma*EBrems-pgamma*pbeta*dplabg(3)
            radgamp(3)=-pgamma*pbeta*EBrems+pgamma*dplabg(3)
            ! Copy dplabg(1) and (2) components to radgamp unchanged
            ! to avoid having to interface to the entire phiout
            ! common block in C++ code.
            do i = 1, 2
               radgamp(i) = dplabg(i)
            enddo
         else
            EBrems=0D0
            radgamE=0D0
            do i = 1, 3
               radgamp(i)=0D0 
            enddo
         endif

         ! N is from the pyjets common block
         ! Determine the number of tracks to print.
         ! If there was radiative emission there is one
         ! more track than from the pythia record.
         if(mcRadCor_EBrems.gt.0.) then
            nrtrack = N + 1
         else
            nrtrack = N
         endif

         if((msti(1).ge.91).and.(msti(1).le.94)) then
            msti(16)=0
         endif

         F2 = py6f2
         F1 = py6f1
         R = py6r
         sigma_rad = mcRadCor_Sigrad
         SigRadCor = mcRadCor_sigcor

         ! Print ASCII output if requested
         if(printascii .ne. 0) then
            ! Write the event header
            write(asciilun, 32) 0, IEV, genevent, msti(1), msti(12),
     +         msti(16), pari(34), msti(15), pari(33), pari(53),
     +         VINT(309), VINT(307), trueX, trueW2, trueNu,
     +         VINT(313), pari(14), pari(15), pari(16),
     +         pari(18),  pari(22), sngl(py6f2), sngl(py6f1),
     +         py6r, mcRadCor_Sigrad, mcRadCor_sigcor, EBrems,
     +         VINT(319), VINT(45), nrtrack 
32          format((I4,1x),(I10,1x),3(I4,1x),(I10,1x),f9.6,1x,
     +         I12,1x,
     +         2(f12.6,1x),7(f18.11,3x),12(f19.9,3x),I12,$/)
            write(asciilun, *)'============================================'
            ! Write track records
            ! N is from the pyjets common block
            do i = 1, N
               ! Check that the parent index K(I,3) of this particle is in
               ! a sensible range.
               if(K(I,3).le.nrtrack) then
                  write(asciilun,34) I,K(I,1),K(I,2),K(I,3),K(I,4),K(I,5),
     +               P(I,1),P(I,2),P(I,3),P(I,4),P(I,5),
     +               V(I,1),V(I,2),V(I,3)
               endif
            enddo
            ! Write radiative photon
            if(mcRadCor_EBrems.gt.0.) then
               write(asciilun,34) nrtrack, 55, 22, 1, 0, 0,
     +            sngl(dplabg(1)),sngl(dplabg(2)),sngl(-radgamp),
     +            sngl(radgamE), 0., 0., 0., 0.
            endif
34          format(2(I6,1x),I10,1x,3(I8,1x),8(f15.6,1x),$/)
            write(asciilun,*)'=============== Event finished ==============='
         endif ! printascii
         lastgenevent=NGEN(0,3)
         
         ! Increment event count and return
         IEV = IEV + 1
         generate = 0

      end function generate
      ! ================================================================
      
      
      ! ================================================================
      ! Perform end-of-run actions.
      ! Print output to screen.
      ! Close ASCII file if it was requested.
      ! ================================================================
      subroutine finish()
      
         include 'pythia.inc'
      
         ! Print cross sections.
         call PYSTAT(1)
         call PYSTAT(4)

         write(*, *)"The charm mass used is: ", PMAS(4,1)

         ! Print the Pythia cross section which is needed to get an absolut 
         ! normalisation the number is in microbarns
         write(*, *)'==================================================='
         write(*, *)'Pythia total cross section normalisation:',
     +      pari(1)*1000, ' microbarn'
         write(*, *)'Total Number of generated events', MSTI(5)
         write(*, *)'Total Number of trials', NGEN(0,3)
         write(*, *)'==================================================='
      
         if(printascii .ne. 0) then
            close(asciilun)
         endif

         ! Check pdf status       
         call PDFSTA

      end subroutine finish
      ! ================================================================
      
      
      ! ================================================================
      ! Reset publicly accessible values to default values.
      ! ================================================================
      subroutine reset()
         genevent = 0
         trueX = 0.D0
         trueW2 = 0.D0
         trueNu = 0.D0
         F2 = 0.D0
         F1 = 0.D0
         R = 0.D0
         sigma_rad = 0.D0
         SigRadCor = 0.D0
         EBrems = 0.D0
         radgamE = 0.D0
         radgamp = 0.D0
      end subroutine reset
      ! ================================================================
      
      end module pythia6
! ======================================================================

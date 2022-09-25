MODULE DOUBLE

IMPLICIT NONE 

integer,Parameter  :: DP=Kind(1d0) 

END MODULE DOUBLE

program aovmain

!! Purpose:<br>
!!`Aov` is a tool for period analysis of irregularly sampled time series, 
!!such as stellar and planetary transit light curves. It relies on state-of-art
!!statistically poweful, robust and fast algorithms implemented in
!!my `AovPer` package, partly based on my original research and
!!publications. The references are given in concerned routine descriptions.<br>

!!Usage:<br>
!!The example of calling the periodogram/frequency spectrum routines is 
!!provided in the form of `AovDrv` routine in
!!`AovPer` module. The observations file required by the main program
!!is an ASCII file named by `STAR`
!!parameter, containing 2 or 3 columns (time, values and optional errors/weights,
!!depending on `ERRIN`).<br>

!!Control of the code operation is by means of a configuration/parameter 
!!file `aov.par` and/or by command line arguments. Their default 
!!values are set in `AovPrms` module. For test purposes we
!!provide a routine simulating test data and sample output SIMUL1_reference.txt.<br>
!!The basic input parameter is the number of model parameters NPAR,
!!corresponding to NPAR phase bins or NPAR/2 Fourier harmonics.
!!I was using 3<=NPAR<=9 and I do not recommend exceeding these limits
!!in general period surveys. For transient search use much higher values 
!!of NPAR, so that transit width ~1/NPAR. In theory, the CPU time
!!is proportional to ~NPAR, but as NPAR increases the spectral
!!resolution also increases, so it is advisable to decrease
!!the frequency step by 1/NPAR (done automatically by the code),
!!which causes CPU time scaling like ~NPAR^2. There
!!are possible additional CPU savings, but these are not
!!worth worrying about at this stage. NPAR=5 which
!!yields oversampling of the frequency spectrum worked well
!!for near sinusoidal oscillations. FSAMPL needs to be set
!!=>3. The FSAMPL=4 choice means 4 frequencies covering a spectral line.
!!This is in my opinion close to optimal number as the code fits a 
!!parabola to find peak frequency up to a fraction of frequency increment.
!!But if the speed is not a problem one can use e.g. FSAMPL=7.<br>
!!The fastest are `AOV` or `AOVW`, depending on `ERRIN`, yet they
!!suffer from an artificial noise/residuals due to sharp steps between phase bins.
!!`AMH` and `AMHW` employ a smooth Fourier model hence are statistically better, yet
!!about 10 times slower. All four algorithms complexity scales proportional
!!to number of observations, frequencies and model parameters.<br>

!!Methods:<br>
!!The algorithms employ either bins and sub-bins or projection onto
!!ortogonal polynomials yielding as fast implementation as possible when
!!irregular sampling with large gaps make FFT impractical.
!!Use of orthogonal functions, either phase binning or 
!!trigonometric polynomials results in statistically near optimal
!!behaviour in that cross-talk/correlation between model elements is minimal.
!!Because of that the algorithms remain fairly robust against outliers. 
!!More importantly use of orthogonal functions permits 
!!application of the classical Fisher AnoVa statistics to
!!evaluate frequency spectrum and reliability of any detected periods. 
!!Depending on `SPTYPE` the spectrum obeys either
!!Fisher-Snedecor F or beta probability distributions.<br>

!!Experience:<br>
!!The original periodogram/frequency spectrum routines were used on
!!over 10^5 stars, as part of f95 TATRY package in OGLE and LCVSS
!!mass photometry surveys and in EROS c package 
!!All Star Searcher, among others used to analyse 318 light curves
!!of suspected variable stars in the field M31C. For the Cepheids
!!and the eclipsing binaries. Tests against periods set by manual 
!!inspection of periodograms and light curves were performed 
!!on 25 000 light curves from ASAS South. The comparision involved the most
!!significant periods only. The periods always 
!!agreed or were 2om or 3om or om/2 sub-/harmonics of the published frequency om.
!!For Cepheids NPAR=5 worked well. For eclipsing binaries with proximity 
!!effects (e.g. W UMA) or for sinusoidal variables NPAR=3 was better,
!!however 2om ambiguities occured.

!!Trials demonstrated that the aov detection criterion
!!used in the present code works better than the Laefler-Kinman,
!!PDM and Stetson criteria for eclipsing or pulsating variable stars.
!!However, picking right NPAR for your problem requires some
!!trial runs before decision on the final NPAR. For some applications,
!!such as transit search it sis recommended to perform several runs,
!!e.g. for NPAR=15, 30 and 60, to scan a range of transit widths.<br>

!! History:<br>
!! -The period search/frequency spectrum routines were developed since 
!! 1989 (AOV phase binning routines) and 1996 (AOVMH Fourier series routines). <br>
!! -Substantial revision/upgrade of the AOV routines, involving sub-bins 
!! and transit search was performed in 2003-2004 and included f95 and c 
!! variants of all routines.<br>
!! -Subsequent revision in 2018 removed small inconsistencies and combined 
!! AOV and AOVTR routines into one AOV performing both ordinary phase folding
!! and binning or transit search. Present f95 and c packages are compatible 
!! to a high degree.
!!
!!Copyrights:<br>
!!The whole package is subject to copyrights by its author, 
!!(C) A. Schwarzenberg-Czerny 2003-2018,  alex@camk.edu.pl 
!!Its distribution is free, except that distribution of modifications is
!!prohibited.
!!Please acknowlege use of routines by citation of papers listed in the
!!corresponding routine headers. <br>

!! Compilation:
!!AOV is compiled simply by:
!!
!!gfortran -o aov aovconst.f90 aovprms.f90 aovsub.f90 aovobs.f90 aovspec.f90 aovper.f90 aov.f90
!!
!! For debugging/testing use:
!!gfortran -g -Og -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2003 -pedantic -fbacktrace -o aov aovconst.f90 aovprms.f90 aovsub.f90 aovobs.f90 aovspec.f90 aovper.f90 aov.f90

!! Test run on simulated data: >./aov 'STAR="SIMUL"' 
!! Next to run on the same data from a file >./aov 'STAR="SIMUL1"'

  use DOUBLE
  use aovconst
  use mprms
  use aovper 
  use aovsub
  use aovobs
  use aovspec

  implicit none
! variables
  logical :: simul
  real(SP) :: dfm,thm,qual
  integer  :: ib,j,imet,metcnt
  real(DP) :: value,par1,par2
! To create output files in the desired directory
  character(1000)            :: outputFileName
  integer, parameter       :: outFileUnit = 100000
  CHARACTER(*), PARAMETER  :: fileplace = "/work/saikruba/zdcf/TEST_5/PDM/test1/RN4/aov_output/R2/"

  real(TIME),allocatable:: f(:),fb(:)
  type(SPEC_T) :: spec
  type(OBS_T)  :: obs

  io=prog_args('aov')
  if (io/=0) then
    print *,'Usage:'
    print *,'./aov  [[option=value]...] starname'
    print *,'where option=FR0|FRU|FRS|NH2|METHOD|SPTYPE|NCOV|DIR|EXT'
    print *,'METHOD=aov|atr|amh|psp|aovw|atrw|amhw|pspw'
    print *,'and value & starname must not contain ''='' sign'
    print *,'bracket with "" or '' any value containing SPACE'
    stop
  endif
  NCOV=MAX(1,NCOV)
  FRU=MAX(FR0,FRU)
  NPAR=MAX(2,NPAR)
  NBEST=MAX(1,NBEST)
  METHOD=toupper(METHOD)
  SPTYPE=toupper(SPTYPE)
  EPOCH=toupper(EPOCH)
  allocate(fb(NBEST))
  
  value = 0
print *,
print *,'  AOV/f90 Periodogram Routines, version '//VER
print *,'  by (C) Alex Schwarzenberg-Czerny alex@camk.edu.pl'
print *
Call print_args(6)
print *

! get data ...
  simul=(trim(STAR)=='SIMUL') 
  if(simul) then
! ...from simulation 
    obs=test_obs(100) ! also outputs simulated data to SIMUL1.dat)
  else
! ... or from a file. NOTE: appends file prefix/extension, if any 
    msg=trim(DIR)//trim(STAR)//trim(EXT)
    obs=get_obs(trim(msg)) 
  endif
  
  write(*,*)"The DIR:",msg

! process data/observations
  if ( minval(obs%t(2:)-obs%t(:obs%no-1)) < 0. ) &
      write(*,*) 'Warning: time not sorted in '//trim(STAR)
  write(*,*) 'Are errors given? ',obs%weights
  
  if (obs%no<=NPAR) then
    Print *, 'Too few observations, abandoning '//trim(STAR)
  else 
    metcnt=1
! Set up frequency grid and get periodogram
    spec = fgrid(obs%t,NPAR/2,FR0,FRU,FRS)
    allocate(f(spec%nfr),stat=io)
    if(io/=0) then
      Print *,"Failed memory allocation for f"
      stop
    endif
    forall(j=1:spec%nfr) f(j)=spec%FRS*(j-1)+spec%FR0
    if ( (obs%tmax - obs%tmin) * spec%frs > 0.3 .and.VERBOSE) &
      Print *, " AOV: warning: undersampling in frequency"
!    do imet=1,8
!      method=trim(meths(imet))
      if (metcnt==1) print '(a,2g15.8,i8)','grid: ',spec%FR0,spec%FRS,spec%nfr

! Calculate periodogramme 
      io=aovdrv(METHOD,sptype,obs%t,obs%v,obs%w,NPAR,NCOV,spec)
      par1 = (obs%no-1)
      par2 = (obs%no-NPAR)
      value = par1/par2
     
! Write result periodogram into file
     
      open(3,file=trim(fileplace)//trim(STAR)//trim(METHOD)//'.res',status='unknown',iostat=io)
      if (io/=0) write(*,*) 'Failed to open '//trim(STAR)//'.res'
      write(3,'(g15.8,f12.8)') (f(j),value*spec%th(j),j=1,spec%nfr)
      close(3)
! find & remove peaks
      do ib=1,NBEST     
        call peakrm(obs%tmax - obs%tmin,spec,fb(ib),thm,dfm)
        if (ib==1) qual=log10(thm)
      enddo
!  Print summary 
      if(metcnt==1) Print *,"meth star  no_obs quality best_frequencies"
      print '(1x,a4,1x,a6,i5,f8.3,10f10.5)', METHOD, &
         trim(STAR),obs%no,qual,(fb(ib),ib=1,NBEST)
      metcnt=metcnt+1
!    end do
    deallocate(f)
    call clr_spec(spec)
  endif
  call clr_obs(obs)

print *,'AOV finished'
print *

end program aovmain

 






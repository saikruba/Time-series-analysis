! gfortran zdcf_subQRM_t2.f90 -o zsqrMt2.out -Wall -L/usr/lib/pgplot -L`pwd` -lpgplot -lpng -lz -L/usr/X11R6/lib -lX11 ; zsqrMt2.out

MODULE DOUBLE

IMPLICIT NONE 

integer,Parameter  :: DP=Kind(1d0) 

END MODULE DOUBLE


program letsmakealc4

use double 
real    ::  x1,y1, fcut
real    ::  beta,gamma, DTsamp, tmult, pfunc(7), temp4, temp5, temp6
real    ::  logfbreak, fL, R, Q, temp(5000), sum3, value2, sum2, temp1, temp2,temp3, sum4,sum5, sum6
real    ::  time1,time2,time3,time4,time5,time6,time7,time8,fwhm
real    ::  av1,av2,av3,av4,av5,av6,av7,av8
real    ::  t1_avg(1000),t2_avg(1000),t3_avg(1000),t4_avg(1000)
real    ::  t5_avg(1000),t6_avg(1000),t7_avg(1000),t8_avg(1000)
real    ::  xline1(20),xline2(20),xline3(20),xline4(20),xline5(20),xline6(20),xline7(20),xline8(20),xline9(20)
real                              :: slope(20),checkR(20)
real, allocatable                 :: t(:),dtime(:),dflux(:)
real(DP)                          :: t1(20,1000),t2(20,1000),t3(20,1000),t4(20,1000)
real(DP)                          :: r2(20,1000),r4(20,1000),r6(20,1000),r8(20,1000)
real(DP)                          :: t5(20,1000),t6(20,1000),t7(20,1000),t8(20,1000)
real(DP)                          :: t_1(20,1000),t_2(20,1000),t_3(20,1000),t_4(20,1000)
real(DP)                          :: r_2(20,1000),r_4(20,1000),r_6(20,1000),r_8(20,1000)
real(DP)                          :: t_5(20,1000),t_6(20,1000),t_7(20,1000),t_8(20,1000)
real(DP)                          :: t1good(20,1000),t2good(20,1000),t3good(20,1000),t4good(20,1000)
real(DP)                          :: r2g(20,1000),r4g(20,1000),r6g(20,1000),r8g(20,1000)
real(DP)                          :: t5good(20,1000),t6good(20,1000),t7good(20,1000),t8good(20,1000)
real(DP)                          :: r2b(20,1000),r4b(20,1000),r6b(20,1000),r8b(20,1000)
real(DP)                          :: t1bad(20,1000),t2bad(20,1000),t3bad(20,1000),t4bad(20,1000)
REAL(DP)                          :: t5bad(20,1000),t6bad(20,1000),t7bad(20,1000),t8bad(20,1000)
real                              :: xx1(2),yy1(2), shift(50000), zerofreq(1000),fwhm_avg(1000),var(1000),stdev(1000)
real, dimension(:,:), allocatable :: fwhm2
real(8) ::  myshape, duration
integer ::  idumseed, typemod, Mfactor, N1, option, NSampd,BigN, nbin, IER, filename, filenum
integer ::  i,j, k,ma, a1, a2,a3,b3,l ,m,n,l2
character(80) outputFileName
real(DP) :: DTsim,sampt,corr(50000)
integer, parameter :: outFileUnit = 100000
CHARACTER(*), PARAMETER           :: fileplace = "/work/saikruba/zdcf/QR/M/acf_out/"
real, dimension(:,:),allocatable  :: acf,shiftd
integer                           :: no1(20),no2(20),no3(20),no4(20),no5(20),no6(20),no7(20),no8(20)
integer                           :: val1,val2,val3,val4,val5,val6,val7,val8
integer                           :: valu1,valu2,valu3,valu4,valu5,valu6,valu7,valu8
integer                           :: valr4,valr2
real, dimension(20) :: t4upline99,t4lowline99,t4upline98,t4lowline98,t4lowline97,t4upline97
real, dimension(20) :: t8upline99,t8lowline99,t8upline98,t8lowline98,t8lowline97,t8upline97
real, dimension(20) :: r4upline99,r4lowline99,r4upline98,r4lowline98,r4lowline97,r4upline97
real, dimension(20) :: r8upline99,r8lowline99,r8upline98,r8lowline98,r8lowline97,r8upline97
real, dimension(20) :: upline99,lowline99,upline98,lowline98,lowline95,upline95
integer :: upindex1,upindex2,upindex3,lowindex1,lowindex2,lowindex3
!==================================================================
real                  :: sumt1,sumt2,sumt3,sumt4,sumt5,sumt6,sumt7,sumt8
real,dimension(1000)  :: stdevt1, stdevt2,stdevt3,stdevt4,stdevt5,stdevt6,stdevt7,stdevt8
real,dimension(1000)  :: vart1,vart2,vart3,vart4,vart5,vart6,vart7,vart8
real                  :: svalt1,svalt2,svalt3,svalt4,svalt5,svalt6,svalt7,svalt8
!==================================================================
CHARACTER(len=10)                   :: folder(7)
real                                :: pratio
CHARACTER(len=80)                   :: fileplace1
integer                             :: nfolder
 

folder(1) = "01"
folder(2) = "1"
folder(3) = "10"
folder(4) = "100"
folder(5) = "103"
folder(6) = "104"
folder(7) = "105"



allocate(acf(10000,10000))
allocate(shiftd(10000,10000))
allocate(fwhm2(10000,10000))
 filename = 0

do nfolder = 1,7


 pratio = 10.**(nfolder-2)

 WRITE(fileplace1,"(A,A6)")"/work/saikruba/zdcf/QR/M/M",folder(nfolder)
 
 print*,"Fileplace:",fileplace1
 print*,"pratio:",pratio

upline99  = 0.0
upline98  = 0.0
upline95  = 0.0
lowline99 = 0.0
lowline95 = 0.0
lowline98 = 0.0

slope  = 0.0
checkR = 0.0

temp4 = 0.0
temp3 = 0.0
temp2 = 0.0
temp1 = 0.0

 av1 = 0.0
 av2 = 0.0
 av3 = 0.0
 av4 = 0.0


zerofreq  = 0.0

acf       = 0.0

fwhm      = 0.0

shiftd    = 0.0

logfbreak = 0.0

beta      = 0.0

gamma     = 0.0

duration  = 0.0

 c1       = 0
 
 c2       = 0 

 c3       = 0
 
 c4       = 0 

 c5       = 0

 filenum = 0

do loop2  = 1,14

temp5 = 0.0

temp6 = 0.0
!For PGLINE

xline1(loop2) = 4.0
!First minimum
xline2(loop2) = 7.5

xline3(loop2) = 11.0
!First maximum
xline4(loop2) = 14.5

xline5(loop2) = 18.1

xline6(loop2) = 21.8
!First minimum
xline7(loop2) = 25.4
!second maximum
xline8(loop2) = 29.0

xline9(loop2) = 0.0

zerofreq(loop2) = 0

     no1(loop2) = 0
     
     no2(loop2) = 0

     no3(loop2) = 0

     no4(loop2) = 0

     no5(loop2) = 0

     no6(loop2) = 0

     no7(loop2) = 0

     no8(loop2) = 0

val1  = 0
val2  = 0
val3  = 0
val4  = 0
val5  = 0
val6  = 0
val7  = 0
val8  = 0
valu1 = 0
valu2 = 0
valu3 = 0
valu4 = 0
valu5 = 0
valu6 = 0
valu7 = 0
valu8 = 0
valr2 = 0
valr4 = 0

 do ma=1,1000

  fwhm2(loop2,ma)=0.0
  !Set idumseed to some negative odd integer. Change it EVERY SINGLE TIME
  !you call the simlc routine
  corr = 0.0
 
  shift = 0.0

  filenum = filenum+1 
   
  filename = filename+1
   
  idumseed = -1-(filename*2)

!c---  typemod = 1 for unbroken-PL PSD; 
!c---  typemod = 2 for a sharply-broken PL PSD;
!c --- typemod = 3 for a very slow bend PSD
!c---  typemod = 4 for a bending power-law that's not as sharp as #3
!c---  typemod = 8 for a Lorentzian 

  typemod = 1

!c----  beta = high-frequency power-law slope: P(f) \propto f^-(beta)
!c----  gamma = low-frequency power-law slope: P(f) \propto f^-(gamma)
!c----  Keep them positive for red noise; we add in the negative later 
!c---- For unbroken PL PSD: logfbreak is ignored, but please pay attention
!c---- to the value in the denominator of the function in "myshape"!        
!c---- For a Lorentzian, free parms are R, Q, and fL

 N1 = 250

 sampt = 1.0 

 DTsim    = 86400.0 * 1.0
 
! DTsim    = 1.0
 
 Duration = DTsim * N1
      

 fL = 8.0e-7
 
 beta = 0.2 + (loop2*0.2) 

 slope(loop2) = beta
 
 gamma= beta      

 logfbreak = -6.0

 ampl = 1.5e4
 
 Q = 30 

! R = 0.75

  fcut = (10.**(logfbreak))


  R = sqrt((pratio*ampl)*((fL/fcut)**(-beta))*((fL*3.141)/(2*Q)))

  
  checkR(loop2) = R
       
!c---- Units of logfbreak are Hz, e.g., most AGN PSDs have breaks
!c---- 10^-6 to 10^-4 Hz, and many BHXRBs have breaks 10^-1 to 10 Hz.       
!c----  For unbroken PL PSD: Please set logfbreak to -6.0

 
!c----- For unbroken PL PSD (typemod=1): 
!c-----                      "ampl" denotes the value of the PSD, 
!c-----                      in units of /Hz, at f = 10^-6 Hz;

!c-----                      values of 10^3 - 10^4 /Hz are typical for AGN

!c---  For broken-PL PSD (typemod=2):  
!c---  "ampl" denotes amplitude at the break frequency in units of /Hz
!c---  Most broken-PL PSDs for both AGN and BHXRBs have
!c     (break freq)*(Ampl.@break) ~ 0.01 to 0.03
!c---  e.g., for a PSD with a break at f=10^-4 Hz, ampl ~ 100 /Hz 
!c---  e.g., for a PSD with a break at f=10^-6 Hz, ampl ~ 10000 /Hz 

!c---- For a very-slowly-bending PSD (typemod=3), or a 
!c     not-as-slowly-bending power-law (typemod=4), ampl does not 
!c     have the same meaning, and it will typically be in the
!c     neighborhood of 10^-3 - 10^-2.


!c--- Set desired simulation time and duration in units of SECONDS 
!c---  using N1 = Duration/DTsim


!c--- "M-factor" -- you want to create light curves that contain some
!c--- red noise leaked power. Do this by creating a light curve that
!c--- probes the PSD a factor Mfactor (at LEAST 5-10; 50-100 is better)
!c---- lower in frequency than (1/duration).
!c--- Simlc will create one light curve of length (Mfactor*duration)
!c---- You can then slice it up into Mfactor separate simulated light curves.
!c---- Mfactor is integer-valued.   
 Mfactor = 5.0

 tmult   = 1.0
      
!C---- Collect the PSD parameters and put them into one array called "pfunc"
!c         pfunc(1)=logfbreak
!c         pfunc(2)=beta
!c         pfunc(3)=gamma
!c         pfunc(7)=ampl

pfunc(5) = fL
 
pfunc(6) = Q

pfunc(7) = R

pfunc(1) = logfbreak
      
pfunc(2) = beta
      
pfunc(3) = gamma
      
pfunc(4) = ampl
      
!c-- create big lightcurve time column: Has N1*Mfactor entries

 BigN=Mfactor*N1

 allocate(t(BigN))
 allocate(dtime(BigN))
 allocate(dflux(BigN))

 t(1) = 0.0

 do i=2,BigN 

  t(i)= t(1) + (DTsim*(i-1))

 end do

!c--- Just in case we want to verify the big array:
!c       open(unit=4,file='tempmons',status='unknown');
!c       do i=1,BigN 

 !c       write(4,*) i,  t(i)          
!c        end do
!c       close(unit=4)


 30   call simlc(t, BigN,    &
            DTsim,pfunc,idumseed,   &
             typemod,tmult,dtime,dflux)   

!c--- dtime & dflux are the outputted light curve time column, 
!c---    flux column for the just-simulated lightcurve;
!c---    time resolution DTsim;  total entries = Mfactor*N1 
        

!c-- Output in units of DAYS (because it's typical for long-term AGN monit.)
!c--But edit "dtime(i)/86400.0d0" to "dtime(i)" if you want SECONDS.
        dtime = dtime/86400.0

!      call DCF(dtime,dflux,N1,sampt,ma,acf,nbin,shiftd,shift,corr,fwhm2,loop2)

       write (outputFileName, '(a, I0, a)') 'QR',filenum,'.dat'

       open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

          do i=1,N1
             write(outFileUnit,*) dtime(i),dflux(i)

          end do
     
       close(unit=outFileUnit)


       call zdcf(dtime,dflux,nbin,N1,ma,acf,shiftd,fwhm2,loop2)
    
      print*,"%-----------------------NBIN",nbin

      nbin = 83 

       do i = 1,nbin
 
          shift(i) = shiftd(ma,i)
          
          corr(i)  = acf(ma,i)


       end do


       write (outputFileName, '(a, I0, a)') 'acfQR',filenum,'.dat'

       open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

         do i=1,nbin-1
          write(outFileUnit,*) shiftd(ma,i),acf(ma,i)
        end do 
     
       close(unit=outFileUnit)


       call parameters(shift,corr,nbin,N1,ma,t1,t2,t3,t4,t5,t6,t7,t8,r2,r4,r6,r8,loop2)
                    

deallocate(t)
deallocate(dtime)
deallocate(dflux) 

!ccc--end the ma loop:    
end do

  write (outputFileName, '(a, I0, a)') 'zcheckvalQR',nfolder,'.dat'

  open(unit=1, file=fileplace//outputFileName,Access = 'append',status='unknown')

   do i = 1, ma-1  

   write(1,*)slope(loop2),t4(loop2,i),r4(loop2,i),t8(loop2,i),r8(loo2,i)

  end do

 close(unit=1)

do m = 1,nbin-1

 do l2 = 1,ma-2

  do l =l2+1,ma-1

   if(acf(l2,m) .gt. acf(l,m)) then
 
     temp5 = acf(l2,m)

!     temp6 = shiftd(l,m)

     acf(l2,m) = acf(l,m)

!     shiftd(l,m) = shiftd(l+1,m)

     acf(l,m) = temp5

!     shiftd(l+1,m) = shiftd(l,m)  

   end if 

  end do

 end do

 end do


upindex1  = 999

lowindex1 = 2

upindex2  = 998

lowindex2 = 3

upindex3  = 997

lowindex3 = 4


 write (outputFileName, '(a, I0, a)') 'acforderQR',loop2,'.dat'

      open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

        do i=1,nbin-1
         write(outFileUnit,*) shiftd(upindex2,i),acf(upindex1,i),acf(lowindex1,i),acf(upindex2,i),&
         acf(lowindex2,i),acf(upindex3,i),acf(lowindex3,i)
       end do 
     
 close(unit=outFileUnit)



   do i=1,ma-1
      
        if(t1(loop2,i) .gt. 0.0) then

         no1(loop2) = no1(loop2)+1

         t_1(loop2,no1(loop2)) = t1(loop2,i)
        
  
        end if

        if(t2(loop2,i) .gt. 0.0) then

         no2(loop2) = no2(loop2)+1

         t_2(loop2,no2(loop2)) = t2(loop2,i)

         r_2(loop2,no2(loop2)) = r2(loop2,i)

  
        end if

 
        if(t3(loop2,i) .gt. 0.0) then

         no3(loop2) = no3(loop2)+1

         t_3(loop2,no3(loop2)) = t3(loop2,i)

       
  
        end if


 
        if(t4(loop2,i) .gt. 0.0) then

         no4(loop2) = no4(loop2)+1

         t_4(loop2,no4(loop2)) = t4(loop2,i)

         r_4(loop2,no4(loop2)) = r4(loop2,i)
  
        end if

        if(t5(loop2,i) .gt. 0.0) then

         no5(loop2) = no5(loop2)+1

         t_5(loop2,no5(loop2)) = t5(loop2,i)

        
  
        end if

  
 
        if(t6(loop2,i) .gt. 0.0) then

         no6(loop2) = no6(loop2)+1

         t_6(loop2,no6(loop2)) = t6(loop2,i)

          r_6(loop2,no6(loop2)) = r6(loop2,i)
  
        end if

        if(t7(loop2,i) .gt. 0.0) then

         no7(loop2) = no7(loop2)+1

         t_7(loop2,no7(loop2)) = t7(loop2,i)

         
  
        end if

 
        if(t8(loop2,i) .gt. 0.0) then

         no8(loop2) = no8(loop2)+1

         t_8(loop2,no8(loop2)) = t8(loop2,i)

         r_8(loop2,no8(loop2)) = r8(loop2,i)

         
  
        end if


     end do

 write (outputFileName, '(a, I0, a)') 'zcheckvaluesMT4_',loop2,'.dat'

      open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

        do i=1,no4(loop2)
         write(outFileUnit,*)slope(loop2),t_4(loop2,i),r_4(loop2,i)
       end do 
     
 close(unit=outFileUnit)



 write (outputFileName, '(a, I0, a)') 'zcheckvaluesMT8_',loop2,'.dat'

      open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

        do i=1,no8(loop2)
         write(outFileUnit,*)slope(loop2),t_8(loop2,i),r_8(loop2,i)
       end do 
     
 close(unit=outFileUnit)

! ========================= checking the t -parameters ===========================

    do j = 1, no1(loop2)
       
      if(t_1(loop2,j) >= 3.0 .and. t_1(loop2,j) <= 5.0) then
        
        val1 = val1+1 

        t1good(loop2,val1) = t_1(loop2,j)
        
      else
        valu1 = valu1 +1
        t1bad(loop2,valu1) = t_1(loop2,j)
 
      end if
   

    end do

    do j = 1, no2(loop2)
       
      if(t_2(loop2,j) >= 6.0 .and. t_2(loop2,j) <= 9.0) then
        
        val2 = val2+1 
        t2good(loop2,val2) = t_2(loop2,j)
        r2g(loop2,val2) = r_2(loop2,j)

        if(r2g(loop2,val2) <= 0.5) then     
              valr2 = valr2+1
        end if


      else
        valu2 = valu2 +1
        t2bad(loop2,valu2) = t_2(loop2,j)
 
        

      end if
   

    end do

    do j = 1, no3(loop2)
       
      if(t_3(loop2,j) >= 10.0 .and. t_3(loop2,j) <= 13.0) then
        
        val3 = val3+1 

        t3good(loop2,val3) = t_3(loop2,j)

      else
        valu3 = valu3 +1
        t3bad(loop2,valu3) = t_3(loop2,j)
 
      end if
   

    end do

    do j = 1, no4(loop2)
       
      if(t_4(loop2,j) >= 12.0 .and. t_4(loop2,j) <= 18.0) then
        
        val4 = val4+1 

        t4good(loop2,val4) = t_4(loop2,j)
        r4g(loop2,val4) = r_4(loop2,j)

        if(r4g(loop2,val4) >= 0.5) then     
              valr4 = valr4+1
        end if

      else
        valu4 = valu4 +1
        t4bad(loop2,valu4) = t_4(loop2,j)
        r4b(loop2,valu4) = r_4(loop2,j)  

      end if
   

    end do


    do j = 1, no5(loop2)
       
      if(t_5(loop2,j) >= 16.0 .and. t_5(loop2,j) <= 21.0) then
        
        val5 = val5+1 

        t5good(loop2,val5) = t_5(loop2,j)

      else
        valu5 = valu5 +1
        t5bad(loop2,valu5) = t_5(loop2,j)
 
      end if
   

    end do

    do j = 1, no6(loop2)
       
      if(t_6(loop2,j) >= 19.0 .and. t_6(loop2,j) <= 25.0) then
        
        val6 = val6+1 

        t6good(loop2,val6) = t_6(loop2,j)
        r6g(loop2,val6) = r_6(loop2,j)     

      else
        valu6 = valu6 +1
        t6bad(loop2,valu6) = t_6(loop2,j)
        r6b(loop2,valu6) = r_6(loop2,j) 

      end if
   

    end do

    do j = 1, no7(loop2)
       
      if(t_7(loop2,j) >= 21.0 .and. t_7(loop2,j) <= 29.0) then
        
        val7 = val7+1 

        t7good(loop2,val7) = t_7(loop2,j)

     else
        valu7 = valu7 +1
        t7bad(loop2,valu7) = t_7(loop2,j)
 
      end if
   

    end do

    do j = 1, no8(loop2)
       
      if(t_8(loop2,j) >= 24.0 .and. t_8(loop2,j) <= 33.0) then
        
        val8 = val8+1 

        t8good(loop2,val8) = t_8(loop2,j)
        r8g(loop2,val8) = r_8(loop2,j)


      else
        valu8 = valu8 +1
        t8bad(loop2,valu8) = t_8(loop2,j)
        r8b(loop2,valu8) = r_8(loop2,j) 

 
      end if
   

    end do

! ================================================================================


 write (outputFileName, '(a, I0, a)') 'TR4goodQRM',loop2,'.dat'

      open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

        do k=1,val4
         write(outFileUnit,*) t4good(loop2,k),r4g(loop2,k)
        end do 
     
 close(unit=outFileUnit)


 write (outputFileName, '(a, I0, a)') 'TR8goodQRM',loop2,'.dat'

      open(outFileUnit, file=trim(fileplace1)//"/"//outputFileName,status="unknown")

        do k=1,val8
         write(outFileUnit,*) t8good(loop2,k),r8g(loop2,k)
       end do 
     
 close(unit=outFileUnit)


  write (outputFileName, '(a, I0, a)') 'zCFQRbadM',nfolder,'.dat'

  open(unit=2, file=fileplace//outputFileName,Access = 'append',status='unknown')

    write(2,*)slope(loop2),valu1,valu2,valu3,valu4,valu5,valu6,valu7,valu8
  
 close(unit=2)



  write (outputFileName, '(a, I0, a)') 'zCFQRgoodM',nfolder,'.dat'

  open(unit=3, file=fileplace//outputFileName,Access = 'append',status='unknown')
  
   write(3,*)slope(loop2),val1,val2,val3,val4,val5,val6,val7,val8 
  
 close(unit=3)


  write (outputFileName, '(a, I0, a)') 'zCFQRgoodMR',nfolder,'.dat'

  open(unit=3, file=fileplace//outputFileName,Access = 'append',status='unknown')
  
   write(3,*)slope(loop2),valr2,valr4
  
 close(unit=3)

! ================================================================================


temp1 = 0.0

temp2 = 0.0

do k = 1,no4(loop2)-1
   
  do l = k+1 , no4(loop2)
                 
   if (t_4(loop2,k) > t_4(loop2,l)) then
      
     temp1 = t_4(loop2,k)
      
     t_4(loop2,k) = t_4(loop2,l)
      
     t_4(loop2,l) = temp1

    end if
    
   end do
   
 end do


do k = 1,no4(loop2)-1
   
  do l = k+1 , no4(loop2)
                 
   if (r_4(loop2,k) > r_4(loop2,l)) then

     temp2 = r_4(loop2,k)
      
     r_4(loop2,k) = r_4(loop2,l)
      
     r_4(loop2,l) = temp2
      
    end if
    
   end do
   
 end do



temp1 = 0.0

temp2 = 0.0

do k = 1,no8(loop2)-1
   
  do l = k+1 , no8(loop2)
                 
   if (t_8(loop2,k) > t_8(loop2,l)) then
      
     temp1 = t_8(loop2,k)
      
     t_8(loop2,k) = t_8(loop2,l)
      
     t_8(loop2,l) = temp1
      
    end if
    
   end do
   
 end do

do k = 1,no8(loop2)-1
   
  do l = k+1 , no8(loop2)
                 
   if (r_8(loop2,k) > r_8(loop2,l)) then

     temp2 = r_8(loop2,k)
      
     r_8(loop2,k) = r_8(loop2,l)
      
     r_8(loop2,l) = temp2
      
    end if
    
   end do
   
 end do

upindex1  = (no4(loop2)-1)

upindex2  =(no4(loop2)-2)

upindex3  = (no4(loop2)-3)

 t4upline99(loop2) = (t_4(loop2,upindex1))

 t4lowline99(loop2)= (t_4(loop2,lowindex1))

 t4upline98(loop2) = (t_4(loop2,upindex2))

 t4lowline98(loop2) = (t_4(loop2,lowindex2))

 t4upline97(loop2) = (t_4(loop2,upindex3))

 t4lowline97(loop2) = (t_4(loop2,lowindex3))

 r4upline99(loop2) = (r_4(loop2,upindex1))

 r4lowline99(loop2)= (r_4(loop2,lowindex1))

 r4upline98(loop2) = (r_4(loop2,upindex2))

 r4lowline98(loop2) = (r_4(loop2,lowindex2))

 r4upline97(loop2) = (r_4(loop2,upindex3))

 r4lowline97(loop2) = (r_4(loop2,lowindex3))


upindex1  = (no8(loop2)-1)

upindex2  = (no8(loop2)-2)

upindex3  = (no8(loop2)-3)


 t8upline99(loop2) = (t_8(loop2,upindex1))

 
 t8lowline99(loop2)= (t_8(loop2,lowindex1))


 t8upline98(loop2) = (t_8(loop2,upindex2))


 t8lowline98(loop2) = (t_8(loop2,lowindex2))


 t8upline97(loop2) = (t_8(loop2,upindex3))


 t8lowline97(loop2) = (t_8(loop2,lowindex3))


 r8upline99(loop2) = (r_8(loop2,upindex1))

 
 r8lowline99(loop2)= (r_8(loop2,lowindex1))


 r8upline98(loop2) = (r_8(loop2,upindex2))


 r8lowline98(loop2) = (r_8(loop2,lowindex2))


 r8upline97(loop2) = (r_8(loop2,upindex3))


 r8lowline97(loop2) = (r_8(loop2,lowindex3))


  write (outputFileName, '(a, I0, a)') 'zCFQRMT4_conf',nfolder,'.dat'

  open(unit=4, file=fileplace//outputFileName,Access = 'append',status='unknown')
  
   write(4,*)slope(loop2),t4upline99(loop2),t4lowline99(loop2),t4upline97(loop2),t4lowline97(loop2)
  
 close(unit=4)


  write (outputFileName, '(a, I0, a)') 'zCFQRMR4_conf',nfolder,'.dat'

  open(unit=5, file=fileplace//outputFileName,Access = 'append',status='unknown')
  
   write(5,*)slope(loop2),r4upline99(loop2),r4lowline99(loop2),r4upline97(loop2),r4lowline97(loop2)
  
 close(unit=5)

  write (outputFileName, '(a, I0, a)') 'zCFQRMT8_conf',nfolder,'.dat'

  open(unit=6, file=fileplace//outputFileName,Access = 'append',status='unknown')
  
   write(6,*)slope(loop2),t8upline99(loop2),t8lowline99(loop2),t8upline97(loop2),t8lowline97(loop2)
  
 close(unit=6)

  write (outputFileName, '(a, I0, a)') 'zCFQRMR8_conf',nfolder,'.dat'

  open(unit=7, file=fileplace//outputFileName,Access = 'append',status='unknown')
  
   write(7,*)slope(loop2),r8upline99(loop2),r8lowline99(loop2),r8upline97(loop2),r8lowline97(loop2)
  
 close(unit=7)


! loop2 ends
end do


  write (outputFileName, '(a, I0, a)') 'zCFQRM',nfolder,'.dat'

  open(unit=8, file=fileplace//outputFileName,status='unknown')
  
    do j = 1, loop2-1

    write(8,*)slope(j),no1(j),no2(j),no3(j),no4(j),no5(j),no6(j),no7(j),no8(j) 
   
    end do
 
  close(unit=8)

!nfolder ends
end do

deallocate(fwhm2)
deallocate(acf)
deallocate(shiftd)

contains


!////////////////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////////////////

 subroutine zdcf(ta,a,nbins,N1,ma,acf,shiftd,fwhm2,loop2)
  
    integer ::  mpnts,mbins,minbin
!!c****************************************************************************c
!c                                                                           c
!c   SETTING THE MAXIMAL NO. OF POINTS IN THE LIGHT CURVES (mpnts)           c
!c   AND THE MAXIMAL NO. OF BINS IN THE ACF / CCF (mbins).                   c
!c   ATTENTION: when not using alloc, mpnts should be also defined with      c
!c              the same value in the function alcbin.                       c
!c                                                                           c
      parameter (mpnts = 1550, mbins = 500)
!c****************************************************************************c
      parameter (minbin=mpnts)
      integer :: i,j,n
      integer :: na,nb,nbins,used,unused
      integer :: minpts
      logical :: autocf,no0lag,nwdata,unsmpl
!      real,allocatable :: ta(:), tb(:) , a(:) , b(:)
      integer , intent(in) :: N1,ma,loop2
      real, dimension(:,:), allocatable, intent(inout):: fwhm2
      real :: corr1,fwhm1,corr3,fwhm3,corr2(420000)
      real, dimension(:,:),allocatable,intent(inout) :: acf,shiftd
      real    :: ta(mpnts),tb(mpnts)
      real    :: a(mpnts),b(mpnts)
      real    :: erra(mpnts),errb(mpnts)
      real    :: tdcf(mbins),sigtm(mbins),sigtp(mbins),dcf(mbins)
      real    :: sigdcm(mbins),sigdcp(mbins)
      logical :: wuseda(mpnts),wusedb(mpnts)
      integer :: inbin(mbins),enough
      parameter (enough=11)
!c Areas for Monte Carlo estimates of the measurement errors
      integer :: ncarlo
      real    :: acarlo(mpnts),bcarlo(mpnts)
      real    :: avz(mbins)
!c Work areas for clcdcf
      real    :: wa(minbin),wb(minbin)
      real    :: wta(minbin),wtb(minbin)
      real    :: werra(minbin),werrb(minbin)
      integer :: wai(minbin,mbins),wbj(minbin,mbins)
      real    :: vtau(minbin,mbins)
!c   
      character(72) :: infil1,infil2,outfil,prefix
      integer       :: oerr,cerr,lprfix
      integer       :: seed
      real          :: expz,sigz,eps
      parameter (eps=1e-7)
!c   
      integer :: inpi
      character(1) :: inpc
      real    :: fishe,fishs
      real    :: rand
      real    :: z,r
      data seed/123/     
!c   
!c Reading program parameters

!c   
      print *,'ZDCF V1.2 begins.'
!      print *,'Auto-correlation or cross-correlation? (1/2):'
!      read *,inpi 
!      inpi = 1 for ACF 
      inpi = 1
      autocf = (inpi .eq. 1)

!c===========================================================
!     we dont need the output files ( directly plot the acf)   
!      print *,'Enter output files prefix:'
!      read '(a72)',prefix
!      do i = 1,len(prefix)
!         if (prefix(i:i) .eq. ' ') then
!            lprfix = i-1
!            goto 200
!         endif
!      enddo
!200   continue
!c===========================================================
   
      print *,'Uniform sampling of light curve? (y/n):'
!      read '(a1)',inpc
      inpc = 'y'
      unsmpl = (inpc .eq. 'y' .or. inpc .eq. 'Y')

!      print *,'Enter minimal number of points per bin (0 for default):'
!      read *,minpts
!c     if minpts = 0, use default
!      default value set for minpts

      minpts = 0 
 
      if (minpts .le. 0) minpts = enough
!c  
!     zero lag points ------------- required ??  
      print *,'Omit zero-lag points? (y/n):'
!      read '(a1)',inpc
      inpc = 'n'
      no0lag = (inpc .eq. 'y' .or. inpc .eq. 'Y') 
!c   


!      setting ncarlo to 100
!      print *,'How many Monte Carlo runs for error estimation?'
!      read *,ncarlo
       ncarlo = 100
      
      if (ncarlo.le.1) ncarlo = 0

!=========================================================================   
!c Reading the observed data  

tb = ta
b  = a
na = N1
nb = N1

erra = 0.0 
errb = 0.0

! OPEN(unit=1,file='sine_20s.dat', status='unknown')
 
  DO j=1,na

!   read(1,*)ta(j),a(j),erra(j)

   tb(j) = ta(j)

   b(j) = a(j)

  END DO

! CLOSE(1)

errb = 0.0
erra = 0.0 
 
!      print *,'Enter name of 1st light curve file:'
!      read '(a72)',infil1
!      if (.not. autocf) then
!         print *,'Enter name of 2nd light curve file:'
!         read '(a72)',infil2
!      endif


!   call redobs(na,nb,ta,tb,a,b,erra,errb,mpnts, &
!                  autocf,infil1,infil2)
!==========================================================================

    
!c   
!c Writing the condensed light curve data  --- removed
!c   
      print *,'ZDCF PARAMETERS:'
      print *,'Autocorrelation?  ',autocf
      print *,'Uniform sampling? ',unsmpl
      print *,'Omit zero lags?   ',no0lag
      print *,'Minimal # in bin: ',minpts
      print *,'# of MonteCarlo:  ',ncarlo 
      print *,'MonteCarlo seed:  ',seed 

!================================================================================
!c Estimating the effects of the measurement errors by Monte Carlo simulations
!c NOTE: need to insert values for na and nb ( since readobs subroutine is removed) 

  
      i = rand(seed)
      if (ncarlo .gt. 1) then
         do i = 1,mbins
            avz(i) = 0.0
         enddo
         nwdata = .true.
         do i = 1,ncarlo
            call simerr(a,acarlo,erra,na)
            if (autocf) then
               call clcdcf(ta,acarlo,erra,na, &
!c                        i  i      i    i
                          tb,acarlo,errb,nb, &
!c                        i  i      i    i
                  autocf,no0lag,unsmpl,used,unused,nwdata, &
!c                i      i      i      o    o      i/o    
                  minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp, &
!c                i/o    o    o     o     o   o      o    
                  inbin,minbin,nbins,mbins, &
!c                i    i      i/o   i
                  wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb, &
!c                w   w   w  w  w     w     w   w   w      w
                  vtau)
!c                o
            else
               call simerr(b,bcarlo,errb,nb)
!c   
!c calculating the discrete correlation function for Monte Carlo error approx.
!c   
               call clcdcf(ta,acarlo,erra,na, &
!c                        i  i i    i
                          tb,bcarlo,errb,nb, &
!c                        i  i i    i
                  autocf,no0lag,unsmpl,used,unused,nwdata, &
!c                i      i      i      o    o      i/o    
                  minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp, &
!c                i/o    o    o     o     o   o      o      
                  inbin,minbin,nbins,mbins, &
!c                i    i      i/o   i     
                  wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb, &
!c                w   w   w  w  w     w     w   w   w      w
                  vtau)
!c                o
            endif
            if (i .eq. 1) then
               print *
               print *,nbins,' bins actually used, ',unused, &
                   ' inter-dependent pairs discarded.'
            endif
!c   
!c The summing and averaging is done in z-space.
!c   
            do n = 1,nbins
               r = dcf(n)
               if (r .lt. -1.0+eps) then
                  r = -1.0+eps
               else if (r .gt.  1.0-eps) then
                  r =  1.0-eps
               endif   
               avz(n) = avz(n)+log((1.+r)/(1.-r))/2.
            enddo
         enddo
         do i = 1,nbins
            z = avz(i)/(ncarlo+0.)
            n = inbin(i)
            dcf(i) = tanh(z) 
            r = dcf(i)
            if (r .lt. -1.0+eps) then
               r = -1.0+eps
            else if (r .gt. 1.0-eps) then
               r =  1.0-eps
            endif
            z = log((1.+r)/(1.-r))/2.0
            sigz = fishs(r,n)
            expz = fishe(r,n)
            sigdcm(i) = r-tanh(expz-sigz)
            sigdcp(i) = tanh(expz+sigz)-r
         enddo
      else
!c   
!c calculating the discrete correlation function w/o Monte Carlo error approx.
!c   
         nwdata = .true.
         call clcdcf(ta,a,erra,na, &
!c                  i  i i    i
                    tb,b,errb,nb, &
!c                  i  i i    i
                    autocf,no0lag,unsmpl,used,unused,nwdata, &
!c                  i      i      i      o    o      i/o   
                    minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp, &
!c                  i/o    o    o     o     o   o      o      
                    inbin,minbin,nbins,mbins, &
!c                  o     i      i/o   i,   
               wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb, &
!c             w   w   w  w  w     w     w   w   w      w
              vtau)
!c             o
         print *
         print *,nbins,' bins actually used, ',unused, &
                ' inter-dependent pairs discarded.'
      endif
!c 
!===============================================================================  
!c printing the results: tdcf, dcf, sigdcf
!c   
!      outfil = 'out3.dcf'
!      open (unit=13,file=outfil,status='UNKNOWN',iostat=oerr,err=901)
!      print *
!      print '(1p,6(1x,a10),'' ('',a4,'')'')', &
!           ' tau      ','-sig(tau) ','+sig(tau) ',' dcf      ', &
!           ' -err(dcf)',' +err(dcf)','#bin'
!      print *
!      do i = 1,nbins
!         print '(1p,6(1x,e10.3),'' ('',i4,'')'')', &
!              tdcf(i),sigtm(i),sigtp(i), &
!              dcf(i),sigdcm(i),sigdcp(i),inbin(i) 
!         write (13,'(1p,6(1x,e10.3),1x,i4)')  &
!               tdcf(i),sigtm(i),sigtp(i), &
!              dcf(i),sigdcm(i),sigdcp(i),inbin(i)
!      enddo
!      close (unit=13,status='KEEP',iostat=cerr,err=902)
!      print *
!      print *,outfil(:lprfix+4),' written...'
!c   
!      print *,'ZDCF ended.'
!      stop
! 901  print *,'ERROR ON OPENING FILE - ERRCOD=',oerr
!      stop
! 902  print *,'ERROR ON CLOSING FILE - ERRCOD=',cerr
!      stop

       do i = 1,nbins
!          write(outFileUnit,*)shift(i),corr(i)
          shiftd(ma,i) = tdcf(i)

          acf(ma,i) = dcf(i)
        end do
         

         corr1 = 1.0
         fwhm1 = 0.0
         corr3 = 0.0
         fwhm3 = 0.0
         corr2(ma) = 0.5
      

!         print*,"fwhm2:",fwhm2(loop2,ma)

         do i = 2,nbins
!            print*,"corr:",corr(i)
            if(dcf(i) .le. corr2(ma)) then
               corr3 = dcf(i)
               fwhm3 = tdcf(i)
               corr1 = dcf(i-1)
               fwhm1 = tdcf(i-1)
               exit
            end if
         end do

        fwhm2(loop2,ma) = (((fwhm3-fwhm1)/(corr3-corr1))*(corr2(ma)-corr1))+fwhm1
      


! end of program
      end subroutine zdcf        

! removed the subroutine for reading observation
!c   
!c Generating a "true" signal by substracting gaussian noise from observations.

!c The error erra is the ABSOLUTE error in a.
!c   
      subroutine simerr(a,acarlo,erra,na)
      implicit none
      integer na
      real a(na),acarlo(na),erra(na)
      integer i
      real rndnrm
!c   
      do i = 1,na
         acarlo(i) = a(i)-erra(i)*rndnrm()
      enddo      
      return
      end
!c/////////////////////////////////////////////////////////////////////////////
!c   
!c Allocating points to bins
!c   
      subroutine alcbin(ta,na,tb,nb,wuseda,wusedb, &
                       no0lag,autocf,minpts, &
                       tdcf,sigtm,sigtp, &
                       wai,wbj,inbin,nbins,mbins,minbin,unsmpl, &
                       vtau)
      implicit none
!c****************************************************************************c
      real rslutn
      parameter (rslutn=0.001)
!c****************************************************************************c
      integer na,nb,nbins,mbins,minbin,minpts
      real ta(na),tb(nb)
      integer inbin(mbins)
      real tdcf(mbins),sigtm(mbins),sigtp(mbins)
      integer wai(minbin,mbins),wbj(minbin,mbins)
      real vtau(minbin,mbins)
      logical no0lag,autocf
      logical wuseda(na),wusedb(nb)
!!c-malloc(
!c    integer malloc
!c    pointer (iwaidx,waidx)
!c    pointer (iwbidx,wbidx)
!c    pointer (iwtau,wtau)
!c    pointer (iidx,idx)
!c    integer idx(1),waidx(1),wbidx(1)
!c    real wtau(1)
!c-malloc)
!c-no-malloc(
      integer mpnts
!c   ATTENTION: mpnts here should equal the value defined in the MAIN!!!
      parameter(mpnts=1550)
      integer idx(mpnts*mpnts)
      integer waidx(mpnts*mpnts),wbidx(mpnts*mpnts)
      real wtau(mpnts*mpnts) 
!c-no-malloc)
      integer i,j,p,np,pfr,pmax,incr,inb,nnegtv,plo,phi
      real tij,tolrnc
      logical frstym,unsmpl
      integer min0,max0
      data frstym/.true./
!c   
!c Allocating the dynamic work areas
!c   
!c-malloc(
!c    iwaidx = malloc(na*nb*4)
!c    iwbidx = malloc(na*nb*4)
!c    iwtau  = malloc(na*nb*4)
!c    iidx   = malloc(na*nb*4)
!c-malloc)
!c   
!c Calculate all the time lag points
!c   
      if (frstym) then
         frstym = .false.
         print  &
           '(''Binning with minimum of '',i3,'// &
           ''' points per bin and resolution of '',1p,g8.2,'' .'')', &
           minpts,rslutn
      endif
      np = 0
      if (autocf) then
         do i = 1,na
            do j = i,nb
               tij = tb(j)-ta(i)
               if (no0lag .and. tij .eq. 0.0) goto 11
               np = np+1
               wtau(np) = tij
               waidx(np) = i
               wbidx(np) = j
11             continue
            enddo
         enddo
      else
         do i = 1,na
            do j = 1,nb
               tij = tb(j)-ta(i)
               if (no0lag .and. tij .eq. 0.0) goto 12
               np = np+1
               wtau(np) = tij
               waidx(np) = i
               wbidx(np) = j
12             continue
            enddo
 2       enddo
      endif
!c   
!c Sort according to increasing time lag
!c   
      call srtind(wtau,np,idx,1)
!c   
!c calculating the tolerance level for lags to be considered the same
!c   
      tij = wtau(idx(np))-wtau(idx(1))
      tolrnc = tij*rslutn
!c   
!c Looping on bins
!c   
      nbins = 0
!c   
!c If binned CCF: binning from median time-lag upwards and backwards!
!c   
      if (autocf .or. unsmpl) then
         pfr = 1
         pmax = np
         incr = 1
      else
         pfr = np/2
         pmax = 1
         incr = -1
      endif
20    continue
      inb = 0
      tij = incr*1e30
      nbins = nbins+1
      if (nbins .gt. mbins) then
         print *,'Not enough work area - increase mbins in MAIN!'
         stop
      endif
      tdcf(nbins) = 0.0
!c   
!c Initialize used point flag vectors
!c   
      do i = 1,na
         wuseda(i) = .false.
      enddo
      do i = 1,nb
         wusedb(i) = .false.
      enddo
!c   
!c Collect points into bins that contain at least "enough" points,
!c but do not break points with the same lag (up to the tolerance)
!c into separate bins
!c   
      do i = pfr,pmax,incr
         p = idx(i)
!c Check whether bin is full
         if ( ( incr*(wtau(p)-tij-incr*tolrnc) .gt. 0.0 .and.  &
               inb .ge. minpts ) .or.  &
            (unsmpl .and. incr*(wtau(p)-tij-incr*tolrnc) .gt. 0.0) .or. & 
            (i .eq. pmax) ) then
!c   
!c Bin is full: Calculating tau and its std (before proceeding to the next 
!c bin, at label 20)  
!c   
            inbin(nbins) = inb
            tdcf(nbins) = tdcf(nbins)/inb
            if (unsmpl) then
               sigtm(nbins) = 0.0
               sigtp(nbins) = 0.0
               pfr = i
!c If not enough points in bin, ignore it
               if (inb .lt. minpts) nbins = nbins-1
               if (pfr .ne. pmax) goto 20
            else
!c If the last point is alone in its bin, we ignore it (to avoid subsequent
!c divisions by zero) and get immediately out of loop (label 40)
               if (inb .le. 1) then
                  nbins = nbins-1
                  goto 40
               endif
!c   
!c Finding the 0.3413 (+/-1sig) points above and below the mean
!c   
               if (incr .eq. 1) then
                  plo = pfr
                  phi = i-1
               else
                  plo = i+1
                  phi = pfr
               endif
               do p = plo,phi
                  if (wtau(idx(p)) .ge. tdcf(nbins)) then 
                     j = p
                     goto 50
                 endif
               enddo
50             continue
               p = max0(j-nint(float(j-plo)*0.3413*2.0),plo)
               sigtm(nbins) = tdcf(nbins)-wtau(idx(p))
               p = min0(j+nint(float(phi-j)*0.3413*2.0),phi)
               sigtp(nbins) = wtau(idx(p))-tdcf(nbins)
               pfr = i
               if (pfr .ne. pmax) goto 20
            endif
!c If no more points - get out of loop (label 40)
            goto 40
         endif
!c   
!c Adding another point to the bin...
!c   
         if (wuseda(waidx(p)) .or. wusedb(wbidx(p))) goto 30
         inb = inb+1
         if (inb.gt. minbin) then
            print *,'Not enough work area - increase minbin in MAIN!'
            stop
         endif
         wuseda(waidx(p)) = .true.
         wusedb(wbidx(p)) = .true.
         tij = wtau(p)
         tdcf(nbins) = tdcf(nbins)+tij
         vtau(inb,nbins) = tij
         wai(inb,nbins) = waidx(p)
         wbj(inb,nbins) = wbidx(p)
30       continue
      enddo
!c   
!c Binning is finished
!c   
40    continue
      if (.not. (autocf .or. unsmpl .or. incr .eq. 1)) then
         pfr = np/2+1
         pmax = np
         incr = 1
         nnegtv = nbins
         goto 20
      endif
!c   
!c If CCF (and NOT uniform sampling): Sort the bins into increasing
!c chronological order: The nnegtv negative bins are at the beginning but at 
!c reverse order.
!c   
      if (.not. (autocf .or. unsmpl)) then
         do i = 1,nnegtv/2
            j = nnegtv+1-i
            inb = inbin(i) 
            inbin(i) = inbin(j)
            inbin(j) = inb
            tij = tdcf(i)
            tdcf(i) = tdcf(j)
            tdcf(j) = tij
            tij = sigtp(i)
            sigtp(i) = sigtp(j)
            sigtp(j) = tij
            tij = sigtm(i)
            sigtm(i) = sigtm(j)
            sigtm(j) = tij
            do p = 1,max0(inbin(i),inbin(j))
               tij = vtau(p,i)
               vtau(p,i) = vtau(p,j)
               vtau(p,j) = tij
               inb = wai(p,i)
               wai(p,i) = wai(p,j)
               wai(p,j) = inb
               inb = wbj(p,i)
               wbj(p,i) = wbj(p,j)
               wbj(p,j) = inb
            enddo
         enddo
      endif
!c   
!c Freeing the allocated work areas
!c   
!c-malloc(
!c    call free(iwaidx)
!c    call free(iwbidx)
!c    call free(iwtau)
!c    call free(iidx)
!c-malloc)
!c   
!c Emergency exit - something is very wrong...
!c   
      if (nbins .eq. 0) then
         print *,'alcbin: No bin contains enough points - stopping!'
         stop
      endif
      return
      end
!c////////////////////////////////////////////////////////////////////////////
!c Calculating the discrete correlation function.
!c POSITIVE lag values mean b lags after a.
!c This implementation requires rather big work areas in the interest of a 
!c faster algorithm (and is therefore suitable for Monte Carlo simulations). 
!c   
      subroutine clcdcf(ta,a,erra,na, &
!c                     i  i i    i
                       tb,b,errb,nb, &
!c                     i  i i    i
                       autocf,no0lag,unsmpl,used,unused,nwdata, &
!c                     i      i      i      o    o      i/o    
                       minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp, &
!c                     i/o    o    o     o     o   o      o      
                       inbin,minbin,nbins,mbins, &
!c                     o     i      i/o   i     
                       wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb, &
!c                     w   w   w  w  w     w     w   w   w      w
                       vtau)
!c                     w
      implicit none
      save
      integer i,j,k,n,used
      integer ibin
      integer na,nb,nbins,minbin,mbins,unused,minpts
      real ta(na),tb(nb)
      real a(na),b(nb)
      real erra(na),errb(nb)
      real xa,xb
      real wa(minbin),wb(minbin)
      real wta(minbin),wtb(minbin)
      real werra(minbin),werrb(minbin)
      integer wai(minbin,mbins),wbj(minbin,mbins)
      real vtau(minbin,mbins)
      real tdcf(mbins),sigtm(mbins),sigtp(mbins),dcf(mbins)
      real sigdcm(mbins),sigdcp(mbins)
      integer inbin(mbins)
      real expa,expb,vara,varb,vnorm,expbin,varbin,z,expz,sigz
      real sqrt,tanh,fishe,fishs,eps
      logical no0lag,nwdata,autocf,unsmpl
      logical wuseda(na),wusedb(nb)
      parameter (eps=1e-7)
!c   
!c If new data (i.e. NOT another Monte Carlo run) - allocate pairs to bins
!c   
      if (nwdata) then
!c Allocating the lags to the bins
         call alcbin(ta,na,tb,nb,wuseda,wusedb,no0lag,autocf,minpts, &
                    tdcf,sigtm,sigtp,wai,wbj,inbin,nbins,mbins,minbin, &
                    unsmpl, &
                    vtau)
!c Counting the unused points
         used = 0
         do ibin = 1,nbins
            used = used + inbin(ibin)
         enddo
         if (autocf) then
            if (no0lag) then 
               unused = na*(na-1)/2-used
            else
               unused = na*na/2-used
            endif
         else
            unused = na*nb-used
         endif
      endif
!c   
!c After allocating pairs to bins: calculating the dcf
!c   
      do ibin = 1,nbins
!c Collecting the points of the bin
         n = inbin(ibin)
         do k = 1,n
            i = wai(k,ibin)
            j = wbj(k,ibin)
            wta(k) = ta(i)
            wtb(k) = tb(j)
            wa(k) = a(i)
            wb(k) = b(j)
            werra(k) = erra(i)
            werrb(k) = errb(j)
         enddo
         expa = 0.0
         expb = 0.0
         vara = 0.0
         varb = 0.0
         do i = 1,n
            xa = wa(i)
            xb = wb(i)
            expa = expa+xa
            expb = expb+xb
            vara = vara+xa**2
            varb = varb+xb**2
         enddo
         expa = expa/n
         expb = expb/n
         vara = (vara-n*expa**2)/(n-1.)
         varb = (varb-n*expb**2)/(n-1.)
!c   
         expbin = 0.0
         varbin = 0.0
!c If normalization factor is 0 ...
         vnorm = vara*varb
         if (vnorm .le. 0.0) then
            vnorm = 0.0
            expbin = 0.0
         else
            vnorm = sqrt(vnorm)
            do i = 1,n
               expbin = expbin + (wa(i)-expa)*(wb(i)-expb)
            enddo
!c Dividing by (n-1) for an unbiased estimator of the correlation coefficient
!c cf Barlow / Statistics, p. 80
            expbin = expbin/vnorm/(n-1.)
         endif
         dcf(ibin) = expbin
!c   
!c Calculating the +/- 1 Sigma limits from Fisher's z
!c   
         if      (expbin .ge.  1.0-eps) then
            expbin =  1.0-eps
         else if (expbin .le. -1.0+eps) then 
            expbin = -1.0+eps
         endif
!c   
!c NOTE: This error estimation is by "bootstrapping": fishe & fishs give
!c     the true E(z) and S(z) when the TRUE correlation coefficient is 
!c     given. We are using the empirical r itself, similarily to the
!c     common poissonian estimate of n +/- sqrt(n)
!c   
        z = log((1.+expbin)/(1.-expbin))/2.0
        sigz = fishs(expbin,n)
        expz = fishe(expbin,n)
        sigdcm(ibin) = expbin-tanh(expz-sigz)
        sigdcp(ibin) = tanh(expz+sigz)-expbin
!c   
      enddo
      if (nwdata) then
         nwdata = .false.
      endif
      return
      end
!c/////////////////////////////////////////////////////////////////////////////

!c//////////////////////////////////////////////////////////////////////
      SUBROUTINE srtind (X, N, IPERM, KFLAG)
!c**LIBRARY   SLATEC
!c**CATEGORY  N6A1B, N6A2B
!c**TYPE      SINGLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!c**KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!c**AUTHOR  Jones, R. E., (SNLA)
!c         Rhoads, G. S., (NBS)
!c         Wisniewski, J. A., (SNLA)
!c**DESCRIPTION
!C
!c SRTIND returns the permutation vector IPERM generated by sorting
!c the array X and, optionally, rearranges the values in X.  X may
!c be sorted in increasing or decreasing order.  A slightly modified
!c quicksort algorithm is used.
!C
!c IPERM is such that X(IPERM(I)) is the Ith value in the rearrangement
!c of X.  IPERM may be applied to another array by calling IPPERM,
!c SPPERM, DPPERM or HPPERM.
!C
!c The main difference between SRTIND and its active sorting equivalent
!c SSORT is that the data are referenced indirectly rather than
!c directly.  Therefore, SRTIND should require approximately twice as
!c long to execute as SSORT.  However, SRTIND is more general.
!C
!c Description of Parameters
!c    X - input/output -- real array of values to be sorted.
!c        If ABS(KFLAG) = 2, then the values in X will be
!c        rearranged on output; otherwise, they are unchanged.
!c    N - input -- number of values in array X to be sorted.
!c    IPERM - output -- permutation array such that IPERM(I) is the
!c            index of the value in the original order of the
!c            X array that is in the Ith location in the sorted
!c            order.
!c    KFLAG - input -- control parameter:
!c          =  2  means return the permutation vector resulting from
!c                sorting X in increasing order and sort X also.
!c          =  1  means return the permutation vector resulting from
!c                sorting X in increasing order and do not sort X.
!c          = -1  means return the permutation vector resulting from
!c                sorting X in decreasing order and do not sort X.
!c          = -2  means return the permutation vector resulting from
!c                sorting X in decreasing order and sort X also.
!c**REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!c               for sorting with minimal storage, Communications of
!c               the ACM, 12, 3 (1969), pp. 185-187.
!C
!c   .. Scalar Arguments ..
      INTEGER KFLAG, N
!c   .. Array Arguments ..
      REAL X(*)
      INTEGER IPERM(*)
!c   .. Local Scalars ..
      REAL R, TEMP
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
!c   .. Local Arrays ..
      INTEGER IL(21), IU(21)
!c   .. Intrinsic Functions ..
      INTRINSIC ABS, INT

      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         print *,'SRTIND: N<1'
         stop 
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         print *,'SRTIND: Illegal sort control parameter'
         stop 
      ENDIF

!c   Initialize permutation vector

      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE

!c   Return if only one value is to be sorted

      IF (NN .EQ. 1) RETURN

!c   Alter array X to get decreasing order if needed

      IF (KFLAG .LE. -1) THEN
         DO 20 I=1,NN
            X(I) = -X(I)
   20    CONTINUE
      ENDIF

!c   Sort X only

      M = 1
      I = 1
      J = NN
      R = .375E0

   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF

   40 K = I

!c   Select a central element of the array and save it in location L

      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)

!c   If first element of array is greater than LM, interchange with LM

      IF (X(IPERM(I)) .GT. X(LM)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J

!c   If last element of array is less than LM, interchange with LM

      IF (X(IPERM(J)) .LT. X(LM)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)

!c      If first element of array is greater than LM, interchange
!c      with LM

         IF (X(IPERM(I)) .GT. X(LM)) THEN
            IPERM(IJ) = IPERM(I)
            IPERM(I) = LM
            LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT

!c   Find an element in the second half of the array which is smaller
!c   than LM

   60 L = L-1
      IF (X(IPERM(L)) .GT. X(LM)) GO TO 60

!c   Find an element in the first half of the array which is greater
!c   than LM

   70 K = K+1
      IF (X(IPERM(K)) .LT. X(LM)) GO TO 70

!c   Interchange these elements

      IF (K .LE. L) GO TO 50

!c   Save upper and lower subscripts of the array yet to be sorted

      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 90

!c   Begin again on another portion of the unsorted array

   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)

   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1

  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (X(IPERM(I)) .LE. X(LM)) GO TO 100
      K = I

  110 IPERM(K+1) = IPERM(K)
      K = K-1

      IF (X(LM) .LT. X(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100

!c   Clean up

  120 IF (KFLAG .LE. -1) THEN
         DO 130 I=1,NN
            X(I) = -X(I)
  130    CONTINUE
      ENDIF

!c   Rearrange the values of X if desired

     IF (KK .EQ. 2) THEN

!c      Use the IPERM vector as a flag.
!c      If IPERM(I) < 0, then the I-th value is in correct location

         DO 150 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0) THEN
               INDX = ISTRT
               INDX0 = INDX
               TEMP = X(ISTRT)
  140          IF (IPERM(INDX) .GT. 0) THEN
                  X(INDX) = X(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 140
               ENDIF
               X(INDX0) = TEMP
            ENDIF
  150    CONTINUE

!c      Revert the signs of the IPERM values

         DO 160 I=1,NN
            IPERM(I) = -IPERM(I)
  160    CONTINUE

      ENDIF

      RETURN
      END
!c/////////////////////////////////////////////////////////////////////////////
!//////////////////////////////////////////////////////////////////////////////

    subroutine parameters(shift,corr,nbin,nsamp,ma,t1,t2,t3,t4,t5,t6,t7,t8,r2,r4,r6,r8,loop2)

       real(DP), intent(in)  :: corr(5000)
       real , intent(in)     :: shift(5000)
       integer , intent(in)  :: nsamp,nbin,ma,loop2
       integer               :: FLAG1,FLAG3,temp2,FLAG5,FLAG7,FLAG9
       integer               :: i,j,k,l
       integer               :: index1,index2,index3,index4,index5,index6,index7,index8,index9
!       integer , dimension(100)  , intent(inout)   :: count1,count2,count3,count4,count5
       real(DP) ,dimension(20,1000) ,intent(out) :: r2,r4,r6,r8
       real(DP) ,dimension(20,1000)              :: r1,r3,r5,r7,r9
       real(DP) ,dimension(20,1000), intent(out) :: t1,t2,t3,t4,t5,t6,t7,t8
       real(DP)                                   :: t9(20,1000)
       real                         :: tempr2, tempr6, tempr4, tempr8

       character(80) outputFileName
       integer, parameter :: outFileUnit = 100000
 
       
       FLAG1 = 0
       FLAG3 = 0 
       FLAG5 = 0
       FLAG7 = 0
       FLAG9 = 0
       index1 = 0
       index2 = 0
       index3 = 0
       index4 = 0
       index5 = 0
       index6 = 0
       index7 = 0
       index8 = 0
       index9 = 0
! t1 - first zero crossing time
       t1(loop2,ma) = 0.0 
! t2 - first minimum reaching time
       t2(loop2,ma) = 0.0 
! t3 - second zero crossing time
       t3(loop2,ma) = 0.0 
! t4 - first maximum reaching time
       t4(loop2,ma) = 0.0
! t5 - third zero crossing time
       t5(loop2,ma) = 0.0 
! t6 - second minimum reaching time
       t6(loop2,ma) = 0.0 
! t7 - fourth zero crossing time
       t7(loop2,ma) = 0.0 
! t8 - second maximum reaching time
       t8(loop2,ma) = 0.0  
       t9(loop2,ma) = 0.0 
       r1(loop2,ma) = 0.0 
       r2(loop2,ma) = 0.0 
       r3(loop2,ma) = 0.0 
       r4(loop2,ma) = 0.0
       r5(loop2,ma) = 0.0 
       r6(loop2,ma) = 0.0 
       r7(loop2,ma) = 0.0 
       r8(loop2,ma) = 0.0
       r9(loop2,ma) = 0.0  
   
       temp2 = 1
 
!      write(*,*)"loop:",ma

! starting from tau=0, finding the first zero crossing point 

      ! write (outputFileName, '(a, I0, a)') 'Flag_trigger.dat'
      ! open(outFileUnit, file=fileplace//outputFileName,Access = 'append',status="unknown")
      ! write(outFileUnit,*)"loop2:,ma:",loop2,ma
           
        do 10 i = temp2,nbin

!         write(outFileUnit,*)"testing now at ",i,shift(i),corr(i)
       
         if(corr(i) .le. 0.0) then

           if(FLAG1 .eq. 0) then
 
            FLAG1 = 1
 
!            count1(loop2) = count1(loop2)+1
            
            t1(loop2,ma) = shift(i)

            r1(loop2,ma) = corr(i)

            index1 = i

            temp2 = i

            GOTO 10 
 
               
           elseif(FLAG5 .eq. 0 .and. FLAG3 .eq. 1) then

                 FLAG5 = 1

!                 count3(loop2) = count3(loop2)+1
             
                 t5(loop2,ma) = shift(i)

                 r5(loop2,ma) = corr(i)

                 index5 = i

                 temp2 = i

                 GOTO 10 


           elseif(FLAG9 .eq. 0 .and. FLAG7 .eq. 1) then

                 FLAG9 = 1

!                 count5(loop2) = count5(loop2)+1
             
                 t9(loop2,ma) = shift(i)

                 r9(loop2,ma) = corr(i)

                 index9 = i

                 temp2 = i

                 GOTO 10 


           end if


!------------------------------------------------------------
       elseif(corr(i) .ge. 0.0) then
          
             if(FLAG1 .eq. 1 .and. FLAG3 .eq. 0) then
               
                FLAG3 = 1

!                count2(loop2) = count2(loop2)+1
            
                t3(loop2,ma) = shift(i)

                r3(loop2,ma) = corr(i)
                
                index3 = i

                temp2 = i

                GOTO 10 

        
             elseif(FLAG5 .eq. 1 .and. FLAG7 .eq. 0)then
 
                 FLAG7 = 1

!                 count4(loop2) = count4(loop2)+1
             
                 t7(loop2,ma) = shift(i)

                 r7(loop2,ma) = corr(i)

                 index7 = i

                 temp2 = i

                 GOTO 10 



             end if

!----------------------

       end if

10     continue



           if(FLAG1 .eq. 1 .and. FLAG3 .eq. 1) then

             tempr2 = corr(index1)
              
              do j = index1,index3

               if(corr(j) .le. tempr2) then
       
                 index2 = j

                 tempr2 = corr(j)

                 t2(loop2,ma) = shift(j)

                 r2(loop2,ma) = corr(j)

               end if
    

              end do

           end if   
                                      
               
          if(FLAG5 .eq. 1 .and. FLAG7 .eq. 1) then


             tempr6 = corr(index5)
              
              do j = index5,index7

               if(corr(j) .le. tempr6) then
       
                 index6 = j
 
                 tempr6 = corr(j)

                 t6(loop2,ma) = shift(j)

                 r6(loop2,ma) = corr(j)

               end if

              end do  
                          

           end if


           if(FLAG3 .eq. 1 .and. FLAG5 .eq. 1) then

             tempr4 = corr(index3)
              
              do j = index3,index5

               if(corr(j) .ge. tempr4) then
       
                 index4 = j

                 tempr4 = corr(j)

                 t4(loop2,ma) = shift(j)

                 r4(loop2,ma) = corr(j)

               end if
    

              end do

          end if

        
          if(FLAG7 .eq. 1 .and. FLAG9 .eq. 1)then

             tempr8 = corr(index7)
              
              do j = index7,index9

               if(corr(j) .ge. tempr8) then
       
                 index8 = j

                 tempr8 = corr(j)

                 t8(loop2,ma) = shift(j)

                 r8(loop2,ma) = corr(j)

               end if
     

              end do 
                          

           end if




  
  end subroutine parameters


end program letsmakealc4


!cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine simlc(t,npoints,bsize,pfunc,idum, &
        typemod,tmult, dtime,dflux)

use double
!c  Uses myshape (function), four1 (Num.Rec. FFT),  gasdev

implicit none

real    :: tmax, pfunc(7),nyqf,nyqt,newt,var,tbin,utbin,minf
real    :: per(4200000),gasdev,pow,f
real    :: t(4200000), dtime(4200000),a,dflux(4200000)
real    :: tmult
real(DP) :: DTsim,bsize
real(8) :: myshape 
integer :: nf,i,ii,idum,j,npoints,k,typemod
integer :: NumSampd

!c---  Warning: the variable "nyqf" here is mis-used as nyquist time
        
utbin=0.

tbin=bsize

nyqf=tbin*2.0

print*,'SIMLC:npts,t(npts),t(1)', npoints,t(npoints),t(1)

print*,'SIMLC:DTsim = ',bsize

print*,'SIMLC:pfunc=', pfunc

tmax=nyqf*2**(int(log10((t(npoints)-t(1))/nyqf)/log10(2.0))+tmult)

print*,'tmax(s) = ', tmax

minf=1/tmax

nf=nint(tmax/nyqf)

print*,'SIMLC:Num of Freqs = ',nf

print*,'SIMLC:idum = ',idum

print*,'SIMLC:minf = ',minf

per(1)=0.0

per(2)=0.0

do i=3,((2*nf)+1),2

 f=(real(i-1)/2.0)*minf

!c          print*,'SIMLC:pfunc=',pfunc,typemod,f
 pow=myshape(f,pfunc,typemod)
!c          print*,'SIMLC:pow=  ',pow

 pow=pow*(2*nf/nyqf)

 per(i)=gasdev(idum)*sqrt(0.5*pow)

 per(i+1)=gasdev(idum)*sqrt(0.5*pow)

 ii=(4*nf)-(i-2)

 per(ii)=per(i)     

 per(ii+1)=-per(i+1)

enddo

per((2*nf)+1)=0.0  
        
call four1(per,2*nf,-1)     

k=0

j=1

do i=1,4*nf-1,2 

 k=k+1

 per(i)=per(i)/real(2*nf)

 newt=((i-1)/2)*(nyqf/2.0)

 dflux(k)=per(i)

 dtime(k)=newt
        
!c--   Trim from [2^X] down to Nsampd:

 if (k.ge.npoints) go to 444 

enddo

 444  print*,'k (trimmed to nsampd) = ',k

!c       open(unit=5,file='sampdk',status='unknown')
!c       do i=1,k

!c       write(5,*),dtime(i),dflux(i)
!c       end do
!c       close(unit=5)

end subroutine simlc
!c=================================gnu
!c=================================
subroutine four1(data,nn,isign)

integer :: isign,nn

real    :: data(2*nn)

integer :: i,istep,j,m,mmax,n

real    :: tempi,tempr

double precision :: theta,wi,wpi,wpr,wr,wtemp
      
n=2*nn

j=1

do 11 i=1,n,2

 if (j.gt.i) then

  tempr=data(j)

  tempi=data(j+1)

  data(j)=data(i)

  data(j+1)=data(i+1)

  data(i)=tempr

  data(i+1)=tempi

 endif

 m=n/2

 1 if ((m.ge.2).and.(j.gt.m)) then

    j=j-m

    m=m/2

    goto 1

   endif

   j=j+m
   
 11   enddo

   mmax=2
 2    if (n.gt.mmax) then
      istep=2*mmax
      theta=6.28318530717959d0/(isign*mmax)
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      do 13 m=1,mmax,2
        do 12 i=m,n,istep
          j=i+mmax
          tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
          tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
          data(j)=data(i)-tempr
          data(j+1)=data(i+1)-tempi
          data(i)=data(i)+tempr
          data(i+1)=data(i+1)+tempi
 12     enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
 13   enddo
      mmax=istep
      goto 2
      endif
      return
      end subroutine four1


      function gasdev(idum)
      integer idum
      real gasdev
      integer iset
      real fac,gset,rsq,v1,v2,ran2
      save iset,gset
      data iset/0/
      if (iset.eq.0) then
 1    v1=2.*ran2(idum)-1.
      v2=2.*ran2(idum)-1.
      rsq=v1**2+v2**2
      if (rsq.ge.1..or.rsq.eq.0.) goto 1
      fac=sqrt(-2.*log(rsq)/rsq)
      gset=v1*fac
      gasdev=v2*fac
      iset=1
      else
      gasdev=gset
      iset=0
      endif
      return
      end function gasdev


      function ran2(idum)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1, &
        ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211, &
        ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-6,rnmx=1.-eps)

      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
 11     continue
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if (iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)
      return
      end function ran2
        
!c===========================================================
!c====================================================
!c pfunc(1)=log10(break freq)     or fL (lin.,not log!) for Lorentz.
!c pfunc(2)=beta                  or Q                  for Lorentz.
!c pfunc(3)=gamma                 or zero               for Lorentz.
!c pfunc(7)=ampl                  or R                  for Lorentz.

!c "function myshape" edited to "real*8 function myshape" in ~2013
             real*8 function myshape(f,pfunc,typemod)
             real pfunc(7),ampl,f,beta,gamma,fcut,PPsn,temp
             real fL,Q,R
             integer typemod
                          
!c---- Inputs    f = frequency, in linear space e.g., 2.0e-6
!c----           pfunc array  (pfunc(1)=log10(fcut), etc.)

!c---  typemod = 1 for unbroken-PL PSD (I call it "UNBR") 
!c---  typemod = 2 for a sharply-broken PL PSD ("SING" = singly-broken PL)
!c --- typemod = 3 for a very slow bend PSD ("SLOW")
!c---  typemod = 4 for a bending power-law that's not as sharp as #3 ("SLOX")
!c---  typemod = 8 for Lorentzian (as defined in Belloni02, Nwk00, KP03; NOT Marko07)
             fcut = 10.**pfunc(1)
             beta = pfunc(2)
             gamma= pfunc(3)
             ampl=pfunc(4)
             R = pfunc(7)
             Q=pfunc(6)
             fL=pfunc(5)
             
!c---- PPsn: Power due to Poisson noise, in units of /Hz
!c---- Usually equal to [2*Src / (Src + Bkgd)^2] * [DTsamp/DTbin]
!c----- Src = NET count rate    bkgd = bkgd count rate
!c---  See Eqn A2 of S. Vaughan et al 2003 MNRAS 345 1271
!c---- PPsn = equals 2/Src for continuous sampling and zero bkgd
             PPsn = 0.0
             temp = 0.0
             
             if (typemod.eq.1) then
                myshape=  ampl* ((f/1.0E-6)**(-beta))  
                temp =  (2.*R*R*Q*fL/3.142)/( (fL**2.)+(2.*Q*(f-fL))**2.)        
             end if

             if (typemod.eq.2) then
                  if (f.le.fcut) then
                   myshape =  ampl * ((f/fcut)**(-gamma))
                  else
                   myshape =  ampl * ((f/fcut)**(-beta))
                  endif
             end if

       if (typemod.eq.3) then
       myshape=(ampl * f**(-gamma))/( ( 1 + f/fcut)**(beta-gamma))
       end if
       if (typemod.eq.4) then
       myshape=(ampl * f**(-gamma))/(  1 + (f/fcut)**(beta-gamma))
       end if

       if (typemod.eq.8) then
       myshape=(2.*R*R*Q*fL/3.142)/( (fL**2.)+(2.*Q*(f-fL))**2.)   
 !      temp = ampl* ((f/1.0E-6)**(-beta))      
  
       end if


       myshape=myshape+PPsn+temp
             return
             end function myshape

!c=================================

      FUNCTION ran3(idum)

      INTEGER idum
      INTEGER MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
      SAVE iff,inext,inextp,ma

      DATA iff /0/

      if(idum.lt.0.or.iff.eq.0)then
         iff=1
         mj=MSEED-iabs(idum)
         mj=mod(mj,MBIG)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
 12         continue
 13      continue
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return

      end function ran3

!c==================================================

      function initrand()

!C
!C   INITRAND:  generate a random number seed from the system clock
!C
!C------------------------------------------------------------------------------

!C   Internal variables

      integer(4) initrand

      initrand = time()
      end function initrand



!c   
!c Fisher's small sample approximation for E(z) (Kendall + Stuart Vol. 1 p.391)
!c   
      real function fishe(r,n)
      implicit none
      integer n
      real r
!c   
      fishe = log((1.+r)/(1.-r))/2.+   &
             r/2./(n-1.)*(1.+(5.+r**2)/4./(n-1.)+ &
                          (11.+2*r**2+3*r**4)/8./(n-1.)**2)
      return
      end
!c/////////////////////////////////////////////////////////////////////////////
!c   
!c Fisher's small sample approximation for s(z) (Kendall + Stuart Vol. 1 p.391)
!c   
      real function fishs(r,n)
      implicit none
      integer n
      real r,sqrt
!c   
      fishs = 1./(n-1.)*(1.+(4.-r**2)/2./(n-1.)+  &
                        (22.-6*r**2-3*r**4)/6./(n-1.)**2)
      fishs = sqrt(fishs)
      return
      end
!c//////////////////////////////////////////////////////////////////////
!c   Generating a standard Gaussian distr. using the Box-Muller method
!c   The existence of a subroutine rand(x) is assumed.
      real function rndnrm()
      implicit none
      save y1,y2
      real x1,x2,y1,y2,a
      real rand
      real pi2
      parameter (pi2 = 6.283185307)
      logical odd
      data odd/.true./
!c
      if (odd) then
         x1 = rand(0)
         x2 = rand(0)
         x2 = x2*pi2
         a =  sqrt(-2.*log(x1))
         y1 = a*cos(x2)
         y2 = a*sin(x2)
         rndnrm = y1
      else 
         rndnrm = y2
      endif 
      odd = .not. odd
      return 
      end

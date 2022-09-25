!gfortran rmspdm.f90 -o rmp.out

MODULE DOUBLE

IMPLICIT NONE 

integer,Parameter  :: DP=Kind(1d0) 

END MODULE DOUBLE

program rmsfreqplot

use double 

character(80) outputFileName
integer, parameter :: outFileUnit1 = 100000
integer, parameter :: outFileUnit2 = 10000
CHARACTER(*), PARAMETER       :: fileplace1 = "/work/saikruba/zdcf/Q_test/Q/M_data/aov_output/"
CHARACTER(*), PARAMETER       :: fileplace2 = "/work/saikruba/zdcf/R_test/R_noise/aov_outputR/"
REAL , DIMENSION(1520)        :: x,y,xr,yr,thetac,Psucc
REAL                          :: temp
REAL , DIMENSION(1520,200)    :: time,flux,timer,fluxr
integer                       :: num,k,l,count1,count2,start

 
N      = 302
num    = 0
thetac = 0.0
start  = 800
time   = 0.0
timer  = 0.0
flux   = 0.0
fluxr  = 0.0
! load data in 2D array for 100 values

do loop1 = 1501,1600

start = start+1
!print*,start
x= 0.0

y= 0.0

xr =0.0

yr =0.0

 num = num+1

 write (outputFileName, '(a, I0, a)') 'Q',loop1,'AOV.res'
 
  open(outFileUnit1, file=fileplace1//outputFileName,status="unknown")
   
   do i=1,N
   
     read(outFileUnit1,*)x(i),y(i)
   
   end do
 
  close(unit=outFileUnit1)

 write (outputFileName, '(a, I0, a)') 'R',start,'AOV.res'
 
  open(outFileUnit2, file=fileplace2//outputFileName,status="unknown")
   
   do i=1,N
   
     read(outFileUnit2,*)xr(i),yr(i)
   
   end do
 
  close(unit=outFileUnit2)

do loop2= 1,N

time(loop2,num)  = x(loop2)
flux(loop2,num)  = y(loop2)
timer(loop2,num) = xr(loop2)
fluxr(loop2,num) = yr(loop2)

end do
 
end do

! determine the critical theta minimun : va. sort the R array; b. find the 99 percent value (lower value)
thetac =0.0
  
 open(100, file=fileplace1//'rmsprob.dat',status="unknown")

do loop3 = 1, N
  
temp = 0.0

do k = 1 , num-1
   
  do l = k+1 , num
                 
   if (fluxr(loop3,k) > fluxr(loop3,l)) then
      
     temp = fluxr(loop3,k)
      
     fluxr(loop3,k) = fluxr(loop3,l)
      
     fluxr(loop3,l) = temp
            
    end if
    
   end do
   
  ! print*,fluxr(loop3,k)
   
 end do

 thetac(loop3) = fluxr(loop3,2)
! print*,thetac(loop3)

 count1 =0

 do loop4 = 1,num
 
 if(flux(loop3,loop4) < thetac(loop3)) then
   
   count1 = count1 + 1   
  
 end if 

 end do

 Psucc(loop3) = count1

 write(100,*)time(loop3,1),Psucc(loop3),thetac(loop3) 
 
end do
 
 close(100)

end program rmsfreqplot




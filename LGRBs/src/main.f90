MODULE MODELparameters
  IMPLICIT NONE
  INTEGER, PARAMETER :: nvar=4
  DOUBLE PRECISION, DIMENSION(nvar,nvar) :: MVNCOV
  DOUBLE PRECISION, DIMENSION(nvar,nvar) :: INVMVNCOV
  DOUBLE PRECISION, DIMENSION(nvar) :: X
  DOUBLE PRECISION :: dummy1     
  DOUBLE PRECISION :: dummy2     
  DOUBLE PRECISION :: dummy22    
  DOUBLE PRECISION :: dummy3     
  DOUBLE PRECISION :: dummy32    
  DOUBLE PRECISION :: aEL,bEL    
  DOUBLE PRECISION :: conepkzmean
  DOUBLE PRECISION :: loglumean  
  DOUBLE PRECISION :: logepkzmean
  DOUBLE PRECISION :: logeisomean
  DOUBLE PRECISION :: logt90zmean
END MODULE MODELparameters

MODULE Zparameters
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: zmin=0.1d0,zmax=1.d2
  DOUBLE PRECISION, PARAMETER :: z0=0.993d0,z1=3.8d0 
  DOUBLE PRECISION, PARAMETER :: g0=3.3d0,g1=0.0549d0,g2=-4.46d0
  DOUBLE PRECISION :: gamma1
  DOUBLE PRECISION :: gamma2
  DOUBLE PRECISION :: zplus1,logzplus1,lumdisterm,dvdzOzplus1
END MODULE Zparameters

MODULE COSMOparameters
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: omega_DE=0.7d0
  DOUBLE PRECISION, PARAMETER :: omega_DM=0.3d0
END MODULE COSMOparameters

MODULE constants
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: c=3.0d5,H=7.1d1 
  DOUBLE PRECISION, PARAMETER :: Mpc2cm=3.09d24  
  DOUBLE PRECISION, PARAMETER :: pi=dacos(-1.0d0)
  DOUBLE PRECISION, PARAMETER :: log10Mpc2cmSq4pi=dlog10(4.*pi)+2*dlog10(Mpc2cm)
  DOUBLE PRECISION, PARAMETER :: sqrt2=dsqrt(2.0d0)
  DOUBLE PRECISION, PARAMETER :: CoverH=C/H
  DOUBLE PRECISION, PARAMETER :: log10e=0.434294481903259
  DOUBLE PRECISION, PARAMETER :: sqrtPiO2=1.2533141373155
  DOUBLE PRECISION, PARAMETER :: sqrt2Pi=2.50662827463100
END MODULE constants

MODULE detection
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: signif=4.d0
  DOUBLE PRECISION :: stdevthresh
  DOUBLE PRECISION :: meanthresh
  DOUBLE PRECISION :: x0
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: dummy4
  DOUBLE PRECISION, PARAMETER :: logp1024min=4.92d0,logp1024max=6.318167895318538d0
  DOUBLE PRECISION, PARAMETER :: logepkmin=-2.915056638230699d0,logepkmax=5.4093868613659435d0
  DOUBLE PRECISION, PARAMETER :: logphminmaxdiff=logp1024max-logp1024min
  DOUBLE PRECISION :: logepkzmin,logepkzmax
  DOUBLE PRECISION :: log10pbol
END MODULE detection

MODULE GRBworld
  IMPLICIT NONE
  INTEGER, PARAMETER :: npar=16
  DOUBLE PRECISION, DIMENSION(npar) :: STCV
  DOUBLE PRECISION, DIMENSION(npar) :: STDEV
  DOUBLE PRECISION, DIMENSION(npar,npar) :: COV
  DOUBLE PRECISION, DIMENSION(npar,npar) :: COR
  DOUBLE PRECISION, DIMENSION(npar,2) :: SR
END MODULE GRBworld

MODULE OBSGRBDATA
  IMPLICIT NONE
  INTEGER, PARAMETER :: ndata=1366
  INTEGER :: idata
  TYPE :: GRBDATA
    INTEGER :: trigger
    DOUBLE PRECISION :: logpbol,logepk,logsbol,logt90,prob
  END TYPE GRBDATA
  TYPE(GRBDATA), DIMENSION(ndata) :: GRB
END MODULE OBSGRBDATA

MODULE Bandmodel
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: alpha=-1.1,beta=-2.3
END MODULE Bandmodel

PROGRAM TOYMODEL
  USE GRBworld
  USE OBSGRBDATA
  USE detection
  USE Zparameters, ONLY: z0,z1,g0,g1,g2,gamma1,gamma2
  USE constants, ONLY: log10Mpc2cmSq4pi
  USE Bandmodel
  IMPLICIT NONE
  INTEGER :: i,j
  INTEGER, PARAMETER :: iter=10**5,ireport=iter/10**3,idelay=10**5
  DOUBLE PRECISION, PARAMETER :: scalefactor=0.8*2.38**2/dble(npar)
  DOUBLE PRECISION, EXTERNAL :: LogLikelihood
  DOUBLE PRECISION :: zprob,PbolEpk2P1024ph_default
  DOUBLE PRECISION :: dummy
  CHARACTER(LEN=72) :: labels
  CHARACTER(LEN=30), DIMENSION(npar) :: LABEL
  gamma1=(1.+z0)**(g0-g1)
  gamma2=gamma1*(1.+z1)**(g1-g2)
  open(unit=11,file='../in/BATSE_1366_LGRB_P1024ph_Epk_Sch23ph.txt',status='old')
  open(unit=2121,file='BATSE_1366_LGRB_Pbol_Epk_Sbol(0.001,20000).txt')
  write(2121,'(8A30)') 'trigger','Log(Pbol[0.001,20000]Kev)','Log(Sbol[0.001,20000]Kev)','Log(Epk)',&
  'Log(EPR1024)','Log(EFR)','Log(FPR1024)','Log(T90)'
  read(11,*) labels
  do i=1,ndata
    read(11,*) GRB(i)%trigger,GRB(i)%logpbol,GRB(i)%logepk,GRB(i)%logsbol,GRB(i)%logt90
    GRB(i)%logpbol=GRB(i)%logpbol-PbolEpk2P1024ph_default(GRB(i)%logepk)
    GRB(i)%logsbol=GRB(i)%logsbol-PbolEpk2P1024ph_default(GRB(i)%logepk)
    write(2121,'(I30,7E30.6)') GRB(i)%trigger,GRB(i)%logpbol,GRB(i)%logsbol,GRB(i)%logepk,&
    GRB(i)%logepk-GRB(i)%logpbol,GRB(i)%logepk-GRB(i)%logsbol,GRB(i)%logsbol-GRB(i)%logpbol,&
    GRB(i)%logt90
  end do
  close(11)
  close(2121)
  open(unit=13,file='../in/iniparam.txt',status='old')
  read(13,*) labels
  read(13,*) labels,(STCV(j),j=1,npar)  
  write(*,*) 'MEAN of the parameters:'
  write(*,*) labels,(STCV(j),j=1,npar)  
  read(13,*) labels,(STDEV(j),j=1,npar)  
  write(*,*) 'Standard Deviations of the parameters:'
  write(*,*) labels,(STDEV(j),j=1,npar)  
  read(13,*) labels
  write(*,*) 'Correlation Matrix:'
  do i=1,npar
    read(13,*) LABEL(i),(COR(i,j),j=1,npar)
    write(*,*) LABEL(i),(COR(i,j),j=1,npar)
  end do
  write(*,*) 'Covariance Matrix:'
  do i=1,npar
    do j=1,npar
      COV(i,j)=STDEV(i)*STDEV(j)*COR(i,j)*scalefactor
    end do
    write(*,*) (COV(i,j),j=1,npar)
  end do
  close(13)
  SR(1,1)=-huge(scalefactor)
  SR(2,1)=-huge(scalefactor)
  SR(3,1)=-huge(scalefactor)
  SR(4,1)=-huge(scalefactor)
  SR(5,1)=-huge(scalefactor)
  SR(6,1)=-huge(scalefactor)
  SR(7,1)=-huge(scalefactor)
  SR(8,1)=-huge(scalefactor)
  SR(9,1)=-0.999d0 
  SR(10,1)=-0.999d0
  SR(11,1)=-0.999d0
  SR(12,1)=-0.999d0
  SR(13,1)=-0.999d0
  SR(14,1)=-0.999d0
  SR(15,1)=-huge(scalefactor)
  SR(16,1)=-huge(scalefactor)
  SR(1,2)=huge(scalefactor)
  SR(2,2)=huge(scalefactor)
  SR(3,2)=huge(scalefactor)
  SR(4,2)=huge(scalefactor)
  SR(5,2)=huge(scalefactor)
  SR(6,2)=huge(scalefactor)
  SR(7,2)=huge(scalefactor)
  SR(8,2)=huge(scalefactor)
  SR(9,2)=0.999d0 
  SR(10,2)=0.999d0
  SR(11,2)=0.999d0
  SR(12,2)=0.999d0
  SR(13,2)=0.999d0
  SR(14,2)=0.999d0
  SR(15,2)=huge(scalefactor)
  SR(16,2)=huge(scalefactor)
  call AMH(npar,iter,ireport,idelay,STCV,COV,SR,LogLikelihood,LABEL)
END PROGRAM TOYMODEL

INCLUDE 'LogLikelihood.f90'
INCLUDE 'AMH_TO.f90'
INCLUDE 'choldc.f90'
INCLUDE 'gasdevran2.f90'
INCLUDE 'MVNRND.f90'
INCLUDE 'ran2.f90'
INCLUDE 'samcovmat.f90'
INCLUDE 'posdef.f90'
INCLUDE 'erfcc.f90'
INCLUDE 'inversematrix.f90'
INCLUDE 'lubksb.f90'
INCLUDE 'ludcmp.f90'
INCLUDE 'determinant.f90'
INCLUDE 'qrombPbol.f90'
INCLUDE 'qrombEpk.f90'
INCLUDE 'PbolEpk2P1024ph_default.f90'
INCLUDE 'PbolEpk2P1024ph.f90'
INCLUDE 'qromo.f90'
INCLUDE 'midexp.f90'
INCLUDE 'polint.f90'
INCLUDE 'ldisWickram.f90'

!  Amir Shahmoradi, Suanday 11:19 PM, April 29, 2012, IFS, UTEXAS.
FUNCTION LogLikelihood(PARAM)
  USE MODELparameters
  USE Zparameters, ONLY: zmin,zmax
  USE constants, ONLY: sqrt2,sqrt2Pi
  USE GRBworld, ONLY: npar
  USE OBSGRBDATA
  USE detection
  IMPLICIT NONE
  INTEGER :: i,j,posdef
  DOUBLE PRECISION :: LogLikelihood,determinant
  DOUBLE PRECISION :: lumsigma,epkzsigma,conepkzsigma,normfac
  DOUBLE PRECISION, DIMENSION(npar) :: PARAM
  DOUBLE PRECISION, EXTERNAL :: modelintz,midexp,probatz
  DOUBLE PRECISION :: modelint
  loglumean=PARAM(1)
  logepkzmean=PARAM(2)
  logeisomean=PARAM(3)
  logt90zmean=PARAM(4)
  lumsigma=10.d0**PARAM(5)
  epkzsigma=10.d0**PARAM(6)
  bEL=PARAM(9)*epkzsigma/lumsigma
  aEL=PARAM(2)-PARAM(1)*bEL
  conepkzsigma=dsqrt(1.d0-PARAM(9)**2)*epkzsigma
  meanthresh=PARAM(15)
  stdevthresh=10.d0**PARAM(16)
  dummy2=sqrt2*lumsigma     
  dummy22=sqrt2*conepkzsigma
  dummy3=sqrt2Pi*lumsigma
  dummy32=sqrt2Pi*conepkzsigma
  dummy4=sqrt2*stdevthresh
  x0=meanthresh-signif*stdevthresh-logphminmaxdiff-logp1024min
  x1=meanthresh+signif*stdevthresh+logphminmaxdiff-logp1024min
  do i=1,nvar
    MVNCOV(i,i)=10.**PARAM(i+4)
  end do
  MVNCOV(1,2)=MVNCOV(1,1)*MVNCOV(2,2)*PARAM(9)
  MVNCOV(2,1)=MVNCOV(1,2)
  MVNCOV(1,3)=MVNCOV(1,1)*MVNCOV(3,3)*PARAM(10)
  MVNCOV(3,1)=MVNCOV(1,3)
  MVNCOV(1,4)=MVNCOV(1,1)*MVNCOV(4,4)*PARAM(11)
  MVNCOV(4,1)=MVNCOV(1,4)
  MVNCOV(2,3)=MVNCOV(2,2)*MVNCOV(3,3)*PARAM(12)
  MVNCOV(3,2)=MVNCOV(2,3)
  MVNCOV(2,4)=MVNCOV(2,2)*MVNCOV(4,4)*PARAM(13)
  MVNCOV(4,2)=MVNCOV(2,4)
  MVNCOV(3,4)=MVNCOV(3,3)*MVNCOV(4,4)*PARAM(14)
  MVNCOV(4,3)=MVNCOV(3,4)
  do i=1,nvar
    MVNCOV(i,i)=MVNCOV(i,i)**2
  end do
  INVMVNCOV=MVNCOV
  if (posdef(INVMVNCOV,nvar)==0) then
    write(*,*) 'Covariance matrix not positive definite..cycling..'
    LogLikelihood=-huge(LogLikelihood)
    RETURN
  end if
  normfac=determinant(nvar,nvar,MVNCOV)
  if (normfac<=0) then
    write(*,*) 'Covariance determinant is <=0',normfac,posdef(MVNCOV,nvar)
    do i=1,nvar
      write(*,*) (MVNCOV(i,j),j=1,nvar)
    end do
    STOP
  end if
  normfac=39.478417604357434*dsqrt(normfac)
  call inversematrix(nvar,nvar,MVNCOV,INVMVNCOV)
  call qromo(modelintz,zmin,zmax,modelint,midexp)
  if (modelint<=0.0d0) then
    write(*,*) 'model_integral (variable modelint in loglikelihood.f90) is zero:', modelint
    write(*,*) 'Press Enter to continue...'
    read(*,*)
  end if
  GRBPROBINTEG: do idata=1,ndata
    call qromo(probatz,zmin,zmax,GRB(idata)%prob,midexp)
    if (GRB(idata)%prob<0.0d0) then
      write(*,*) 'GRB_prob is negative:', GRB(idata)%prob
      write(*,*) 'Press Enter to continue...'
      read(*,*)
    end if
    GRB(idata)%prob=dlog(GRB(idata)%prob/(modelint*normfac))
  end do GRBPROBINTEG
  LogLikelihood=sum(GRB(1:ndata)%prob)
  if (LogLikelihood==0.d0) LogLikelihood=-huge(LogLikelihood)
END FUNCTION LogLikelihood

FUNCTION modelintz(z)
  USE MODELparameters, ONLY: dummy1,dummy2,loglumean
  USE COSMOparameters
  USE Zparameters
  USE constants, ONLY: sqrt2,sqrt2Pi,log10Mpc2cmSq4pi
  USE GRBworld, ONLY: npar
  USE OBSGRBDATA
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: z
  DOUBLE PRECISION :: lumdisMPc,lumdisMPcSq,erfcc
  DOUBLE PRECISION :: modelintz,ldisWickram    
  DOUBLE PRECISION, EXTERNAL :: Lisointegration
  zplus1=z+1.0d0
  logzplus1=dlog10(zplus1)
  lumdisMPc=ldisWickram(z)
  lumdisMPcSq=lumdisMPc*lumdisMPc
  lumdisterm=log10Mpc2cmSq4pi+dlog10(lumdisMPcSq)
  dummy1=lumdisterm-loglumean
  logepkzmin=logepkmin+logzplus1
  logepkzmax=logepkmax+logzplus1
  call qrombPbol(Lisointegration,x0+lumdisterm,x1+lumdisterm,modelintz)
  modelintz=modelintz+0.5d0*erfcc((x1+dummy1)/dummy2)
  dvdzOzplus1=lumdisMPcSq/(zplus1**3*dsqrt(omega_DM*zplus1**3+omega_DE))
  if (z<z0) then
    modelintz=modelintz*dvdzOzplus1*zplus1**g0
  else if (z>=z0 .and. z<z1) then
    modelintz=modelintz*gamma1*dvdzOzplus1*zplus1**g1
  else if (z>=z1) then
    modelintz=modelintz*gamma2*dvdzOzplus1*zplus1**g2
  else
    write(*,*) 'invalid redshift input in Likelihood integration, z:', z
    write(*,*) 'Press Enter to continue...'
    read(*,*)
  end if
END FUNCTION modelintz

FUNCTION Lisointegration(loglum)
  USE MODELparameters, ONLY: dummy2,dummy22,dummy3,aEL,bEL,conepkzmean,loglumean
  USE Zparameters, ONLY: lumdisterm
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: loglum,Lisointegration,erfcc,efficiency
  DOUBLE PRECISION, EXTERNAL :: EpkzProbGivenLiso
  conepkzmean=aEL+bEL*loglum
  log10pbol=loglum-lumdisterm
  call qrombEpk(EpkzProbGivenLiso,logepkzmin,logepkzmax,Lisointegration)
  efficiency=0.5d0+0.5d0*(1.d0-erfcc((logp1024min+log10pbol-meanthresh)/dummy4))
  Lisointegration=Lisointegration+efficiency*&
  (0.5d0*erfcc((logepkzmax-conepkzmean)/dummy22)+1.0d0-0.5d0*erfcc((logepkzmin-conepkzmean)/dummy22))
  Lisointegration=Lisointegration/(dummy3*dexp(((loglum-loglumean)/dummy2)**2))
END FUNCTION Lisointegration

FUNCTION EpkzProbGivenLiso(xx)
  USE MODELparameters, ONLY: dummy22,dummy32,conepkzmean
  USE Zparameters, ONLY: logzplus1
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: EpkzProbGivenLiso,erfcc,xx,efficiency,PbolEpk2P1024ph
  efficiency=0.5d0+0.5d0*(1.d0-erfcc((PbolEpk2P1024ph(xx-logzplus1,log10pbol)-meanthresh)/dummy4))
  EpkzProbGivenLiso=efficiency/(dummy32*dexp(((xx-conepkzmean)/dummy22)**2))
END FUNCTION EpkzProbGivenLiso

FUNCTION probatz(z)
  USE MODELparameters
  USE COSMOparameters
  USE Zparameters
  USE constants, ONLY: sqrt2,sqrt2Pi,log10Mpc2cmSq4pi
  USE GRBworld, ONLY: npar
  USE OBSGRBDATA
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: probatz,erfcc,PbolEpk2P1024ph,ldisWickram  
  DOUBLE PRECISION :: lumdisMPc,lumdisMPcSq
  DOUBLE PRECISION :: expterm,efficiency,logp1024ph
  DOUBLE PRECISION, INTENT(IN) :: z
  zplus1=z+1.0d0
  logzplus1=dlog10(zplus1)
  lumdisMPc=ldisWickram(z)
  lumdisMPcSq=lumdisMPc*lumdisMPc
  lumdisterm=log10Mpc2cmSq4pi+dlog10(lumdisMPcSq)
  dummy1=lumdisterm-loglumean
    X(1)=GRB(idata)%logpbol+dummy1
    X(2)=GRB(idata)%logepk+logzplus1-logepkzmean
    X(3)=GRB(idata)%logsbol+lumdisterm-logzplus1-logeisomean
    X(4)=GRB(idata)%logt90-0.66d0*logzplus1-logt90zmean
    expterm=0.5d0*dot_product(X,matmul(INVMVNCOV,X))
    if (GRB(idata)%logpbol<x0) then
      efficiency=0.d0
    else if (GRB(idata)%logpbol<x1) then
      logp1024ph=PbolEpk2P1024ph(GRB(idata)%logepk,GRB(idata)%logpbol)
      efficiency=0.5d0+0.5d0*(1.d0-erfcc((logp1024ph-meanthresh)/dummy4))
    else if (GRB(idata)%logpbol>=x1) then
      efficiency=1.d0
    else
      write(*,*) 'Wrong efficiency limit in LogLikelihood.f90, x0,x1:',x0,x1
      STOP
    end if
    dvdzOzplus1=lumdisMPcSq/(zplus1**3*dsqrt(omega_DM*zplus1**3+omega_DE))
    if (z<z0) then
      probatz=dvdzOzplus1*zplus1**g0*efficiency/dexp(expterm)
    else if (z>=z0 .and. z<z1) then
      probatz=gamma1*dvdzOzplus1*zplus1**g1*efficiency/dexp(expterm)
    else if (z>=z1) then
      probatz=gamma2*dvdzOzplus1*zplus1**g2*efficiency/dexp(expterm)
    else
      write(*,*) 'invalid redshift input in Likelihood integration, z:', z
      write(*,*) 'Press Enter to continue...'
      read(*,*)
    end if
END FUNCTION probatz

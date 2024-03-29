!  Amir Shahmoradi, March 22, 2012, 2:21 PM, IFS, UTEXAS
SUBROUTINE MVNRND(idum,n,MU,COV,X)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(INOUT) :: idum
  DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: MU
  DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: X
  DOUBLE PRECISION, DIMENSION(n,n), INTENT(INOUT) :: COV
  DOUBLE PRECISION, DIMENSION(n) :: DIAG,VECTOR
  DOUBLE PRECISION :: gasdevran2
  INTEGER :: i,j
  call choldc(COV,n,n,DIAG)
  do i=1,n
    VECTOR(i)=gasdevran2(idum)
    X(i)=VECTOR(i)*DIAG(i)
  end do
  do i=2,n
    X(i)=X(i)+dot_product(COV(i,1:i-1),VECTOR(1:i-1))
  end do
  X=X+MU
END SUBROUTINE MVNRND

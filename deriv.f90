MODULE MOCHILA 
IMPLICIT NONE 
REAL :: X, X0, X1, DX
INTEGER :: I
INTEGER,PARAMETER :: NPOINTS=1000
REAL, PARAMETER :: H=0.001

CONTAINS 

REAL FUNCTION F(X)
REAL, INTENT (IN) :: X 

F = SIN(X)

END FUNCTION

REAL FUNCTION DFX(X)
REAL, INTENT (IN) :: X 

DFX=( F(X+H)-F(X) ) / H

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE MOCHILA 


PROGRAM DERIVAD
USE MOCHILA 
IMPLICIT NONE

OPEN(FILE="OUTPUT.DAT",UNIT=10)

PRINT*,"INTRODUCE LOS LIMTES IZQ. Y DER. [A,B]"
READ*,X0,X1

DX = (X1-X0) / REAL(NPOINTS)

DO I=1, NPOINTS
    X=X0+DX*I
    WRITE(10,*) X, F(X), DFX(X)
ENDDO

CLOSE(10)
STOP
END PROGRAM DERIVAD
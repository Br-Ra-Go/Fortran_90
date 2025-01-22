
! EJERCICIO 5.A "BURDEN, NUMERICAL ANALYSIS" PAG. 173

MODULE MOCHILA 
    IMPLICIT NONE 
    REAL :: T , H, W, A, B, Y, L, M, P, Z
    INTEGER :: I
    INTEGER, PARAMETER :: NPOINTS=10
    REAL, PARAMETER :: DELTA=0.01
    
    CONTAINS 
    
        REAL FUNCTION F(T,Y)
        REAL, INTENT (IN) :: T, Y
    
            F =  -(Y+1)*(Y+3)
    
        ENDFUNCTION
    
        REAL FUNCTION YI(T)
        REAL, INTENT (IN) :: T
    
            YI = -3 + (2)/(1+EXP(-2*T))

        ENDFUNCTION
    
        REAL FUNCTION DFX(T)
        REAL, INTENT (IN) :: T 
    
            DFX=( YI(T+DELTA)-YI(T) ) / DELTA
    
        END FUNCTION
    
        REAL FUNCTION DDFX(T)
        REAL, INTENT (IN) :: T 
    
            DDFX=( DFX(T+DELTA)-DFX(T) ) / DELTA
    
        END FUNCTION
    
    END MODULE MOCHILA 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    PROGRAM ECU_DIF
    USE MOCHILA 
    IMPLICIT NONE 
    
    OPEN(FILE="OUTPUT_C_5.DAT", UNIT=11)
    
    ! INTRODUCE EL INTERVALO [A,B] Y LA CONDICIÓN INICIAL (Y0)
    
    A = 0.0
    B = 2.0
    Y = -2.0
    
    H = (B-A)/REAL(NPOINTS)
    T = A
    W = Y
    Z = YI(T)-W
    
    !DERIVADA DE LA ECUACIÓN DIFERNECIAL 
    
    L = -4-2*Y
    
    ! DERIVADA DE LA DERIVADA DE LA FUNCION ORIGINAL
    
    M = DDFX(B)
    
    WRITE(11,*) "# t_i               w_i            y_i=y(t_i)      ACTUAL_ERROR        ERROR BOUND"
    WRITE(11,*) T, W, YI(T), Z
    
    DO I = 1, NPOINTS
        W = W + (H*F(T,W))
        T = A + (I*H)
        Z = YI(T)-W
        P = ((H*M)/(2.0*L))*(EXP(L*(T-A))-1)
        WRITE(11,*) T, W, YI(T), Z, P 
    ENDDO
    
    
    !CLOSE(10)
    CLOSE(11)
    STOP
    END PROGRAM ECU_DIF
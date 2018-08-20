
     C B L A S _ D S Y M M  EXAMPLE PROGRAM

     INPUT DATA
       M=3  N=5
       ALPHA=  0.14  BETA=  1.62
       SIDE = CblasRight  UPLO = CblasLower  
       LAYOUT = CblasRowMajor  
       ARRAY A   LDA=5
            1.000  
           -2.320     4.330  
            3.000    -5.000     6.000  
            7.100     8.000     9.210    10.100  
           11.000    12.450    13.000   -14.700    15.000  
       ARRAY B   LDB=5
            1.000     2.000     3.000     4.000     5.000  
            1.000     2.000     3.000     4.000     5.000  
            1.000     2.000     3.000     4.000     5.000  
       ARRAY C   LDC=5
            1.100     1.100     1.100     1.100     1.100  
            1.100     1.100     1.100     1.100     1.100  
            1.100     1.100     1.100     1.100     1.100  

     OUTPUT DATA
       ARRAY C   LDC=5
           14.208    13.765    17.580     4.250    14.536  
           14.208    13.765    17.580     4.250    14.536  
           14.208    13.765    17.580     4.250    14.536  
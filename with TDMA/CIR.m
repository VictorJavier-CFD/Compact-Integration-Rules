%==========================================================================
%
% CIR:  Computes the integral of a function f = f(xi) using a Compact
% Integration Rule of m-th order accuracy. A compact integration rule
% solves a linear system in the form of
%                             _
%                            |     
%                          M | f dxi = Dxi Q f
%                           _|  
%                            
% where M and Q are banded matrices. Since the Compact Integration Rule
% implemented only uses 4th-, 6th-order accuracy, M is a tridiagonal 
% matrix. 
%
% In case of convergence issues (e.g. matrix M is not diagonal dominant), 
% use a successive over-relaxation method. If the linear system to solve 
% is written as Ax = b with A = D + L + U being D = diag(A), L = lower(A), 
% and U = upper(A), then this iterative method starts from a guess value, 
% i.e. x^0 = 0, and then apply the formula:
%
%      [  1            ]           1 - SOR          
%      [ --- D + L + U ] x^k = b + ------- D x^{k-1}    k = 1, 2, 3, ...   
%      [ SOR           ]             SOR                     
%
% where 0 < SOR <= 1 until convergence. The stopping criterion is based on
% reduce the residual so that resid_CCS <= residMIN_CCS. The residual is 
% defined as resid_CCS := ||b - A x^k||_2 / ||b - A x^0||_2 being || ||_2 
% the L2-norm. 
%
%
% Copyright © 2022 Llorente Lázaro, Víctor Javier
% Last Update: July 25, 2022
% Website: https://sites.google.com/view/vjllorente
% Contact: victor.javier.llorente@gmail.com
%
%--------------------------------------------------------------------------
%
% ------------
% Description:
% ------------
%
%  [ mcf(1) mef(1)                                                   ] [   vecint(1)    ]       [ qm(1,1) qm(1,2) qm(1,3) qm(1,4) qm(1,5)                                                                         ] [   vecfun(1)    ]
%  [ mwf(2) mcf(2) mef(2)                                            ] [   vecint(2)    ]       [ qm(2,1) qm(2,2) qm(2,3) qm(2,4)                                                                                 ] [   vecfun(2)    ]
%  [        mwf(3) mcf(3) mef(3)                                     ] [                ]       [         qm(3,2) qm(3,3) qm(3,4) qm(3,5)                                                                         ] [                ]
%  [               ...    ...    ...                                 ] [      ...       ] = Dxi [         ...     ...     ...     ...     ...                                                                     ] [      ...       ]
%  [                      ...    ...         ...                     ] [                ]       [                 ...     ...     ...     ...     ...                                                             ] [                ]
%  [                             mwf(nnod-2) mcf(nnod-2) mef(nnod-2) ] [ vecint(nnod-2) ]       [                                           qm(nnod-2,nnod-3) qm(nnod-2,nnod-2) qm(nnod-2,nnod-1) qm(nnod-2,nnod) ] [ vecfun(nnod-1) ]
%  [                                         mwf(nnod-1) mcf(nnod-1) ] [ vecint(nnod-1) ]       [                         qm(nnod-1,nnod-4) qm(nnod-1,nnod-3) qm(nnod-1,nnod-2) qm(nnod-1,nnod-1) qm(nnod-1,nnod) ] [  vecfun(nnod)  ]
%                                                                                           [   bcf(1)    ]
%                                                                                           [   bcf(2)    ]
%                                                                                           [             ]
%                                                                                         = [     ...     ]
%                                                                                           [             ]
%                                                                                           [ bcf(nnod-2) ]
%                                                                                           [ bcf(nnod-1) ]
%
% ----------------
% MATLAB function: 
% ----------------
%   vecint = CIR( vecfun, Dxi, mth )
%   Dependencies: 
%     1. TDMA.m
%
% ------
% INPUT:
% ------
%   Double, (nnod)x1 array        :: vecfun       - function to integrate
%   Integer                       :: Dxi          - interval size
%   Integer                       :: mth          - CIR order of accuracy
%
% -------
% OUTPUT:
% -------
%   Double, (nnod-1)x1 array      :: vecint       - integrals
%
% ----------------
% LOCAL VARIABLES:
% ----------------
%   Double, (nnod-1)x1 array      :: mcf          - diagonal of M matrix
%   Double, (nnod-2)x1 array      :: mwf, mef     - lower and upper diagonal of M matrix (Note: mwf(1) = mef(nnod-1) = 0)
%   Double, (nnod-1)x(nnod) array :: qm           - Q matrix
%   Double, (nnod-1)x1 array      :: bcf          - discrete source
%   Double                        :: SOR_CIR      - successive over-relaxation parameter
%   Double                        :: residMIN_CIR - minimal residual of SOR
%   Integer                       :: iterMAX      - maximum iteration of SOR
%
% -----------
% REFERENCES:
% -----------
%   1. Llorente, V.J. and Pascau, A. (2020) Compact Integration Rules as a 
%      quadrature method with some applications. Computers & Mathematics 
%      with Applications, 79(5), 1241-1265.
%
%==========================================================================
function vecint = CIR( vecfun, Dxi, mth )

    %% determines number of nodes
    nnod = length( vecfun );
    
    %% allocates
    qm = zeros(nnod - 1,nnod);
    mcf = zeros(nnod - 1,1); mef = mcf; mwf = mcf;
    vecint = zeros(nnod - 1,1);
    
    %% build up
    switch mth
        % Fourth order of accuracy
        case 4
            % Matrix M
            % upper diagonal
            mef(1) = 1; mef(2:nnod - 2) = 1 / 10;
            % principal diagonal
            mcf = ones(nnod - 1,1);
            % lower diagonal
            mwf(2:nnod - 2) = mef(2:nnod - 2); mwf(nnod - 1) = mef(1);
            % Matrix Q
            % left boundary
            qm(1,1) = 1 / 3; qm(1,2) = 4 / 3; qm(1,3) = 1 / 3;
            % inner
            for i = 2:nnod-2
                qm(i,i) = 3 / 5; qm(i,i + 1) = qm(i,i);
            end
            % right boundary
            qm(nnod - 1,nnod - 2) = qm(1,3); qm(nnod - 1,nnod - 1) = qm(1,2); qm(nnod - 1,nnod) = qm(1,1);
        % Sixth order of accuracy    
        case 6
            % Matrix M
            % upper diagonal
            mef(1) = 27 / 11; mef(2:nnod - 2) = 11 / 38;
            % principal diagonal
            mcf = ones(nnod - 1,1);
            % lower diagonal
            mwf(2:nnod - 2) = mef(2:nnod - 2); mwf(nnod - 1) = mef(1);
            % Matrix Q
            % left boundary
            qm(1,1) = 281 / 990; qm(1,2) = 1028 / 495; qm(1,3) = 196 / 165; qm(1,4) = -52 / 495; qm(1,5) = 1 / 90;
            % inner
            for i = 2:nnod-2
                qm(i,i - 1) = 3 / 38; qm(i,i) = 27 / 38; qm(i,i + 1) = qm(i,i); qm(i,i + 2) = qm(i,i - 1);
            end
            % right boundary
            qm(nnod - 1,nnod - 4) = qm(1,5); qm(nnod - 1,nnod - 3) = qm(1,4); qm(nnod - 1,nnod - 2) = qm(1,3); qm(nnod - 1,nnod - 1) = qm(1,2); qm(nnod - 1,nnod) = qm(1,1);
    end    
    
    %% discrete source
    bcf = Dxi * qm * vecfun;
    
    %% successive over-relaxation (SOR)
    SOR_CIR = 1.0;
    iterMAX = 1;
    residMIN_CIR = 1.0e-6;
    
    %% TDMA solver with SOR
    mcf_mod = mcf / SOR_CIR;
    matM = diag(mcf,0) + diag(mwf(2:nnod - 1),-1) + diag(mef(1:nnod - 2),1);
    for iter = 1:iterMAX
        bcf_mod = bcf + ( ( 1.0 - SOR_CIR ) / SOR_CIR ) * (mcf .* vecint);
        vecint = TDMA( mwf, mcf_mod, mef, bcf_mod );
        resid_CIR = norm( bcf - matM * vecint, 2 ) / norm( bcf, 2 );
        if resid_CIR <= residMIN_CIR
            break
        end
    end

end
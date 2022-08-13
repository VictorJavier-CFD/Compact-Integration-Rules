%==========================================================================
%
% CIR:  Computes the integral of a function f = f(xi) using a Compact
% Integration Rule of m-th order accuracy. A compact integration rule
% solves a linear system in the form of
%                            /     
%                          M | f dxi = Dxi Q f
%                            /  
% where M and Q are banded matrices. Since the Compact Integration Rule
% implemented only uses 4th-, 6th-order accuracy, M is a tridiagonal 
% matrix. 
%
% In case of convergence issues (e.g. matrix M is not diagonal dominant), 
% use a successive over-relaxation method. If the linear system to solve 
% is written as a_{ij} x_j = b_i, then this iterative method starts from
% a guess value, i.e. x_i^{0} = 0, and then apply the formula:   
%      a_{ii}                                   1 - SOR     
%      ------ x_i^k + sum(a_{ij} x_j^k) = b_i + ------- a_{ii} x_i^{k-1}
%       SOR                                       SOR
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

%**************************************************************************
%  Function to calculate an integral by CIR
%    y   -> |fdx vector 
%    var -> f vector 
%    nn  -> number of nodes
%    hh  -> interval length
%    th  -> order of CIR
%**************************************************************************

function y = CIR(var,nn,hh,th)
  %% Matrix
  Ms = sparse(nn-1,nn-1); Qs = sparse(nn-1,nn);
  %% Builder
  switch th
    case 4
      % Left boundary
      i = 1;
      C = i; E = C + 1; EE = C + 2;
      Ms(C,C) = 1; Ms(C,E) = 1;
      Qs(C,C) = 1/3; Qs(C,E) = 4/3; Qs(C,EE) = 1/3; 
      % Inner
      for i = 2:nn-2 
        C = i; W = C - 1; E = C + 1;
        Ms(C,W) = 1/10; Ms(C,C) = 1; Ms(C,E) = 1/10;
        Qs(C,W) = 0; Qs(C,C) = 3/5; Qs(C,E) = 3/5;    
      end
      % Right boundary
      i = nn-1; 
      C = i; W = C - 1; E = C + 1;
      Ms(C,W) = 1; Ms(C,C) = 1;
      Qs(C,W) = 1/3; Qs(C,C) = 4/3; Qs(C,E) = 1/3;        
    case 6
      % Left boundary
      i = 1;
      C = i; E = C + 1; EE = C + 2; EEE = C + 3; EEEE = C + 4;
      Ms(C,C) = 1; Ms(C,E) = 27/11;
      Qs(C,C) = 281/990; Qs(C,E) = 1028/495; Qs(C,EE) = 196/165; Qs(C,EEE) = -52/495; Qs(C,EEEE) = 1/90; 
      % Inner
      for i = 2:nn-2 
        C = i; W = C - 1; E = C + 1; EE = C + 2;
        Ms(C,W) = 11/38; Ms(C,C) = 1; Ms(C,E) = 11/38;
        Qs(C,W) = 3/38; Qs(C,C) = 27/38; Qs(C,E) = 27/38; Qs(C,EE) = 3/38;    
      end
      % Right boundary
      i = nn-1; 
      C = i; W = C - 1; WW = C - 2; WWW = C - 3; E = C + 1;
      Ms(C,W) = 27/11; Ms(C,C) = 1;
      Qs(C,WWW) = 1/90; Qs(C,WW) = -52/495; Qs(C,W) = 196/165; Qs(C,C) = 1028/495; Qs(C,E) = 281/990;              
  end 
  %% Solver
  M = full(Ms); Q = full(Qs);
  y = linsolve(M,hh*Q*var')';
end
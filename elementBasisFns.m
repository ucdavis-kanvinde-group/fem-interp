function [ N, dN, d2N ] = elementBasisFns( xi_, lmtype )
%Element Basis Functions
% Given the parent space coordinates, returns the values of the basis
% functions, first derivatives, and second derivatives.
%
%Input --
%   xi_:    vector of the parent space coordinates [xi, eta, zta]
%
%   lmtype: string identifier of the element type (i.e. 'QUAD4'). If you
%           specify a 2D element (like QUAD4) then the 3rd component of the
%           xi_ vector will not be used. There is no error or warning if
%           the 3rd component of the xi_ vector is nonzero.
%
%Output --
%   N:   (n x 1) vector of the basis functions for the n nodes
%
%   dN:  (n x 3) are the first derivs, where dN(a,j) is the derivative of 
%        the a-th basis function w.r.t. the j-th direction (i.e. N_a,j)
%
%   d2N: (n x 3 x 3) are the second derivs, where d2N(a,j,k) is the second 
%        derivative of the a-th basis function w.r.t. the j-th and k-th 
%        direction (i.e. N_a,jk)
%
%Note that the d2N output is not memory efficient, since there is symmetry
%that is not accounted for. For example, N_1,ij = N_1,ji ... meaning that
%d2N(1,i,j) = d2N(1,j,i). To be memory efficient, you would normally pack
%this matrix into a rank-2 (n x 6) matrix, where the ordering is:
%[N_a,11 N_a,22 N_a,33 N_a,23 N_a,13 N_a,12].
%
%However, I have instead opted for simplicity in output access, so that it
%provides more easily readable code.
%


switch upper(lmtype)
    case 'QUAD4'
        [N, dN, d2N] = basisQ4( xi_(1), xi_(2) );
    case 'QUAD8'
        [N, dN, d2N] = basisQ8( xi_(1), xi_(2) );
    case 'BRK8'
        [N, dN, d2N] = basisBRK8( xi_(1), xi_(2), xi_(3) );
    case 'BRK20'
        [N, dN, d2N] = basisBRK20( xi_(1), xi_(2), xi_(3) );
    otherwise
        error('undefined element type')
end

return;
end


function [ N, dN, d2N ] = basisQ4( xi, eta )
%QUAD4 Element
%
%N (4x1) are the basis functions, [N1, N2, N3, N4]
%
%dN (4x3) are the first derivs, where dN(a,j) is the a-th basis function
%w.r.t. the j-th direction (i.e. N_a,j)
%
%d2N (4x3x3) are the second derivs, where d2N(a,j,k) is the a-th basis
%function w.r.t. the j-th and k-th direction (i.e. N_a,jk)
%

%basis functions
N = [0.25*(1 - xi)*(1 - eta);
	 0.25*(1 + xi)*(1 - eta);
     0.25*(1 + xi)*(1 + eta);
     0.25*(1 - xi)*(1 + eta)];

%first derivatives
%dN is setup so that dN(a,j) is the a-th basis function, partial derivative
%w.r.t. the j-th direction (i.e. N_a,j)
dN = [-0.25*(1 - eta), -0.25*(1 - xi), 0;
       0.25*(1 - eta), -0.25*(1 + xi), 0;
       0.25*(1 + eta),  0.25*(1 + xi), 0;
      -0.25*(1 + eta),  0.25*(1 - xi), 0];

%d2N is setup so that d2N(:,:,k) is dN(:,:)_,k
d2N(:,:,1) = [0,  0.25, 0;
              0, -0.25, 0;
              0,  0.25, 0;
              0, -0.25, 0];

d2N(:,:,2) = [ 0.25, 0, 0;
              -0.25, 0, 0;
               0.25, 0, 0;
              -0.25, 0, 0];

d2N(:,:,3) = zeros(4,3);

return;
end

function [ N, dN, d2N ] = basisQ8( xi, eta )
%QUAD8 Serendipity Element
%
%N (8x1) are the basis functions, [N1, N2, N3, N4, ...]
%
%dN (8x3) are the first derivs, where dN(a,j) is the a-th basis function
%w.r.t. the j-th direction (i.e. N_a,j)
%
%d2N (8x3x3) are the second derivs, where d2N(a,j,k) is the a-th basis
%function w.r.t. the j-th and k-th direction (i.e. N_a,jk)
%

%basis functions: from "Concepts and applications of finite element
%analysis" / Robert D. Cook [et al.].--4th ed.
N5 = 0.5*(1 - xi*xi)*(1 - eta);
N6 = 0.5*(1 + xi)*(1 - eta*eta);
N7 = 0.5*(1 - xi*xi)*(1 + eta);
N8 = 0.5*(1 - xi)*(1 - eta*eta);

N1 = 0.25*(1 - xi)*(1 - eta) - 0.5*(N8 + N5); 
N2 = 0.25*(1 + xi)*(1 - eta) - 0.5*(N5 + N6);
N3 = 0.25*(1 + xi)*(1 + eta) - 0.5*(N6 + N7);
N4 = 0.25*(1 - xi)*(1 + eta) - 0.5*(N7 + N8);

N = [N1; N2; N3; N4; N5; N6; N7; N8];

%first derivatives
%dN# is setup so that, e.g, dN1(i) = N_1,i
dN5 = [   - xi*(1 - eta), -0.5*(1 - xi*xi)];
dN6 = [0.5*(1 - eta*eta),    -eta*(1 + xi)];
dN7 = [   - xi*(1 + eta),  0.5*(1 - xi*xi)];
dN8 = [-0.5*(1 - eta*eta),   -eta*(1 - xi)];

dN1 = [-0.25*(1 - eta) - 0.5*(dN8(1) + dN5(1)), ...
       -0.25*(1 -  xi) - 0.5*(dN8(2) + dN5(2))];
dN2 = [ 0.25*(1 - eta) - 0.5*(dN5(1) + dN6(1)), ...
       -0.25*(1 +  xi) - 0.5*(dN5(2) + dN6(2))];
dN3 = [ 0.25*(1 + eta) - 0.5*(dN6(1) + dN7(1)), ...
        0.25*(1 +  xi) - 0.5*(dN6(2) + dN7(2))];
dN4 = [-0.25*(1 + eta) - 0.5*(dN7(1) + dN8(1)), ...
        0.25*(1 -  xi) - 0.5*(dN7(2) + dN8(2))];


%dN is setup so that dN(a,j) is the a-th basis function, partial derivative
%w.r.t. the j-th direction (i.e. N_a,j)
dN  = [dN1, 0; dN2, 0; dN3, 0; dN4, 0; dN5, 0; dN6, 0; dN7, 0; dN8, 0];

%second derivatives
%d2N is setup so that d2N(:,:,k) is dN(:,:)_,k
d2N(:,:,1) = [0.5*(1 - eta),  0.25 - 0.5*(eta + xi), 0;
              0.5*(1 - eta), -0.25 - 0.5*(xi - eta), 0;
              0.5*(1 + eta),  0.25 + 0.5*(xi + eta), 0;
              0.5*(1 + eta), -0.25 - 0.5*(eta - xi), 0;
                - (1 - eta),                     xi, 0;
                          0,                  - eta, 0;
                - (1 + eta),                   - xi, 0;
                          0,                    eta, 0];

d2N(:,:,2) = [ 0.25 - 0.5*(eta + xi), 0.5*(1 - xi), 0;
              -0.25 - 0.5*(xi - eta), 0.5*(1 + xi), 0;
               0.25 + 0.5*(xi + eta), 0.5*(1 + xi), 0;
              -0.25 - 0.5*(eta - xi), 0.5*(1 - xi), 0;
                                  xi,            0, 0;
                               - eta,   - (1 + xi), 0;
                                - xi,            0, 0;
                                 eta,   - (1 - xi), 0];

d2N(:,:,3) = zeros(8,3);

return;
end

function [ N, dN, d2N ] = basisBRK8( xi, eta, zta )
%BRK8: 8 Node "Brick" 3D Continuum element
%
%N (8x1) are the basis functions, [N1, N2, N3, N4, ...]
%
%dN (8x3) are the first derivs, where dN(a,j) is the a-th basis function
%w.r.t. the j-th direction (i.e. N_a,j)
%
%d2N (8x3x3) are the second derivs, where d2N(a,j,k) is the a-th basis
%function w.r.t. the j-th and k-th direction (i.e. N_a,jk)
%


%basis functions: from "Finite Element Procedures in Engineering 
%Analysis" / Klaus-Jurgen Bath (1982)
N1 = 0.125*(1 - xi)*(1 - eta)*(1 - zta);
N2 = 0.125*(1 + xi)*(1 - eta)*(1 - zta);
N3 = 0.125*(1 + xi)*(1 + eta)*(1 - zta);
N4 = 0.125*(1 - xi)*(1 + eta)*(1 - zta);
N5 = 0.125*(1 - xi)*(1 - eta)*(1 + zta);
N6 = 0.125*(1 + xi)*(1 - eta)*(1 + zta);
N7 = 0.125*(1 + xi)*(1 + eta)*(1 + zta);
N8 = 0.125*(1 - xi)*(1 + eta)*(1 + zta);

N = [N1; N2; N3; N4; N5; N6; N7; N8];

%first derivatives
%dN# is setup so that, e.g, dN1(i) = 8 * N_1,i (factor of 1/8 left out)
dN1 = [-(1 - eta)*(1 - zta), -(1 - xi)*(1 - zta), -(1 - xi)*(1 - eta)];
dN2 = [ (1 - eta)*(1 - zta), -(1 + xi)*(1 - zta), -(1 + xi)*(1 - eta)];
dN3 = [ (1 + eta)*(1 - zta),  (1 + xi)*(1 - zta), -(1 + xi)*(1 + eta)];
dN4 = [-(1 + eta)*(1 - zta),  (1 - xi)*(1 - zta), -(1 - xi)*(1 + eta)];
dN5 = [-(1 - eta)*(1 + zta), -(1 - xi)*(1 + zta),  (1 - xi)*(1 - eta)];
dN6 = [ (1 - eta)*(1 + zta), -(1 + xi)*(1 + zta),  (1 + xi)*(1 - eta)];
dN7 = [ (1 + eta)*(1 + zta),  (1 + xi)*(1 + zta),  (1 + xi)*(1 + eta)];
dN8 = [-(1 + eta)*(1 + zta),  (1 - xi)*(1 + zta),  (1 - xi)*(1 + eta)];

%dN is setup so that dN(a,j) is the a-th basis function, partial derivative
%w.r.t. the j-th direction (i.e. N_a,j)
dN = 0.125 * [dN1; dN2; dN3; dN4; dN5; dN6; dN7; dN8]; %need to include 1/8

%second derivatives
%d2N is setup so that d2N(:,:,k) is dN(:,:)_,k
d2N(:,:,1) = [ 0,  0.125*(1 - zta),  0.125*(1 - eta);
               0, -0.125*(1 - zta), -0.125*(1 - eta);
               0,  0.125*(1 - zta), -0.125*(1 + eta);
               0, -0.125*(1 - zta),  0.125*(1 + eta);
               0,  0.125*(1 + zta), -0.125*(1 - eta);
               0, -0.125*(1 + zta),  0.125*(1 - eta);
               0,  0.125*(1 + zta),  0.125*(1 + eta);
               0, -0.125*(1 + zta), -0.125*(1 + eta)];

d2N(:,:,2) = [ 0.125*(1 - zta),  0,  0.125*(1 - xi);
              -0.125*(1 - zta),  0,  0.125*(1 + xi);
               0.125*(1 - zta),  0, -0.125*(1 + xi);
              -0.125*(1 - zta),  0, -0.125*(1 - xi);
               0.125*(1 + zta),  0, -0.125*(1 - xi);
              -0.125*(1 + zta),  0, -0.125*(1 + xi);
               0.125*(1 + zta),  0,  0.125*(1 + xi);
              -0.125*(1 + zta),  0,  0.125*(1 - xi)];

d2N(:,:,3) = [ 0.125*(1 - eta),  0.125*(1 - xi),  0;
              -0.125*(1 - eta),  0.125*(1 + xi),  0;
              -0.125*(1 + eta), -0.125*(1 + xi),  0;
               0.125*(1 + eta), -0.125*(1 - xi),  0;
              -0.125*(1 - eta), -0.125*(1 - xi),  0;
               0.125*(1 - eta), -0.125*(1 + xi),  0;
               0.125*(1 + eta),  0.125*(1 + xi),  0;
              -0.125*(1 + eta),  0.125*(1 - xi),  0];

return;
end

function [ N, dN, d2N ] = basisBRK20( xi, eta, zta )
%BRK20: 20 Node "Brick" 3D Continuum element
%
%N (20x1) are the basis functions, [N1, N2, N3, N4, ...]
%
%dN (20x3) are the first derivs, where dN(a,j) is the a-th basis function
%w.r.t. the j-th direction (i.e. N_a,j)
%
%d2N (20x3x3) are the second derivs, where d2N(a,j,k) is the a-th basis
%function w.r.t. the j-th and k-th direction (i.e. N_a,jk)
%


%basis functions: from "Finite Element Procedures in Engineering 
%Analysis" / Klaus-Jurgen Bath (1982).

%first, compute the basis functions for midside nodes
N9  = 0.25*(1 -   xi*xi)*(1 - eta)*(1 - zta);
N10 = 0.25*(1 - eta*eta)*(1 +  xi)*(1 - zta);
N11 = 0.25*(1 -   xi*xi)*(1 + eta)*(1 - zta);
N12 = 0.25*(1 - eta*eta)*(1 -  xi)*(1 - zta);
N13 = 0.25*(1 -   xi*xi)*(1 - eta)*(1 + zta);
N14 = 0.25*(1 - eta*eta)*(1 +  xi)*(1 + zta);
N15 = 0.25*(1 -   xi*xi)*(1 + eta)*(1 + zta);
N16 = 0.25*(1 - eta*eta)*(1 -  xi)*(1 + zta);
N17 = 0.25*(1 - zta*zta)*(1 -  xi)*(1 - eta);
N18 = 0.25*(1 - zta*zta)*(1 +  xi)*(1 - eta);
N19 = 0.25*(1 - zta*zta)*(1 +  xi)*(1 + eta);
N20 = 0.25*(1 - zta*zta)*(1 -  xi)*(1 + eta);

%next, compute the basis functions for the corner nodes. They are the same
%as a BRK8 element, except now they're superposed by the midside nodes
N1 = 0.125*(1 - xi)*(1 - eta)*(1 - zta) - 0.5*(N9  + N12 + N17);
N2 = 0.125*(1 + xi)*(1 - eta)*(1 - zta) - 0.5*(N9  + N10 + N18);
N3 = 0.125*(1 + xi)*(1 + eta)*(1 - zta) - 0.5*(N10 + N11 + N19);
N4 = 0.125*(1 - xi)*(1 + eta)*(1 - zta) - 0.5*(N11 + N12 + N20);
N5 = 0.125*(1 - xi)*(1 - eta)*(1 + zta) - 0.5*(N13 + N16 + N17);
N6 = 0.125*(1 + xi)*(1 - eta)*(1 + zta) - 0.5*(N13 + N14 + N18);
N7 = 0.125*(1 + xi)*(1 + eta)*(1 + zta) - 0.5*(N14 + N15 + N19);
N8 = 0.125*(1 - xi)*(1 + eta)*(1 + zta) - 0.5*(N15 + N16 + N20);

N = [N1; N2; N3; N4; N5; N6; N7; N8; ...
     N9; N10; N11; N12; N13; N14; N15; N16; N17; N18; N19; N20];

%
% first derivatives
%
%dN# is setup so that, e.g, dN1(i) = 8 * N_1,i

%first derivs of midside nodes
dN9  = [ -0.5*xi*(1 - eta)*(1 - zta), -0.25*(1 - xi*xi)*(1 - zta), ...
                                      -0.25*(1 - xi*xi)*(1 - eta)];
dN10 = [ 0.25*(1 - eta*eta)*(1 - zta), -0.5*eta*(1 + xi)*(1 - zta), ...
                                       -0.25*(1 - eta*eta)*(1 + xi)];
dN11 = [ -0.5*xi*(1 + eta)*(1 - zta), 0.25*(1 - xi*xi)*(1 - zta), ...
                                     -0.25*(1 - xi*xi)*(1 + eta)];  
dN12 = [-0.25*(1 - eta*eta)*(1 - zta), -0.5*eta*(1 - xi)*(1 - zta), ...
                                      -0.25*(1 - eta*eta)*(1 - xi)];
dN13 = [ -0.5*xi*(1 - eta)*(1 + zta), -0.25*(1 - xi*xi)*(1 + zta), ...
                                       0.25*(1 - xi*xi)*(1 - eta)];
dN14 = [ 0.25*(1 - eta*eta)*(1 + zta), -0.5*eta*(1 + xi)*(1 + zta), ...
                                       0.25*(1 - eta*eta)*(1 + xi)];
dN15 = [ -0.5*xi*(1 + eta)*(1 + zta), 0.25*(1 - xi*xi)*(1 + zta), ...
                                      0.25*(1 - xi*xi)*(1 + eta)];
dN16 = [-0.25*(1 - eta*eta)*(1 + zta), -0.5*eta*(1 - xi)*(1 + zta), ...
                                       0.25*(1 - eta*eta)*(1 - xi)];
dN17 = [-0.25*(1 - zta*zta)*(1 - eta), -0.25*(1 - zta*zta)*(1 - xi), ...
                                       -0.5*zta*(1 - xi)*(1 - eta)];
dN18 = [ 0.25*(1 - zta*zta)*(1 - eta), -0.25*(1 - zta*zta)*(1 + xi), ...
                                       -0.5*zta*(1 + xi)*(1 - eta)];
dN19 = [ 0.25*(1 - zta*zta)*(1 + eta), 0.25*(1 - zta*zta)*(1 + xi), ...
                                       -0.5*zta*(1 + xi)*(1 + eta)];
dN20 = [-0.25*(1 - zta*zta)*(1 + eta), 0.25*(1 - zta*zta)*(1 - xi), ...
                                       -0.5*zta*(1 - xi)*(1 + eta)];

%first derivs of vertex nodes, without 1/8 term or midside superposition
dN1 = [-(1 - eta)*(1 - zta), -(1 - xi)*(1 - zta), -(1 - xi)*(1 - eta)];
dN2 = [ (1 - eta)*(1 - zta), -(1 + xi)*(1 - zta), -(1 + xi)*(1 - eta)];
dN3 = [ (1 + eta)*(1 - zta),  (1 + xi)*(1 - zta), -(1 + xi)*(1 + eta)];
dN4 = [-(1 + eta)*(1 - zta),  (1 - xi)*(1 - zta), -(1 - xi)*(1 + eta)];
dN5 = [-(1 - eta)*(1 + zta), -(1 - xi)*(1 + zta),  (1 - xi)*(1 - eta)];
dN6 = [ (1 - eta)*(1 + zta), -(1 + xi)*(1 + zta),  (1 + xi)*(1 - eta)];
dN7 = [ (1 + eta)*(1 + zta),  (1 + xi)*(1 + zta),  (1 + xi)*(1 + eta)];
dN8 = [-(1 + eta)*(1 + zta),  (1 - xi)*(1 + zta),  (1 - xi)*(1 + eta)];

%superpose midside nodes into vertex nodes, and include 1/8 term
dN1 = 0.125*dN1 - 0.5*(dN9  + dN12 + dN17);
dN2 = 0.125*dN2 - 0.5*(dN9  + dN10 + dN18);
dN3 = 0.125*dN3 - 0.5*(dN10 + dN11 + dN19);
dN4 = 0.125*dN4 - 0.5*(dN11 + dN12 + dN20);
dN5 = 0.125*dN5 - 0.5*(dN13 + dN16 + dN17);
dN6 = 0.125*dN6 - 0.5*(dN13 + dN14 + dN18);
dN7 = 0.125*dN7 - 0.5*(dN14 + dN15 + dN19);
dN8 = 0.125*dN8 - 0.5*(dN15 + dN16 + dN20);

%dN is setup so that dN(a,j) is the a-th basis function, partial derivative
%w.r.t. the j-th direction (i.e. N_a,j)
dN = [ dN1;  dN2;  dN3;  dN4;  dN5;  dN6;  dN7;  dN8; ...
    dN9; dN10; dN11; dN12; dN13; dN14; dN15; dN16; dN17; dN18; dN19; dN20];

%
% second derivatives
%
%the derivatives are getting quite complicated now, so I chose to introduce
%intermediate derivative vectors. Since the 2nd derivs of an individual 
%basis function results in a symmetric matrix, I chose to store them as 
%vectors similar to how a Cauchy stress tensor is stored in vector form.
%Specifically:
%   d2N#(i) = N_#,ii {for i = 1:3}
%   d2N#(4) = N_#,23 = N_#,32
%   d2N#(5) = N_#,13 = N_#,31
%   d2N#(6) = N_#,12 = N_#,21

%2nd derivs of midside nodes
d2N9  = [ -0.5*(1 - eta)*(1 - zta),       0,               0,  ...
         0.25*(1 - xi*xi), 0.5*xi*(1 - eta), 0.5*xi*(1 - zta)];

d2N10 = [0,            -0.5*(1 + xi)*(1 - zta),                 0, ...
         0.5*eta*(1 + xi), -0.25*(1 - eta*eta), -0.5*eta*(1 - zta)];

d2N11 = [-0.5*(1 + eta)*(1 - zta),         0,                0, ...
         -0.25*(1 - xi*xi), 0.5*xi*(1 + eta), -0.5*xi*(1 - zta)];

d2N12 = [0,           -0.5*(1 - xi)*(1 - zta),                0, ...
         0.5*eta*(1 - xi), 0.25*(1 - eta*eta), 0.5*eta*(1 - zta)];
     
d2N13 = [-0.5*(1 - eta)*(1 + zta),          0,               0, ...
         -0.25*(1 - xi*xi), -0.5*xi*(1 - eta), 0.5*xi*(1 + zta)];

d2N14 = [0,            -0.5*(1 + xi)*(1 + zta),                 0, ...
         -0.5*eta*(1 + xi), 0.25*(1 - eta*eta), -0.5*eta*(1 + zta)];

d2N15 = [-0.5*(1 + eta)*(1 + zta),         0,                0, ...
         0.25*(1 - xi*xi), -0.5*xi*(1 + eta), -0.5*xi*(1 + zta)];

d2N16 = [0,             -0.5*(1 - xi)*(1 + zta),                0, ...
         -0.5*eta*(1 - xi), -0.25*(1 - eta*eta), 0.5*eta*(1 + zta)];

d2N17 = [   0,      0,                   -0.5*(1 - xi)*(1 - eta), ...
         0.5*zta*(1 - xi), 0.5*zta*(1 - eta), 0.25*(1 - zta*zta)];

d2N18 = [   0,      0,                     -0.5*(1 + xi)*(1 - eta), ...
         0.5*zta*(1 + xi), -0.5*zta*(1 - eta), -0.25*(1 - zta*zta)];

d2N19 = [   0,      0,                    -0.5*(1 + xi)*(1 + eta), ...
         -0.5*zta*(1 + xi), -0.5*zta*(1 + eta), 0.25*(1 - zta*zta)];

d2N20 = [   0,      0,                     -0.5*(1 - xi)*(1 + eta), ...
         -0.5*zta*(1 - xi), 0.5*zta*(1 + eta), -0.25*(1 - zta*zta)];

%2nd derivs of vertex nodes, without midside superposition
d2N1 = [0, 0, 0,  0.125*(1 - xi),  0.125*(1 - eta),  0.125*(1 - zta)];
d2N2 = [0, 0, 0,  0.125*(1 + xi), -0.125*(1 - eta), -0.125*(1 - zta)];
d2N3 = [0, 0, 0, -0.125*(1 + xi), -0.125*(1 + eta),  0.125*(1 - zta)];
d2N4 = [0, 0, 0, -0.125*(1 - xi),  0.125*(1 + eta), -0.125*(1 - zta)];
d2N5 = [0, 0, 0, -0.125*(1 - xi), -0.125*(1 - eta),  0.125*(1 + zta)];
d2N6 = [0, 0, 0, -0.125*(1 + xi),  0.125*(1 - eta), -0.125*(1 + zta)];
d2N7 = [0, 0, 0,  0.125*(1 + xi),  0.125*(1 + eta),  0.125*(1 + zta)];
d2N8 = [0, 0, 0,  0.125*(1 - xi), -0.125*(1 + eta), -0.125*(1 + zta)];

%add in the midside superposition
d2N1 = d2N1 - 0.5*(d2N9  + d2N12 + d2N17);
d2N2 = d2N2 - 0.5*(d2N9  + d2N10 + d2N18);
d2N3 = d2N3 - 0.5*(d2N10 + d2N11 + d2N19);
d2N4 = d2N4 - 0.5*(d2N11 + d2N12 + d2N20);
d2N5 = d2N5 - 0.5*(d2N13 + d2N16 + d2N17);
d2N6 = d2N6 - 0.5*(d2N13 + d2N14 + d2N18);
d2N7 = d2N7 - 0.5*(d2N14 + d2N15 + d2N19);
d2N8 = d2N8 - 0.5*(d2N15 + d2N16 + d2N20);

%need to assemble these derivatives into the proper format for output.
%d2N is setup so that d2N(:,:,k) is dN(:,:)_,k

d2N(:,:,1) = [ d2N1(1),  d2N1(6),  d2N1(5);
               d2N2(1),  d2N2(6),  d2N2(5);
               d2N3(1),  d2N3(6),  d2N3(5);
               d2N4(1),  d2N4(6),  d2N4(5);
               d2N5(1),  d2N5(6),  d2N5(5);
               d2N6(1),  d2N6(6),  d2N6(5);
               d2N7(1),  d2N7(6),  d2N7(5);
               d2N8(1),  d2N8(6),  d2N8(5);
               d2N9(1),  d2N9(6),  d2N9(5);
              d2N10(1), d2N10(6), d2N10(5);
              d2N11(1), d2N11(6), d2N11(5);
              d2N12(1), d2N12(6), d2N12(5);
              d2N13(1), d2N13(6), d2N13(5);
              d2N14(1), d2N14(6), d2N14(5);
              d2N15(1), d2N15(6), d2N15(5);
              d2N16(1), d2N16(6), d2N16(5);
              d2N17(1), d2N17(6), d2N17(5);
              d2N18(1), d2N18(6), d2N18(5);
              d2N19(1), d2N19(6), d2N19(5);
              d2N20(1), d2N20(6), d2N20(5)];

d2N(:,:,2) = [ d2N1(6),  d2N1(2),  d2N1(4);
               d2N2(6),  d2N2(2),  d2N2(4);
               d2N3(6),  d2N3(2),  d2N3(4);
               d2N4(6),  d2N4(2),  d2N4(4);
               d2N5(6),  d2N5(2),  d2N5(4);
               d2N6(6),  d2N6(2),  d2N6(4);
               d2N7(6),  d2N7(2),  d2N7(4);
               d2N8(6),  d2N8(2),  d2N8(4);
               d2N9(6),  d2N9(2),  d2N9(4);
              d2N10(6), d2N10(2), d2N10(4);
              d2N11(6), d2N11(2), d2N11(4);
              d2N12(6), d2N12(2), d2N12(4);
              d2N13(6), d2N13(2), d2N13(4);
              d2N14(6), d2N14(2), d2N14(4);
              d2N15(6), d2N15(2), d2N15(4);
              d2N16(6), d2N16(2), d2N16(4);
              d2N17(6), d2N17(2), d2N17(4);
              d2N18(6), d2N18(2), d2N18(4);
              d2N19(6), d2N19(2), d2N19(4);
              d2N20(6), d2N20(2), d2N20(4)];

d2N(:,:,3) = [ d2N1(5),  d2N1(4),  d2N1(3);
               d2N2(5),  d2N2(4),  d2N2(3);
               d2N3(5),  d2N3(4),  d2N3(3);
               d2N4(5),  d2N4(4),  d2N4(3);
               d2N5(5),  d2N5(4),  d2N5(3);
               d2N6(5),  d2N6(4),  d2N6(3);
               d2N7(5),  d2N7(4),  d2N7(3);
               d2N8(5),  d2N8(4),  d2N8(3);
               d2N9(5),  d2N9(4),  d2N9(3);
              d2N10(5), d2N10(4), d2N10(3);
              d2N11(5), d2N11(4), d2N11(3);
              d2N12(5), d2N12(4), d2N12(3);
              d2N13(5), d2N13(4), d2N13(3);
              d2N14(5), d2N14(4), d2N14(3);
              d2N15(5), d2N15(4), d2N15(3);
              d2N16(5), d2N16(4), d2N16(3);
              d2N17(5), d2N17(4), d2N17(3);
              d2N18(5), d2N18(4), d2N18(3);
              d2N19(5), d2N19(4), d2N19(3);
              d2N20(5), d2N20(4), d2N20(3)];

return;
end
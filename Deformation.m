function [z] = Deformation(g, mu, E, I, dx, f)
% Computes the deformations by Az = b equivalence 
% 
% Inputs: 
%         g - gravitational constant
%         mu - mass per unit length of the beam
%         E - young's modulus
%         I - second moment of area of the beam
%         dx - spacing between the points 
%         f - M by 1 column vector representing the point load
% Output: 
%         z - M by 1 column vector of deformations 

    M    = 101;   % Design doc says to let M = 101.
    
    % b vector
    % Solution vector equations 7-11.
    B = [[0; 0];...
        ((-mu*g)-f([3:M-2],:)) .*( (dx^4)/(E*I) );...
        [0; 0]];

    % set up the matrix to compute z, corresponds to the coefficients  
    % of the equations 7-11.
    A = zeros(M);
    A(1,:)   = [1 zeros([1,M-1])];
    A(2,:)   = [-3, 4, -1, zeros([1, M-3])];
    for m = 3:M-2
        A(m,m-2:m+2) = [1, -4, 6, -4, 1];
    end
    A(M-1,:) = [zeros([1, M-3]), 1, -4, 3];
    A(M,:)   = [zeros([1, M-1]) 1];

    % Solve for z by left dividing the coefficient matrix A by the 
    % solution vector v.
    z = A\B;

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Deformation.m
% EAS230
% Professor Sabato, Professor Ali
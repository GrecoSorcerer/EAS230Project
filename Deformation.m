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
    b = [[0 0], ones([M-4,1]).*( ( (-mu*g)-f ).*( (dx.^4)./(E*I) ) ) , [0 0]]';
    
    % set up the matrix to compute z, corresponds to the coefficients  
    % of the equations 7-11.
    A = [[1, zeros([1,M-1])]; [-3, 4, -1, zeros([1, M-3])]; ones([[M-4,M]]).*[1, -4, 6, -4, 1, zeros([M-5, 1])]; [zeros([1, M-3]), 1, -4, 3]; [zeros([1, M-1]), 1]];
    
    % Solve for z by left dividing the coefficient matrix A by the 
    % solution vector v.
    z = A\b;

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Deformation.m
% EAS230
% Professor Sabato, Professor Ali
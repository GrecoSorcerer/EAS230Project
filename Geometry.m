function [a, b, I] = Geometry(cross_section, area, orientation)
% Computes the dimensions and second moment of area of a beam given a
% user-defined cross-section, area, and orientation
%
% Inputs:
%         cross_section - type of beam desired by the user
%         area - desired geometric area of the beam's cross-section
%         orientation - desired orientation (horizontal(2) or vertical(1))
% Output:
%         a - dimension 1 of the beam
%         b - dimension 2 of the beam
%         I - measure of efficiency of a beam's resistance to deformation


    switch cross_section
        case 1
    
            b = sqrt(area/pi);
            a = 2*b;
            if (orientation == 1 || orientation == 2)
                I = (pi/4)*b^4;
            else
                a = -1;
                b = -1;
                I = -1;
            end
    
        case 2
    
            b = sqrt(area/6);
            a = 6*b;
    
            if orientation == 1
                I = (b*a^3)/12;
            elseif orientation == 2
                I = (a*b^3)/12;
            else
                a = -1;
                b = -1;
                I = -1;
            end
    
        case 3
    
            b = sqrt(area/18);
            a = 6*b;
    
            if orientation == 1
                I = ( ( a * ( ( a + 2*b )^3 ) - (a - b) * (a^3) ) / 12 );
            elseif orientation == 2
                I = ((a*b^3)/12) + ((b*a^3)/6);
            else
                a = -1;
                b = -1;
                I = -1;
            end
    
        case 4
    
            b = sqrt(area/12);
            a = 6*b;
    
            if orientation == 1
                y_c = (a+b) - ( ( b * ( ( a + b )^2 ) + (a - b) * b^2 ) ...
                    / (4*a*b) );
    
                I = ( ( b * ( ( a + b )^3 ) + (a - b) * b^3 ) / 3 ) - (2*a*b)*(a + b - y_c)^2;
    
            elseif orientation == 2
                I = ((a*b^3)+(b*a^3))/12;
    
            else
                a = -1;
                b = -1;
                I = -1;
    
            end
    
        case 5
    
            b = sqrt(area/11);
            a = 6*b;
            if (orientation == 1 || orientation == 2)
                I = ( ( b * ( (5*a^2) - (5*a*b) + b^2 ) * ((a^2) - (a*b) + (b^2)) ) ) / ( 12*(2*a - b) ) ;
            else
                a = -1;
                b = -1;
                I = -1;
            end
    
        otherwise
            a = -1;
            b = -1;
            I = -1;
    
    end

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Geometry.m
% EAS230
% Professor Sabato, Professor Ali
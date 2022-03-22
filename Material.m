function [rho, E, sigma] = Material(mat)
% Returns the material properties; density, Youngs modulus, and Modulus 
% of Rupture, for a given material.
%   
% Inputs: 
%       mat - numeric input 1-7 to select material properties
% Output: 
%       rho - Density in Specific Gravity ( )
%         E - Youngs Modulus (Pa)
%     sigma - Modulus of Rupture (Pa)
% 
%  < mat options >
%  1 - White oak       2 - Western white pine  3 - Red maple
%  4 - Particle board  5 - Plywood             6 - Aluminum
%  7 - Steel (Structural Steel ATSM-A36)


    filename = "MaterialData.dat";
    data = importdata(filename)
    
    % Logical indexing to fix table
    dataYS = data(:,3:4);
    dataYS(isnan(data(:,3:4))) = 0;
    data = [data(:,1:2) dataYS(:,1)+dataYS(:,2)];
    
    % Select properties for desired material
    rho   = data(mat, 1);       % density of the beam
    E     = data(mat, 2)*(1e9); % in Pa (from GPa)
    sigma = data(mat, 3)*(1e6); % in Pa (from MPa)

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Material.m
% EAS230
% Professor Sabato, Professor Ali
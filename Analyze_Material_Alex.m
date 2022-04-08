% DECLARATIONS_____________________________________________________________

MAXTRIES = 3;

cs_area = 0.01;  % units in m^2
L  = 3;     % units in m
S_f     = 4;     % unitless
g       = 9.81;  % units in m/s^2
M       = 101;   % unitless



% CROSS SECTION INPUT______________________________________________________

cross_section = Print_CS_Menu(MAXTRIES);

if (cross_section  == -1)
    error("Too many invalid entries!")
end



%ORIENTATION INPUT_________________________________________________________

orientation = Print_O_Menu;

if (orientation  == -1)
    error("Too many invalid entries!")
end



% CALLING Geometry.m and Material.m________________________________________

[a, b, I] = Geometry(cross_section, cs_area, orientation);
[rho_oak, E_oak, sigma_oak] = Material(1);        % White Oak
[rho_pine, E_pine, sigma_pine] = Material(2);     % Western white pine
[rho_maple, E_maple, sigma_maple] = Material(3);  % Red maple
[rho_pb, E_pb, sigma_pb] = Material(4);           % Particle board
[rho_ply, E_ply, sigma_ply] = Material(5);        % Plywood
[rho_Al, E_Al, sigma_Al] = Material(6);           % Aluminum
[rho_St, E_St, sigma_St] = Material(7);           % Steel

%COMPUTATING RECCOMENDED MAX LOAD__________________________________________

deltaX = (L / M - 1);

% initializing the point load array
f_m = zeros([1,M]);
m = 1:M; % indexing array
% computing the location of the point load
f_m(m == (M-1)/2) = F/dx;
f_m = f_m';

% computing the max load

sigmaMax = ( max(a,b) .* ( (f_m .* L) ./ (4*I) )) ./ safety_factor;

% computing mu for all 7 materials______________________________________

mu_oak = rho_oak*cs_area; 
mu_pine = rho_pine*cs_area;
mu_maple = rho_maple*cs_area;
mu_pb = rho_pb*cs_area;
mu_ply = rho_ply*cs_area;
mu_AL = rho_Al*cs_area;
mu_St = rho_St*cs_area;

% creating deformation vectors for all 7 materials______________________

[z_oak] = Deformation(g,mu_oak,E,I,dx,f_m); 
[z_pine] = Deformation(g,mu_pine,E,I,dx,f_m);
[z_maple] = Deformation(g,mu_maple,E,I,dx,f_m);
[z_pb] = Deformation(g,mu_pb,E,I,dx,f_m);
[z_ply] = Deformation(g,mu_ply,E,I,dx,f_m);
[z_Al] = Deformation(g,mu_Al,E,I,dx,f_m);
[z_St] = Deformation(g,mu_St,E,I,dx,f_m);

% creating table of deformation vectors_________________________________

deformTable = [z_oak' z_pine' z_maple' z_pb' z_ply' z_Al' z_St'];

% creating table of materials, max load, failure load, max deformation, weight

materials = ['White Oak';'Western White Pine';'Red Maple';'Particle board';'Plywood';'Aluminum';'Steel'];
max_loads = sigmaMax';
failure_loads = -1; %didnt know how to compute
max_deformations = -1; %didnt know how to compute
weights = -1; %didnt know how to compute

printedTable = [];

% 1) saving deformTable according to chosen orientation and cross_section
% 2) printing printedTable
% 3) plotting z versus x for each material
% 4) plotting max load vs. Young;s modulus for each material on log-log
%__________________________________________________________________________

switch orientation 
    case 1
        switch cross_section 
            case 1
                save Circular_Vertical.dat deformTable -ascii
                % do 2-4 listed above for each case
            case 2
                save Rectangular_Vertical.dat deformTable -ascii
            case 3
                save I-Beam_Verticle.dat deformTable -ascii
            case 4
                save T-Beam_Verticle.dat deformTable -ascii
            case 5
                save L-Beam_Verticle.dat deformTable -ascii
        end

    case 2
        switch cross_section 
            case 1
                save Circular_Horizontal.dat deformTable -ascii
            case 2
                save Rectangular_Horizontal.dat deformTable -ascii
            case 3
                save I-Beam_Horizontal.dat deformTable -ascii
            case 4
                save T-Beam_Horizontal.dat deformTable -ascii
            case 5
                save L-Beam_Horizontal.dat deformTable -ascii
        end
end

%HELPER FUNCTIONS__________________________________________________________

function [cross_section] = Print_CS_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response
    % Returns -1 if exceeds tries
    % Be sure to catch the -1 as an error
    
    if (tries == 0)
        cross_section = -1;
        return;
    end
    
    disp('Choose a cross-section');
    disp('    1 - Circular');
    disp('    2 - Rectangular');
    disp('    3 - I-Beam');
    disp('    4 - T-Beam');
    disp('    5 - L-Beam');

    op = input('Option: ','s');

    if ( isnumeric(op) && ( (op >0) && (op <=5) ) )
        cross_section = op;
    else
        cross_section = Print_CS_Menu(tries - 1);
    end
end

function [orientation] = Print_O_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response.
    % Returns -1 if exceeds tries
    % Be sure to catch the -1 as an error
        
    if (tries == 0)
        orientation = -1;
        return;
    end

    disp('Choose an orientation');
    disp('    1 - Vertical');
    disp('    2 - Horizontal');

    op = input('Option: ', 's');

    if ( isnumeric(op) && ( (op >0) && (op <=2) ) )
        orientation = op;
    else
        orientation = Print_O_Menu(tries - 1); 
    end

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Analyze_Material.m
% EAS230
% Professor Sabato, Professor Ali
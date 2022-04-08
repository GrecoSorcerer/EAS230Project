% DECLARATIONS_____________________________________________________________

MAXTRIES = 3;

cs_area        = 0.01;  % units in m^2
L              = 3;     % units in m
safety_factor  = 4;     % unitless
g              = 9.81;  % units in m/s^2
M              = 101;   % unitless
F              = 0;     % Load in N/m, no default given, assuming 1.

% CROSS SECTION INPUT______________________________________________________

cross_section = Print_CS_Menu(MAXTRIES);

% error handling
if (cross_section  == -1)
    error("Too many invalid entries!")
end

%ORIENTATION INPUT_________________________________________________________

orientation = Print_O_Menu(MAXTRIES);

% error handling
if (orientation  == -1)
    error("Too many invalid entries!")
end

%MATERIAL INPUT____________________________________________________________

material = Print_M_Menu(MAXTRIES);

% error handling
if (material  == -1)
    error(-1,"Too many invalid entries!")
end



%COMPUTATIONS______________________________________________________________

[rho, E, sigma] = Material(material);

[a, b, I] = Geometry(cross_section, cs_area, orientation);

% Compute the change in x
dx  = L / (M -1);

% Initialize point load array
f_m = zeros([1,M]);
m = 1:M; % indexing array
% Compute the point load.
f_m(m == (M-1)/2) = F/dx;
f_m = f_m';

% Compute Theoretical Max force with safety factor
sigmaMax = ( max(a,b) .* ( (f_m .* L) ./ (4*I) )) ./ safety_factor;

mu = rho*cs_area;

[z] = Deformation(g,mu,E,I,dx,f_m);

x = ((m-1)./(M-1)).*L;

%plot(x,z)

%HELPER FUNCTIONS__________________________________________________________ 

function [cross_section] = Print_CS_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response.
    % Returns -1 if exceeds tries
    % Be sure to catch the -1 as an error
    
    % Recursive Exit Condition
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

    op = input('Option: ');
    
    % Recursive loop condition
    if isempty(op)
        op = -1;
    end
    if ( isnumeric(op)) && ( (op >0) && (op <=5) ) 
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
    
    % Recursive Exit Condition
    if (tries == 0)
        orientation = -1;
        return;
    end

    disp('Choose an orientation');
    disp('    1 - Vertical');
    disp('    2 - Horizontal');

    op = input('Option: ');
    if isempty(op)
        op = -1;
    end
    % Recursive loop condition
    if ( isnumeric(op)) && ( (op >0) && (op <=2) ) 
        orientation = op;
    else
        orientation = Print_O_Menu(tries - 1); 
    end

end

function [material] = Print_M_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response.
    % Returns -1 if exceeds tries
    % Be sure to catch the -1 as an error
    
    % Recursive Exit Condition
    if (tries == 0)
        material = -1;
        return;
    end

    disp('Choose a material');
    disp('    1 - White Oak');
    disp('    2 - Western White Pine');
    disp('    3 - Red Maple');
    disp('    4 - Particle board');
    disp('    5 - Plywood');
    disp('    6 - Aluminum');
    disp('    7 - Steel');
    
    op = input('Option: ');
    if isempty(op)
        op = -1;
    end
    % Recursive loop condition
    if ( isnumeric(op) ) && ( (op >0) && (op <=7) ) 
        material = op;
    else
        material = Print_M_Menu(tries - 1);
    end

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Single_Case.m
% EAS230
% Professor Sabato, Professor Ali
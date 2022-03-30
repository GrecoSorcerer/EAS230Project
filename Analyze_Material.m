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



%COMPUTATING RECCOMENDED MAX LOAD__________________________________________

deltaX = (L / M - 1);



%HELPER FUNCTIONS__________________________________________________________

function [cross_section] = Print_CS_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response.
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
% DECLARATIONS_____________________________________________________________

MAXTRIES       = 3;
ORIENTATION    = containers.Map([1,2],{'vertical','horizontal'});
CROSS_SECTION  = containers.Map([1,2,3,4,5], ...
                                {'Circular', 'Rectangular', 'I-Beam',...
                                 'T-Beam',   'L-Beam'});
MATERIAL       = containers.Map([1,2,3,4,5,6,7], ...
                                {'White Oak', 'Western White Pine', ...
                                 'Red Maple', 'Particle board',     ...
                                 'Plywood', 'Aluminum',             ...
                                 'Steel'});

cs_area        = 0.01;  % units in m^2
L              = 3;     % units in m
safety_factor  = 4;     % unitless
g              = 9.81;  % units in m/s^2
M              = 101;   % unitless

% CROSS SECTION INPUT______________________________________________________

cross_section = Print_CS_Menu(MAXTRIES);

if (cross_section  == -1)
    error("Too many invalid entries!")
end


% ORIENTATION INPUT________________________________________________________

orientation = Print_O_Menu(MAXTRIES);

if (orientation  == -1)
    error("Too many invalid entries!")
end


% MATERIAL INPUT___________________________________________________________

material = Print_M_Menu(MAXTRIES);

if (material  == -1)
    error("Too many invalid entries!")
end


% CALLING Geometry.m and Material.m________________________________________


% COMPUTING RECCOMENDED LOADs______________________________________________


% Printing the tables______________________________________________________

% GENERATE PLOTS___________________________________________________________
% HELPER FUNCTIONS_________________________________________________________
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

    op = input('Option: ');

    if ( ~isempty(op) && isnumeric(op) && ( (op >0) && (op <=5) ) )
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

    op = input('Option: ');

    if ( ~isempty(op) && isnumeric(op) && ( (op >0) && (op <=2) ) )
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
    disp('	1 - White Oak');
    disp('	2 - Western White Pine');
    disp('	3 - Red Maple');
    disp('	4 - Particle board');
    disp('	5 - Plywood');
    disp('	6 - Aluminum');
    disp('	7 - Steel');

    op = input('Option: ');
    if isempty(op)
        op = -1;
    end
    % Recursive loop condition
    if ( ~isempty(op) && isnumeric(op) ) && ( (op >0) && (op <=7) ) 
        material = op;
    else
        material = Print_M_Menu(tries - 1);
    end

end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Analyze_Material.m
% EAS230
% Professor Sabato, Professor Ali
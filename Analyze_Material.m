% DECLARATIONS_____________________________________________________________

cross_section = 'none';
orientation = 'none';
ErrorCounter = 0;
cs_area = 0.01; % units in m^2
length = 3; % units in m
S_f = 4; % unitless
g = 9.81; % units in m/s^2
M = 101; % unitless

% CROSS SECTION INPUT______________________________________________________
Print_CS_Menu;

while(cross_section < 1 || cross_section > 5)
Print_CS_Menu;
ErrorCounter = ErrorCounter + 1;
if (ErrorCounter == 3)
    %Call Error Function
end
end

ErrorCounter = 0;

%ORIENTATION INPUT_________________________________________________________

Print_O_Menu;

while(orientation < 1 || orientation > 2)
Print_O_Menu;
ErrorCounter = ErrorCounter + 1;
if (ErrorCounter == 3)
    %Call Error Function
end
end

ErrorCounter = 0;

% CALLING Geometry.m and Material.m________________________________________

[a, b, I] = Geometry(cross_section, cs_area, orientation);
[rho_oak, E_oak, sigma_oak] = Material('oak');
[rho_pine, E_pine, sigma_pine] = Material('pine');
[rho_maple, E_maple, sigma_maple] = Material('maple');

%COMPUTATING RECCOMENDED MAX LOAD__________________________________________

deltaX = (length / M - 1);



function Print_CS_Menu()
disp('\nChoose a cross-section');
disp('    1 - Circular');
disp('    2 - Rectangular');
disp('    3 - I-Beam');
disp('    4 - T-Beam');
disp('    5 - L-Beam\n');
cross-section = input('Option: ', 's');
end

function Print_O_Menu()
disp('\nChoose an orientation');
disp('    1 - Vertical');
disp('    2 - Horizontal');
orientation = input('Option: ','s');
end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Analyze_Material.m
% EAS230
% Professor Sabato, Professor Ali
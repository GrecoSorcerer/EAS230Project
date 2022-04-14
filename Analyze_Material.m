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
F              = 0;     % units in N/Pa


% CROSS SECTION INPUT______________________________________________________

cross_section = Print_CS_Menu(MAXTRIES);

if (cross_section  == -1)
    error("Too many invalid entries!")
end



%ORIENTATION INPUT_________________________________________________________

orientation = Print_O_Menu(MAXTRIES);

if (orientation  == -1)
    error("Too many invalid entries!")
end



% CALLING Geometry.m and Material.m________________________________________
[a,         b,       I]           = Geometry(cross_section, cs_area, orientation);
%{
[rho_oak,   E_oak,   sigma_oak]   = Material(1);        % White Oak
[rho_pine,  E_pine,  sigma_pine]  = Material(2);     % Western white pine
[rho_maple, E_maple, sigma_maple] = Material(3);  % Red maple
%}

% Initialize Materials data materix
Mats = zeros(7,3);
for m = 1:7
    % Populate the Material data matrix with all of our material data.
    % It took way too long to remember I could take multiple outputs from a
    % function and place them as the values in an initialized array.
    Mats(m,:) = Material(m);
end

%COMPUTATING RECCOMENDED MAX LOAD__________________________________________
% Compute the change in x
dx  = L / (M -1);

sigmaMax = zeros(1,7);
F = zeros(1,7);


for material = 1:7

% Calculate max safe stress
sigmaMax(material) = Mats(material,3)/safety_factor;

% Calculate the load
F(material) = ( sigmaMax(material) * ( 4 * I ) ) ...
/ ( max(a,b) * (L) );

end

% Initialize point load array
f_m = zeros([7,M]);
m = 1:M; % indexing array
% Compute the point load.

for material = 1:7
    f_m(material,m == (M+1)/2) = F(material)./dx;
end

f_m = f_m';

% Compute 7 values for mu
mu = Mats(:,1).*cs_area;

% Initialize deformation data matrix
Z_max = zeros(1,7);

Z = zeros(7,M);

for material = 1:7
    % Populate deformation data matrix 
    Z(material,:) = Deformation(g,mu(material,1),Mats(material,2),I,dx,f_m(:,material));
    Z_max(material) = max( abs(Z(material,:)) );
end
%Place Z into output format by transposing it
Z = Z';

file_name = [CROSS_SECTION(cross_section) '_' ORIENTATION(orientation) '_deformation.mat'];
save(file_name,"Z","-mat");

%Printing the table________________________________________________________

fprintf('For a %s cross-section in a %s orientation\n', cross_section, orientation);
disp('          Material   Recommended max load   Failure load   Maximum deformation   Weight');
disp('                                      [N]            [N]                  [mm]     [kg]');

for material = 1:7
fprintf('%s %f %f %f %f\n', MATERIAL(material), F(material)./safety_factor, F(material), Z_max(material), mu(material)*g*L);
end

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

    op = input('Option: ');

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

    op = input('Option: ');

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
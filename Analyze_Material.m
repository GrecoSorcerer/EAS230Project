% DECLARATIONS_____________________________________________________________

MAXTRIES       = 3;
% ORIENTATION   CROSS_SECTION   MATERIAL  are hash maps. They are used to
% automatically recall the string MaterialName associated with a given
% integer 'key'. ex1 MATERIAL(7) returns a value of 'Steel'.
%
% This is used when generating tables and files names. Usually when looping
% or when related to a user related input. 
% ex2 let material = 4 from a user input, 
%         then MATERIAL(material) = 'Particle board'.
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
% Calling the xsection menu
cross_section = Print_CS_Menu(MAXTRIES);
% Error Handling
if (cross_section == -1)
    error("Too many invalid entries!")
end

%ORIENTATION INPUT_________________________________________________________
% Calling the orientation menu
orientation = Print_O_Menu(MAXTRIES);
% Error Handling
if (orientation  == -1)
    error("Too many invalid entries!")
end

% CALLING Geometry.m and Material.m________________________________________

% Get the geometry data by calling the geometry function
[a, b, I] = Geometry(cross_section, cs_area, orientation);

% Initialize Materials data materix
material_data = zeros(7,3);

% Get material data by iteratively calling Material function
for material = 1:7
    [rho, E, sigma] = Material(material);
    material_data(material,:) = [rho, E, sigma];
end

%COMPUTATING RECCOMENDED MAX LOAD__________________________________________

% Compute the change in x
dx  = L / (M -1);

% Init sigma max table
sigmaMax = zeros(1,7);
% Init load force table
F = zeros(1,7);

% Calculate max safe stress and load for each material
for material = 1:7

    % Calculate max safe stress
    sigmaMax(material) = material_data(material,3)/safety_factor;

    % Calculate the load
    F(material) = ( sigmaMax(material) * ( 4 * I ) ) ...
    / ( max(a,b) * (L) );

end

% Initialize point load array
f_m = zeros([7,M]);
m = 1:M; % indexing array

% Compute the point load.
for material = 1:7
    f_m(material,m == (M+1)/2) = F(material)/dx;
end

% flip for what function will expect
f_m = f_m';

% Compute 7 values for mu
mu = material_data(:,1).*cs_area;

% Initialize deformation data matrix
Z_max = zeros(1,7);

% Init Z
Z = zeros(7,M);

% Compute the deformation and max deformationfor each material
for material = 1:7
    % Populate deformation data matrix
    Z(material,:) = Deformation(g,mu(material,1),material_data(material,2),I,dx,f_m(:,material));
    Z_max(material) = max( abs(Z(material,:)) );
end

%Place Z into output format by transposing it
Z = Z';

% Save the data
file_name = [CROSS_SECTION(cross_section) '_' ORIENTATION(orientation) '_deformation.mat'];
save(file_name,"Z","-mat");

%Printing the table________________________________________________________

fprintf('For a %s cross-section in a %s orientation\n', CROSS_SECTION(cross_section), ORIENTATION(orientation));
disp('          Material   Recommended max load   Failure load   Maximum deformation   Weight');
disp('                                      [N]            [N]                  [mm]     [kg]');
% Print each row of the data set to the command window.
for material = 1:7
    fprintf('%18s           %12.4f   %12.4f              %8.4f  %7.2f\n', MATERIAL(material), F(material), F(material).*safety_factor, Z_max(material)*1000, mu(material)*g*L);
end

%Creating z vs. x figures__________________________________________________

% Compute the axis part.
x = ((m-1)./(M-1)).*L;

figure(2);
    % Plot the deformations wrt x
    plot(x,Z,LineWidth=2)
    % Set the fig configs
    grid on
    title('Plot of Deformations vs. X');
    ylabel('z [mm]');
    xlabel('x [m]');
    legend('White Oak','Western White Pine','Red Maple',...
    'Particle board','Plywood','Aluminum','Steel');
    axis([ min(x),        max(x),   ...
           min(min(Z))*2, abs(min(min(Z)))*4])

%creating the max load vs. Young's modulus figures_________________________

YM = material_data(:,2);

figure(3);
    % Plot the data using a loglog format
    loglog(F,YM,'d','markerfacecolor','r','markersize',10);
    % Set text to the points and the target locations
    text(F*1.05,YM*0.95,{'White Oak', 'Western White Pine', 'Red Maple', 'Particle board', 'Plywood', 'Aluminum', 'Steel'}, 'FontSize', 8);
    % Set the fig configs
    grid on
    title('Recommended Maximum Load vs. Young''s Modulus');
    xlabel('Recommended Maximum Load [N]');
    ylabel('Young''s modulus [N/m^2]');
    axis([0, max(F)*1.125, ...
          0, max(YM)*1.125]);

%HELPER FUNCTIONS__________________________________________________________
function [cross_section] = Print_CS_Menu(mtries)
    
    tries = 0;
    while tries < mtries
        cross_section = -1;
    
        disp('Choose a cross-section');
        disp('    1 - Circular');
        disp('    2 - Rectangular');
        disp('    3 - I-Beam');
        disp('    4 - T-Beam');
        disp('    5 - L-Beam');
        
        tries = tries + 1;
        op = input('Option: ','s');
        % resolve op to a double. If not a number then will return NaN
        op  = str2double(op);

        if ( ~isempty(op) && isnumeric(op) && ( (op >0) && (op <=5) ) )
            cross_section = op;
            break
        end
    end
end

function [orientation] = Print_O_Menu(mtries)
    
    tries = 0;
    while tries < mtries
        orientation = -1;
    
        disp('Choose an orientation');
        disp('    1 - Vertical');
        disp('    2 - Horizontal');
    
        tries = tries + 1;
        % take input as a string, to avoid errors when user presses a non numeric key
        op = input('Option: ','s');
        % resolve op to a double. If not a number then will return NaN
        op  = str2double(op);
        
        % return if valid
        if ( ~isempty(op) && isnumeric(op) && ( (op >0) && (op <=2) ) )
            orientation = op;
            break
        end
    end
end

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Analyze_Material.m
% EAS230
% Professor Sabato, Professor Ali

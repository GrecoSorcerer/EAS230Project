% DECLARATIONS_____________________________________________________________

MAXTRIES = 3;

% ORIENTATION   CROSS_SECTION   MATERIAL  are hash maps. They are used to
% automatically recall the string MaterialName associated with a given
% integer 'key'. ex1 MATERIAL(7) returns a value of 'Steel'.
%
% This is used when generating tables and files names. Usually when looping
% or when related to a user related input. 
% ex2 let material = 4 from a user input, 
%         then MATERIAL(material) = 'Particle board'.
ORIENTATION    = containers.Map([1,2],{'Vertical','Horizontal'});
CROSS_SECTION  = containers.Map([1,2,3,4,5], ...
                                {'Circular', 'Rectangular', ...
                                 'I-beam',   'T-beam',      ...
                                 'L-beam'});
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
%F              = 0;     % Load in N/m, no default given, assuming 1.

%CALLING MENUS_____________________________________________________________

% Try to display cross_section menu with MAXTRIES
cross_section = Print_CS_Menu(MAXTRIES);

% error handling
if (cross_section  == -1)
    error("Too many invalid entries!")
end

% Try to display orientation menu with MAXTRIES
orientation = Print_O_Menu(MAXTRIES);

% error handling
if (orientation  == -1)
    error("Too many invalid entries!")
end

% Try to display material menu with MAXTRIES
material = Print_M_Menu(MAXTRIES);

% error handling
if (material  == -1)
    error("Too many invalid entries!")
end

%THE BODY__________________________________________________________________

% Calling the Material function for rho E sigma
[rho, E, sigma] = Material(material);

% Calling the Geometry function for a b I
[a, b, I] = Geometry(cross_section, cs_area, orientation);

% Compute the change in x
dx  = L / (M - 1);

% Calculate max safe stress
sigmaMax  = sigma/safety_factor;

% Calculate the load
F = ( sigmaMax * ( 4 * I ) ) ...
/ ( max(a,b) * (L) );

% Initialize point load array
f_m = zeros([1,M]);

m = 1:M; % indexing array

% Compute the point load.
f_m(m == (M+1)/2) = (F)/dx;
% Put point load into format the Deformation function expects
f_m = f_m';

% Calculate mu
mu = rho*cs_area;

% Calculate the deformation of the beam
[z]   = Deformation(g,mu,E,I,dx,f_m);
z_max = max(abs(z));
% Define the independent variable x wrt length L  (x_m)
x = ((m-1)./(M-1)).*L;

% This may need to be change the test uses a differnt file name the the
% design document.
file_name = "Single_Beam.mat";

% Set Material XSection and Orientation to save, and use later for string
% concatenations.
Beam_Material = MATERIAL(material);
Beam_XSection = CROSS_SECTION(cross_section);
Orientation   = ORIENTATION(orientation);

save(file_name,"Beam_Material", "Beam_XSection",   ...
               "Orientation", "a","b", "I", "rho", ...
               "E", "sigma", "cs_area", "L",       ...
               "sigmaMax", "x", "z");

% fig1 figure(1) handle
fig1 = ...
figure(1);
    
    % draw the plot of deformation
    plot(x,z,'g', ...
        'LineWidth',2)
    % Set fig configs
    grid on
    title("The deformation for a beam made of " + Beam_Material + " and a" + newline ...
        + Beam_XSection + " cross-section in a " + Orientation + '.')
    % set the axis so the deformation is less exagerated.
    axis([ min(x),     max(x),   ...
           min(z)*20, max(abs(z))*7])
    xlabel("Length [m]")
    ylabel("Deformation [mm]")

% Print output table
fprintf("For a beam made of %s and a %s cross-section in a %s orientation,\n", Beam_Material, Beam_XSection, Orientation);
fprintf("Recommended max load [N]: %4.3f\n", F);
fprintf("Failure load [N]: %5.3f\n",F*safety_factor);
fprintf("Maximum deformation [mm]: %2.4f\n", z_max*1000);
fprintf("Weight [kg]: %2.1f\n",(mu*g*L));

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

function [material] = Print_M_Menu(mtries)
    
    tries = 0;
    while tries < mtries
        material = -1;
    
        disp('Choose a material');
        disp('	1 - White Oak');
        disp('	2 - Western White Pine');
        disp('	3 - Red Maple');
        disp('	4 - Particle board');
        disp('	5 - Plywood');
        disp('	6 - Aluminum');
        disp('	7 - Steel');
        
        tries = tries + 1;
        % take input as a string, to avoid errors when user presses a non numeric key
        op = input('Option: ', 's');
        % resolve op to a double. If not a number then will return NaN
        op = str2double(op);
                
        % return if valid
        if ( ~isempty(op) && isnumeric(op) ) && ( (op >0) && (op <=7) ) 
            material = op;
            break
        end
    end
end


% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Single_Case.m
% EAS230
% Professor Sabato, Professor Ali
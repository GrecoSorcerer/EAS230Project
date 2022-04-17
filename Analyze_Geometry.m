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

% MATERIAL INPUT___________________________________________________________

material = Print_M_Menu(MAXTRIES);

if (material  == -1)
    error("Too many invalid entries!")
end


% CALLING Geometry.m and Material.m________________________________________

% init geometry data table, row 1-5 are vert, 6-10 are horiz
geometry_data = zeros(10,3); 
for cross_section = 1:5

    % Grab vertical
    geometry_data(cross_section,:) = Geometry(cross_section, cs_area, 1);
    % Grab Horiz
    geometry_data(cross_section+5,:) = Geometry(cross_section, cs_area, 2);
end

% grab material data <rho E sigma>
[rho, E, sigma] = Material(material);


% COMPUTING LOAD AND DEFORMATION___________________________________________

% Compute the change in x
dx  = L / (M - 1);

% Calculate max safe stress
sigmaMax = sigma/safety_factor;

F_data = zeros(1,10);

f_m = zeros(M,10);

m = 1:M; % indexing array

% Calculate the Rec. Max Force
for beam = 1:10
    %Compute the Rec. Max Load with sigmaMax
    F_data(beam) = ( sigmaMax * ( 4 * geometry_data(beam,3) ) ) ...
        / ( max(geometry_data(beam,1:2)) * (L) );
    % Compute the point load
    f_m(m == (M+1)/2, beam) = F_data(beam)/dx;
end

% Compute mu 
mu = rho.*cs_area;

% Init deformations table, col 1-5 are vert, 6-10 are horiz.
Z = zeros(M,10);

for beam = 1:10

    Z(:,beam) = Deformation(g, mu, E, geometry_data(beam,3),dx,f_m(:,beam));
end

% Printing the tables______________________________________________________
x = ((m-1)./(M-1)).*L;
% GENERATE PLOTS___________________________________________________________
% fig1 figure(1) handle
fig1 = ...
figure(1);
    
    % draw the plot of deformation
    subplot(1,2,1);
    plot(x,Z(:,[1:5]),'g', ...
        'LineWidth',2)
    grid on

    subplot(1,2,2);
    plot(x,Z(:,[6:10]),'g', ...
        'LineWidth',2)
    grid on

    title("The deformation for a beam made of " + Beam_Material + " and a" + newline ...
        + Beam_XSection + " cross-section in a " + Orientation + '.')

    % set the axis so the deformation is less exagerated.
    axis([ min(x),     max(x),   ...
           min(Z)*20, max(abs(Z))*7])

    xlabel("Length [m]")
    ylabel("Deformation [mm]")

% HELPER FUNCTIONS_________________________________________________________
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
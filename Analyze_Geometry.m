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

% MATERIAL INPUT___________________________________________________________

% Calling helper function to print displays and get user inputs
material = Print_M_Menu(MAXTRIES);

% Error Handling
if (material  == -1)
    error("Too many invalid entries!")
end

% CALLING Geometry.m and Material.m________________________________________

% init geometry data table, row 1-5 are vert, 6-10 are horiz
geometry_data = zeros(10,3); 
% Get the geometry data for both orientations
for cross_section = 1:5
    
    [a1,b1,I1] = Geometry(cross_section, cs_area, 1);
    [a2,b2,I2] = Geometry(cross_section, cs_area, 2);
    % Grab vertical
    geometry_data(cross_section,:) = [a1,b1,I1];
    % Grab Horiz
    geometry_data(cross_section+5,:) = [a2,b2,I2];
end

% grab material data <rho E sigma>
[rho, E, sigma] = Material(material);

% Compute the change in x
dx  = L / (M - 1);

% Calculate max safe stress
sigmaMax = sigma/safety_factor;

% init the force data variables
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
Z_max = zeros(1,10);

% Calculate Deformation vectors for dataset
for beam = 1:10
    Z(:,beam) = Deformation(g, mu, E, geometry_data(beam,3),dx,f_m(:,beam));
    Z_max(beam) = max( abs( Z(:,beam) ) );
end

%saving deformations matrix__________________________________________________

file_name = [MATERIAL(material) '_deformation.mat'];
save(file_name,"Z","-mat");

% Printing the table______________________________________________________

fprintf('For a beam made of %s with a weight of %f kg,\n',MATERIAL(material),mu*g*L);

disp('  Vertical Geometry   Recommended Max Load   Failure Load   Maximum Deformation');
disp('                                       [N]            [N]                  [mm]');
% Printing table rows to command window
for cs = 1:5
    fprintf('%19s %22.3f %14.4f %21.4f\n', CROSS_SECTION(cs),F_data(cs),...
    F_data(cs).*safety_factor,Z_max(cs)*1000);
end

fprintf('\n');
% Printing table rows to command window
disp('Horizontal Geometry   Recommended Max Load   Failure Load   Maximum Deformation');
disp('                                       [N]            [N]                  [mm]');
for cs = 1:5
    fprintf('%19s %22.3f %14.4f %21.4f\n', CROSS_SECTION(cs),F_data(cs+5),...
    F_data(cs+5).*safety_factor,Z_max(cs+5)*1000);  
end

fprintf('\n');

% GENERATE PLOTS___________________________________________________________

% set x axis indexing array
x = ((m-1)./(M-1)).*L;

fig4 = ...
figure(4);
    
    % draw the plot of deformations________________________________________
    subplot(1,2,1);

        % Plot the deformation data wrt x
        plot(x,Z(:,1:5), ...
            'LineWidth',1)

        % Set the plot configurations
        grid on
        title("Deformation of Beams Made of " +...
               MATERIAL(material) + newline + "Oriented " + ORIENTATION(1) + "ly.");
        xlabel('Length of Beam [m]');
        ylabel('Deformation [mm]');
        legend('Circular','Rectangular','I-Beam',...
        'T-Beam','L-Beam');
        axis([ min(x),        max(x),   ...
               min(min(Z))*2, abs(min(min(Z)))*2])

    subplot(1,2,2);

        % Plot the deformation data wrt x
        plot(x,Z(:,6:10), ...
            'LineWidth',1)
    
        % Set the plot configurations
        grid on
        title("Deformation of Beams Made of " +...
               MATERIAL(material) + newline + "Oriented " + ORIENTATION(2) + "ly.");
        xlabel('Length of Beam [m]');
        ylabel('Deformation [mm]');
        legend('Circular','Rectangular','I-Beam',...
        'T-Beam','L-Beam');
        axis([ min(x),        max(x),   ...
               min(min(Z))*2, abs(min(min(Z)))*2])
    
% Set an initial size and position for this figure
fig4.Position = [100 100 1000 500];

% draw the plot of rec. max load vs. second moment_____________________
fig5 = ...
figure(5);

    % Plot the data in loglog format
    loglog(F_data',geometry_data(:,3),'d','markerfacecolor','c','markersize',10);


    % Set the text at the target positions for all non overlapping texts
    text(F_data([2,4,6:1:10])',geometry_data([2,4,6:1:10],3)*0.95, ...
        { ...
            'Vertical Rectangular', ...
            'Vertical T-Beam', 'Horizontal Circular',...
            'Horizontal Rectangular','Horizontal I-Beam', 'Horizontal T-Beam', ...
            'Horizontal L-Beam' ...
        }, 'FontSize', 8);
    
    % Set the text at the target positions for text that would have
    % overlapped
    text(F_data([1,3,5])',geometry_data([1,3,5],3)*0.87, ...
        { ...
            'Vertical Circular', ...
            'Vertical I-Beam', ...
            'Vertical L-Beam'
        }, 'FontSize', 8);

    % Set the plot configurations
    title('Recommended Maximum Load vs. Second Moment of Area');
    xlabel('Recommended Maximum Load [N]');
    ylabel('Second Moment of Area [mm^4]');
    grid on
    axis([ min(F_data)*0.65, max(F_data)*1.35, min(min(geometry_data))*0.65, abs( max( geometry_data(:,3) ) )*1.35 ])
    
% Set an initial size and position for this figure
fig5.Position = [1200 100 550 550];

% HELPER FUNCTIONS_________________________________________________________
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
% Analyze_Geometry.m
% EAS230
% Professor Sabato, Professor Ali

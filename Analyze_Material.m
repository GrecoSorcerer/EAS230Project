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



%ORIENTATION INPUT_________________________________________________________

orientation = Print_O_Menu(MAXTRIES);

if (orientation  == -1)
    error("Too many invalid entries!")
end



% CALLING Geometry.m and Material.m________________________________________

[a, b, I] = Geometry(cross_section, cs_area, orientation);

% Initialize Materials data materix
Mats = zeros(7,3);

for material = 1:7
    [rho, E, sigma] = Material(material);
    Mats(material,:) = [rho, E, sigma];
end

%COMPUTATING RECCOMENDED MAX LOAD__________________________________________

% Compute the change in x
dx  = L / (M -1);

% Init sigma max table
sigmaMax = zeros(1,7);
% Init load force table
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
    f_m(material,m == (M+1)/2) = F(material)/dx;
end

% flip for what function will expect
f_m = f_m';

% Compute 7 values for mu
mu = Mats(:,1).*cs_area;

% Initialize deformation data matrix
Z_max = zeros(1,7);

% Init Z
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

fprintf('For a %s cross-section in a %s orientation\n', CROSS_SECTION(cross_section), ORIENTATION(orientation));
disp('          Material   Recommended max load   Failure load   Maximum deformation   Weight');
disp('                                      [N]            [N]                  [mm]     [kg]');

for material = 1:7
    fprintf('%18s           %12.4f     %10.4f              %8.4f  %7.2f\n', MATERIAL(material), F(material), F(material).*safety_factor, Z_max(material)*1000, mu(material)*g*L);
end

%Creating z vs. x figures__________________________________________________

%Figure(2);
%plot(z(1),MATERIAL(1),'r',z(2),MATERIAL(2),'g',z(3),MATERIAL(3),'b',z(4),MATERIAL(4),'c',z(5),MATERIAL(5),'m',z(6),MATERIAL(6),'y',z(7),MATERIAL(7),'k')
%grid on
%title('Plot of Deformations vs. X');
%ylabel('x (material)');
%xlabel('z (deformations)');

%creating the max load vs. Young's modulus figures_________________________

YM = Mats(:,2);

loglog(F,YM,'d','markerfacecolor','r','markersize',10);
text(F*0.90,YM*0.85,{'White Oak', 'Western White Pine', 'Red Maple', 'Particle board', 'Plywood', 'Aluminum', 'Steel'}, 'FontSize', 8);
grid on
title('Recommended Maximum Load vs. Young''s Modulus');
xlabel('Recommended Maximum Load');
ylabel('Young''s modulus');
axis([0 10^4 0 10^11.5]);

%HELPER FUNCTIONS__________________________________________________________
function [cross_section] = Print_CS_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response.
    % Returns -1 if tries exceed maximum allowance
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

% Salvatore L Greco <slgreco@buffalo.edu>
% Alexander M Gross <amgross@buffalo.edu>
% Analyze_Material.m
% EAS230
% Professor Sabato, Professor Ali
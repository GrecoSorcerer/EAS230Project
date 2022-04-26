MAXTRIES = 3; % Max tries before throwing error for failed user inputs

% Display the Beams Menu with MAXTRIES failed inputs allowed.
% I use reccursion because its more interesting
option = Print_Beams_Menu(MAXTRIES);
% !! Handling error in switch case
% Run the selected beam simulation script
switch option
    case 1
        Analyze_Geometry;
    case 2
        Analyze_Material;
    case 3
        Single_Case;
    case -1
    % The error case
        error("Too many invalid entries!")
    otherwise
        error("Unexpected error occured!")
end

% Helper function for printing the menu
function [option] = Print_Beams_Menu(mtries)
    % Set option to negative one, if function call fails, this will be the 
    % return value.
    option = -1;
    
    tries = 0;
    while tries < mtries
        disp('Choose one');
        disp('     1 - analyze geometry');
        disp('     2 - analyze material');
        disp('     3 - analyze single case');
    
        tries = tries + 1;
        op = input('Option: ','s');
        op = str2double(op);
    
        % recursive loop condition (if valid set return val, and exit, other wise try again)
        if ( ~isempty(op) && isnumeric(op)) && ( (op >0) && (op <=3) )
            option = op;
            break
        end
    end
end
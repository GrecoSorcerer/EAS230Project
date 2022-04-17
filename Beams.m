MAXTRIES = 3; % Max tries before throwing error for failed user inputs

% Display the Beams Menu with MAXTRIES failed inputs allowed.
% I use reccursion because its more interesting
option = Print_Beams_Menu(MAXTRIES);

% error handling
if (option  == -1)
    error("Too many invalid entries!")
end

% Run the selected beam simulation script
switch option
    case 1
        Analyze_Geometry;
    case 2
        Analyze_Material;
    case 3
        Single_Case;
end


function [option] = Print_Beams_Menu(tries)
    % Recursive function. Calls itself up to tries times, to get a valid
    % response.
    % Returns -1 if exceeds tries
    % Be sure to catch the -1 as an error
    
    % Recursive Exit Condition
    if (tries == 0)
        option = -1;
        return;
    end

    disp('Choose one');
    disp('     1 - analyze geometry');
    disp('     2 - analyze material');
    disp('     3 - analyze single case');

    op = input('Option: ');

    % recursive loop condition (if valid set return val, and exit, other wise try again)
    if ( ~isempty(op) && isnumeric(op)) && ( (op >0) && (op <=3) )
        % If we reach here the input was valid. Because of how matlab
        % functions work option is already set up to be returned when the
        % function closes.
        %
        % Entering this part of the condtional allows the function to exit
        % the recursive loop
        option = op;
    else
        % Recieved unexpected input, call menu func again and approach exit-condition.
        option = Print_Beams_Menu(tries - 1); 
    end

end
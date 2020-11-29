% Input data handling

function [vX,vY,guess,p] = input_data(cInput)
    if length(cInput) < 5
        [infile, row_or_col,index_1,index_2,guess,delim,p] = get_input();
    else
        infile = cInput{1};
        row_or_col = cInput{2};
        index_1 = cInput{3};
        index_2 = cInput{4};
        p = cInput{5};
        guess = [5.640,9.4*280*10^(-5),280];
        delim = ' ';
    end
    [vX,vY] = import_file(infile, row_or_col, index_1, index_2, delim);
end

function [infile, row_or_col,index_1,index_2,guess,delim,p] = get_input()
    disp('Select data file.');
    infile = uigetfile;
    
    prompt1='Enter 0 if the data is organized in rows, and 1 if otherwise.\n';
    row_or_col = input(prompt1);
    
    prompt2='Enter the index of the independent variable to be analyzed.\n';
    index_1 = input(prompt2);
    
    prompt3='Enter the index of the dependent variable to be analyzed.\n';
    index_2 = input(prompt3);

    prompt4='Enter the degree of the highest polynomial to be fitted.\n';
    p = input(prompt4);
    
    %prompt5='Enter an initial value for nonlinear optimzation.';
    guess = [5.640,9.4*280*10^(-5),280];
    
    delim = ' '; % placeholder for user specified delimiter functionality
end

% This function imports and formats the data given by user input from
% get_input.
function [vX,vY] = import_file(infile, row_or_col, index_1, index_2, delim)
    Data = importdata(infile,delim);
    % the rest of the code assumes that the data is formatted with vX as a
    % horizontal vector and vY as a vertical vector. This should perhaps be
    % changed, but isn't really worth the effort.
    if row_or_col == 1
        vX = transpose(Data(:,index_1));
        vY = Data(:,index_2);
    else
        vX = Data(index_1,:);
        vY = transpose(Data(index_2,:));
    end
end
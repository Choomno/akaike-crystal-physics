% This is a script to calculate the ordinary least-squares regression on a
% set of bivariate data, as well as provide the accompanying standard error
% values, adjusted R^2, AICc values, and evidence ratios.
% Written by Choomno Moos, affiliated Portland State University.

% Last changed 4/8/2020.

format long
warning off
clc
clear

% object definitions
% Note that usually each variable definition is prefixed with a character
% corresponding to the data type of the variable. For example, vX is a
% vector, mX is a matrix, and cX is a cell.
vX = [];               % vector of in dependent variable elements

vY = [];               % vector of dependent variable elements

cParameters = {};      % cell to contain model parameter vectors

cStderr = {};          % cell to contain standard error vectors

vAICc = [];            % vector to contain AICc values for a set of models

vAkaike_weights = [];  % vector to contain Akaike weights

cEvidence_ratios = {}; % cell to contain evidence ratios

cModels = [];          % cell to contains information about models

vRsq_adj = [];         % vector to contain adjusted R^2

% input parameters
infile = 'swaminathan.txt';
                       % Name of file containing input data. Note that only
                       % bivariate data will work input tab-delimited .txt 
                       % file with data pairs as columns. See the
                       % accompanying example data.txt file for an example.

p = 3;                % maximum degree of polynomial to be considered


row_or_col = 0;        % 1 indicates data is formatted in columns,
                       % other indicates data is formatted in rows.

index_1 = 1;           % specifies the index (row or column number) of the
                       % independent variable vector from .txt file.

index_2 = 2;           % specifies the index (row or column number) of the
                       % dependent variable vector.
                       
delimiter = ' ';       % specifies delimiter used in .txt file to be read.

figure_xlabel = 'Temperature $(K)$'; %'Pressure $(GPa)$' 

figure_ylabel = 'Parameter $(\AA)$';

figure_title = 'Swaminathan $a$';

% imports data from file. See function for details
[vX, vY] = input_data(infile, row_or_col, index_1, index_2, delimiter);

% control loop performs least-squares regressions for each 1:p degrees
for i = 1:p
    mX = design_matrix(i,vX);
    [mX,vY,vBeta,vStderr,vRsq_adj(i),s_sq] = least_squares(mX,vY,i);
    cParameters{i} = vBeta;
    cStderr{i} = vStderr;
end

cModels = [cParameters; cStderr];

% control loop finds AICc values for each model
for i = 1:p
    mX = design_matrix(i,vX);
    vBeta = cModels{1,i};
    [vAICc(i),sigmahat_sq] = AICc_calc(mX,vY,vBeta);
end

% calculation of Akaike weights
vAkaike_weights = Akaike_weights(vAICc,p); 

% calculation of evidence ratios
cEvidence_ratios = evidence_ratios(vAkaike_weights,p);

% outputs results into a text file
data_table = output_data(cModels, vRsq_adj, vAICc, ...
                          vAkaike_weights, cEvidence_ratios, p);

% the rest produces graphics
syms x;
y1 = cParameters{1,1}(1,1) + cParameters{1,1}(2,1)*x;
y2 = cParameters{1,2}(1,1) + cParameters{1,2}(2,1)*x ...
     + cParameters{1,2}(3,1)*x^2;
y3 = cParameters{1,3}(1,1) + cParameters{1,3}(2,1)*x ...
     + cParameters{1,3}(3,1)*x^2 + cParameters{1,3}(4,1)*x^3;

figure1 = figure('Visible', 'off');
scatter(vX,exp(vY));
hold on
fplot(x,exp(y1),[(vX(1)-10) (vX(length(vX))+10)],'Linewidth',1);
hold on
fplot(x,exp(y2),[(vX(1)-10) (vX(length(vX))+10)],'Linewidth',1);
hold on
fplot(x,exp(y3),[(vX(1)-10) (vX(length(vX))+10)],'Linewidth',1);
xlabel(figure_xlabel,'Interpreter','latex');
ylabel(figure_ylabel,'Interpreter','latex');
title(figure_title,'Interpreter','latex');
lgd = legend('Experimental Data','First-Order Model', ...
             'Second-Order Model','Third-Order Model','Interpreter', ...
             'latex','Location','best');
title(lgd,'Legend','Interpreter','latex');
grid(gca,'minor');
grid on

saveas(gcf,'figure.pdf')

disp('done')

%%%%%% Function definitions %%%%%%%

% Edit this function to specify the way data is read from .txt file. If
% row_or_col is true, data is assumed to be formatted in columns, otherwise
% data is assumed to be formatted in rows. Also requires the indices for
% the independent and dependent variables to be specified.
function [vX,vY] = input_data(infile, row_or_col, index_1, index_2, delim)
    Data = importdata(infile,delim);
    if row_or_col == 1
        vX = transpose(Data(:,index_1));
        vY = Data(:,index_2);
    else
        vX = Data(index_1,:);
        vY = transpose(Data(index_2,:));
    end
end

% Constructs design matrix for regression. Takes vector and number of
% degrees of polynomials to be considered and returns a matrix.
function [mX] = design_matrix(i,vX)
    mX = [ones(length(vX),1)';vX];
    if i > 1
        for j = 1:(i-1)
            mX = [mX;vX.^(j+1)];
        end
    end
    mX = mX';
end

% Performs least squares regression. Takes the design matrix, the dependent
% data vector, and the integer 'i' indicating the degree of the polynomial
% and echoes the first two, as well as returns the parameters and standard
% error as vectors, and the adjusted R^2 and estimated variance as floats.
function [mX,vY,vBeta,vStderr,Rsq_adj,s_sq] = least_squares(mX,vY,i)
    p = i;
    vBeta = (mX'*mX)^-1 * mX' * vY;
    [i,j] = size(mX);
    vRes = vY - mX * vBeta;
    SS_res = vRes' * vRes;
    s_sq = SS_res/(i - j);
    mVarcovar = s_sq * (mX'*mX)^-1;
    % note s^2 is the OLS estimate for variance, sigmahat^2 is the 
    % MLE estimate. sigmahat^2 is biased and s^2 is unbiased.
    vStderr = sqrt(diag(mVarcovar));
    n = length(vY); % same as i
    ybar = sum(vY)/n;
    SS_tot = sum((vY - ybar).^2);
    Rsq = 1 - SS_res/SS_tot;
    Rsq_adj = 1 - ((1 - Rsq) * (n - 1))/(n - p - 1);
end

% Calculates the AICc, as well as the MLE estimate of the variance.
function [AICc,sigmahat_sq] = AICc_calc(mX,vY,vBeta)
    n = length(vY);
    [i,~] = size(mX);
    vRes = vY - mX * vBeta;
    SS_res = vRes' * vRes;
    sigmahat_sq = SS_res/i;
    K = length(vBeta) + 1;
    AICc = n * log(sigmahat_sq) + 2*K*(n/(n-K-1));
end

% Calculates Akaike weights. Takes a vector of AICc values and returns them
% as a vector.
function [vAkaike_weights] = Akaike_weights(vAICc,p)
    vDelta = [];
    vAkaike_weights = [];
    for i = 1:p
        vDelta(i) = vAICc(i) - min(vAICc);
    end
    for i = 1:p
        vAkaike_weights(i) = exp(-0.5 * vDelta(i))/sum(exp(-0.5 * vDelta));
    end
end

% Calculates evidence ratios from Akaike weights. Takes a vector of Akaike
% weights and returns a vector of evidence ratios.
function [cEvidence_ratios] = evidence_ratios(vAkaike_weights,p)
    cEvidence_ratios = {};
    vEvidence_ratios = [];
    for i = 1:p
        for j = 1:p
            vEvidence_ratios(j) = vAkaike_weights(i)/vAkaike_weights(j);
        end
        cEvidence_ratios{i} = vEvidence_ratios;
    end
end

% Outputs regressions and AIC information as multiple .xls excel files.
function [data_table] = output_data(cModels, vRsq_adj, vAICc, ...
                     vAkaike_weights, cEvidence_ratios, p)
    cData_tables = {};
    cRegression = {};
    cAIC = {};
    variable_names = {'Parameters','Standard_error', 'adjusted_R^2'};
    for i = 1:p
        cRegression{1,1} = cModels{1,i};
        cRegression{1,2} = cModels{2,i};
        cRegression{1,3} = vRsq_adj(i);
        cRegression = pad_cell(cRegression); % vectors need to be same size
        data_table = table(cRegression{1,1},cRegression{1,2}, ...
                           cRegression{1,3},'Variablenames', ...
                           variable_names);
        cData_tables{1,i} = data_table;
    end
    cAIC{1,1} = transpose(vAICc); % vectors are made to be columns
    cAIC{1,2} = transpose(vAkaike_weights);
    cAIC{1,3} = cEvidence_ratios{1,1};
    cAIC = pad_cell(cAIC);
    variable_names = {'AICc','Akaike_weights','Evidence_ratios'};
    data_table = table(cAIC{1,1},cAIC{1,2}, cAIC{1,3},'Variablenames', ...
                       variable_names);
    cData_tables{1,p+1} = data_table;
    for i = 1:p
        name = strcat('degree_',string(i), '.xls');
        writetable(cData_tables{1,i}, name);
    end
    writetable(data_table, 'AIC.xls');
end

% the rest are utility functions
function [M] = max_length(cVectors)
    A = [];
    S = size(cVectors);
    for i = 1:S(2)
        A(i) = length(cVectors{1,i});
    end
    M = max(A);
end

function [vOut] = pad_vector(vIn, size)
    vOut = NaN(size,1);
    for i=1:length(vIn)
        vOut(i) = vIn(i);
    end
end

function [cOut] = pad_cell(cIn)
    M = max_length(cIn);
    S = size(cIn);
    for j = 1:S(2)
        cOut{1,j} = pad_vector(cIn{1,j},M);
    end
end
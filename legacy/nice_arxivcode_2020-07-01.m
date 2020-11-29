% This is a script to calculate the ordinary least-squares regression on a
% set of bivariate data, as well as provide the accompanying standard error
% values, adjusted R^2, AICc values, and evidence ratios.

% Choomno Moos
% Last changed 7/1/2020.

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

cModelslin = [];       % cell to contains information about models

vRsq_adj = [];         % vector to contain adjusted R^2

% input parameters
infile = 'swaminathan.txt';
                       % Name of file containing input data. Input .txt 
                       % file with data pairs as columns or rows.

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

% Arbitrary initial value guesses for nonlinear optimization

% It would be prudent to input values based on past experiements as to
% avoid misleading local minima.

% A note for improvement: For the fitting of the rational approximation of
% the Debye model function, multivariate least squares may yield a more
% informed initial value.
vGuess = [5.640,9.4*280*10^(-5),280];

% imports data from file. See function for details
[vX, vY] = input_data(infile, row_or_col, index_1, index_2, delimiter);

% performs least-squares regressions for each 1:p degrees
for i = 1:p
    mX = design_matrix(i,vX);
    [mX,vY,vBeta,vStderr,vRsq_adj(i),s_sq] = least_squares(mX,vY,i);
    cParameters{i} = vBeta;
    cStderr{i} = vStderr;
end

cModelslin = [cParameters; cStderr];

% finds AICc values for each linear model
for i = 1:p
    mX = design_matrix(i,vX);
    vBeta = cModelslin{1,i};
    [vAICc(i),sigmahat_sq] = AICc_calc(mX,vY,vBeta);
end

% nonlinear model fits and output
[vEinsteinCoefficients, EinsteinAIC] = EinsteinFit(vX,vY,vGuess);
[vDebyeCoefficients, DebyeAIC] = DebyeFit(vX,vY,vGuess);
cModelsnonlin = {vEinsteinCoefficients, vDebyeCoefficients};

vAICc = [vAICc, EinsteinAIC, DebyeAIC];

% calculation of Akaike weights
vAkaike_weights = Akaike_weights(vAICc,p+2);

% calculation of evidence ratios
cEvidence_ratios = evidence_ratios(vAkaike_weights,p+2);

% outputs results into a text file
data_table = output_data(cModelslin, cModelsnonlin, vRsq_adj, vAICc, ...
                          vAkaike_weights, cEvidence_ratios, p+2);

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
% and echoes the first two, and returns the parameters and standard
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
function [data_table] = output_data(cModelslin,cModelsnonlin, vRsq_adj, vAICc, ...
                     vAkaike_weights, cEvidence_ratios, p)
    cData_tables = {};
    cRegression = {};
    cAIC = {};
    variable_names = {'Parameters','Standard_error', 'adjusted_R^2'};
    for i = 1:p
        cRegression{1,1} = cModelslin{1,i};
        cRegression{1,2} = cModelslin{2,i};
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
    writetable(data_table, 'polynomial_models.xls');
end

% nonlinear least squares for Einstein model
function [x_1,AICc] = EinsteinFit(vX,vY,vGuess)
    vY = transpose(vY);

    fun1 = @(x) x(1) + x(2) * ((exp(x(3) * vX.^(-1)) - 1).^(-1)) - vY;
    %fun2 = @(x) x(1) + x(2) * ((exp(x(3) * vX.^(-1)) - 1).^(-1))+ x(4) ...
    %    * ((exp(x(5) * vX.^(-1)) - 1).^(-1))  - vY;
    %fun3 = @(x) x(1) + x(2) * ((exp(x(3) * vX.^(-1)) - 1).^(-1)) + x(4) ...
    %    * ((exp(x(5) * vX.^(-1)) - 1).^(-1)) + x(6) ...
    %    * ((exp(x(7) * vX.^(-1)) - 1).^(-1))  - vY;


    options = optimoptions("lsqnonlin","Algorithm",'levenberg-marquardt',...
               "MaxFunctionEvaluations",10^(100), ...
               "MaxIterations",10^(4), ...
               "Display","iter-detailed");
               %"OptimalityTolerance",10^(-10), ...
               %"FunctionTolerance",10^(-10), ...
               %"StepTolerance",10^(-10)
               
    x0 = vGuess;
    [x_1,resnorm] = lsqnonlin(fun1,x0,[],[],options);
    %x0 = [5.640,9.4*280*10^(-5),280,1,1];
    %[x_2,resnorm2] = lsqnonlin(fun2,x0,[],[],options);
    %x0 = [5.640,9.4*280*10^(-5),280,1,1,1,1];
    %x_3 = lsqnonlin(fun3,x0,[],[],options);
    %[x,resnorm,residual,exitflag,output,lambda,jacobian] 
    syms t;
    figure(1)
    scatter(vX,vY);
    hold on
    fplot((x_1(1) + x_1(2) / (exp(x_1(3) / t) - 1)),[1, 180]);
    hold off
    %figure(2)
    %scatter(vX,vY);
    %hold on
    %fplot((x_2(1) + x_2(2) / (exp(x_2(3) / t) - 1) ...
    %    + x_2(4) / (exp(x_2(5) / t) - 1)),[1, 180]);
    %hold off
    %figure(3)
    %scatter(vX,vY);
    %hold on
    %fplot(x_3(1) + x_3(2) / (exp(x_3(3) / t) - 1) + x_3(4) ...
    %    / (exp(x_3(5) / t) - 1) + x_3(6) ...
    %    / (exp(x_3(7) / t) - 1),[1, 180]);
    %hold off

    n = length(vX);
    sigmahat_sq = resnorm/n;
    K = length(x_1) + 1;
    AICc = n * log(sigmahat_sq) + 2*K*(n/(n-K-1));
end

% nonlinear least squares for rationally approximated Debye model
% The rational approximation coefficients are obtained from a paper by
% Henry Thacher titled "Rational Approximations for the Debye Functions".
% Note that the "epsilon" error term for the approximation is neglected
% here.;
function [x_1,AICc] = DebyeFit(vX,vY,vGuess)
    vY = transpose(vY);
    n = size(vX);
    
    fun1 = @(x) x(1) + (x(2)*vX)...
    .*(3953.632*ones(n)-800.6087*(x(3)*vX.^(-1))...
    +85.07724*(x(3)*vX.^(-1)).^(2)-4.43282*(x(3)*vX.^(-1)).^(3)...
    +0.0946173*(x(3)*vX.^(-1)).^(4))...
    .*(3953.632*ones(n)+682.0012*(x(3)*vX.^(-1))...
    +143.155337*(x(3)*vX.^(-1)).^(2)+15.121491*(x(3)*vX.^(-1)).^(3)...
    +(x(3)*vX.^(-1)).^(4)).^(-1);
    
    options = optimoptions("lsqnonlin","Algorithm",'levenberg-marquardt',...
               "MaxFunctionEvaluations",10^(100), ...
               "MaxIterations",10^(4), ...
               "Display","iter-detailed");
               %"OptimalityTolerance",10^(-10), ...
               %"FunctionTolerance",10^(-10), ...
               %"StepTolerance",10^(-10)
    
    x0 = vGuess;
    [x_1,resnorm] = lsqnonlin(fun1,x0,[],[],options);
    %[x,resnorm,residual,exitflag,output,lambda,jacobian]
    syms t;
    figure(1)
    scatter(vX,vY);
    hold on
    fplot(x_1(1) + (x_1(2)/t)*(3953.632-800.6087*(x_1(3)/t)...
        +85.07724*(x_1(3)/t)^(2)-4.43282*(x_1(3)/t)^(3)...
        +0.0946173*(x_1(3)/t)^(4))/(3953.632+682.0012*(x_1(3)/t)...
        +143.155337*(x_1(3)/t)^(2)+15.121491*(x_1(3)/t)^(3)...
        +(x_1(3)/t)^(4)),[1, 180]);
    hold off
    
    n = length(vX);
    sigmahat_sq = resnorm/n;
    K = length(x_1) + 1;
    AICc = n * log(sigmahat_sq) + 2*K*(n/(n-K-1));
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
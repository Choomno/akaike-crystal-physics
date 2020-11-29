% This is a script to calculate the ordinary least-squares regression on a
% set of bivariate data, as well as a number of statistics.


% I'd like to note here that I am currently unaware of the notational
% conventions for use of 's' versus 'sigma_hat' (greek lower case letter
% with a caret above). I know that, explicitly, 's' refers to sample
% variance, and 'sigma_hat' to an estimate of the variance, however both
% symbols appear to represent the same thing when the theory is played out.

% Performs least-squares regressions for each 1:p degrees and assigns the
% results to an output cell.
function [cModelsLin,cX] = ordinary_least_squares(vX,vY,p)
    for i = 1:p
        mX = construct_design_matrix(i,vX);
        
        % Cell to hold design matrices, to be input into AICc_calc in
        % workflow.
        cX{i} = mX;
        
        [vBeta,std_dev,vStderr,R_sq,R_sq_adj,Var_OLS] = ...
            least_squares_calc(mX,vY,i);
        
        % storing all regression information
        cParameters{i} = vBeta;
        vStd_dev(i) = std_dev;
        cStderr{i} = vStderr;
        vR_sq(i) = R_sq;
        vR_sq_adj(i) = R_sq_adj;
        vVar_OLS(i) = Var_OLS;
    end
    % output cell object
    cModelsLin = ...
        {cParameters; cStderr; vStd_dev; vR_sq; vR_sq_adj; vVar_OLS};
end

% Constructs design matrix for regression. The result is of the form
% / 1  x_1 x_1^2 x_1^3 ... \
% | 1  x_2 x_2^2 x_2^3 ... |
% | 1  x_3 x_3^2 x_3^3 ... |
% | .   .   .     .    .   |
% | .   .   .     .     .  |
% \ .   .   .     .      . /
% where i is the number of columns - 1.
function [mX] = construct_design_matrix(i,vX)
    mX = [ones(length(vX),1)';vX];
    if i > 1
        for j = 1:(i-1)
            mX = [mX;vX.^(j+1)];
        end
    end
    mX = mX';
end

% Performs a least squares regression, as well as standard deviation, 
% standard error, R-squared, and adjusted R-squared. Takes the design
% matrix, the dependent data vector, and the integer 'p' indicating the
% degree of the polynomial to be fit. The function echoes the first two 
% inputs and returns the parameters and standard error as vectors, and the
% adjusted R^2 and estimated variance as floats. The reasons behind the
% computations conducted in this function are entirely mathematical, and 
% are very well known.
function [vBeta,std_dev,vStderr,R_sq,R_sq_adj,Var_OLS] = ...
                                                least_squares_calc(mX,vY,p)
    % fit parameters
    vBeta = (mX'*mX)^-1 * mX' * vY;
    
    % design matrix dimensions
    [i,j] = size(mX);
    
    % residuals
    vRes = vY - mX * vBeta;
    
    % residual sum of squares
    RSS = vRes' * vRes;
    
    % sample variance. Note that this is the OLS estimate for variance,
    % s_sq_MLE, seen in the AIC calculation, is the MLE estimate. s_sq_MLE
    % is biased and Var_OLS is unbiased.
    Var_OLS = RSS/(i - j);
    
    % variance-covariance matrix
    mVarcovar = Var_OLS * (mX'*mX)^-1;
    
    % standard error
    vStderr = sqrt(diag(mVarcovar));
    
    n = length(vY); % same as p, code included for readability
    
    % the mean of y
    y_bar = sum(vY)/n;
    
    % total sum of squares
    TSS = sum((vY - y_bar).^2);
    
    % coefficient of determination (R-squared)
    R_sq = 1 - RSS/TSS;
    
    % adjusted coefficient of determination (adjusted R-squared)
    R_sq_adj = 1 - ((1 - R_sq) * (n - 1))/(n - p - 1);
    
    % standard deviation
    % NEEDS TO BE CALCULATED
    std_dev = 0;
end


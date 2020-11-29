% Performs non-linear least squares for Einstein and Debye functions.

% nonlinear least squares for Einstein model
function [cModelsNonlin] = nonlinear_least_squares(vX,vY,guess)
    [Einstein_coefficients,AICc_Einstein] = EinsteinFit(vX,vY,guess);
    [Debye_coefficients,AICc_Debye] = DebyeFit(vX,vY,guess);
    cModelsNonlin = {Einstein_coefficients,AICc_Einstein,...
                     Debye_coefficients,AICc_Debye};
end


% For both of the following funcitons, I'm not knowledgable enough in
% numerics to know what are adequate tolerances for our fitting, here I
% have implemented more or less arbitrary optimization options.
% I have implemented the Levenberg-Marquardt fitting algorithm, which is
% known to be pretty robust. In the end, as long as the graph looks good,
% it can be assumed that global optimization is achieved.
function [x_1,AICc] = EinsteinFit(vX,vY,guess)
    vY = transpose(vY);

    % Residual Einstein function to be fitted. Note that the temperature is
    % replaced by the vector vX, and the real values are given by vY.
    fun1 = @(x) x(1) + x(2) * ((exp(x(3) * vX.^(-1)) - 1).^(-1)) - vY;

    options = optimoptions("lsqnonlin","Algorithm",'levenberg-marquardt');
               %"MaxFunctionEvaluations",10^(100), ...
               %"MaxIterations",10^(4));
               %"Display","iter-detailed", ...
               %"OptimalityTolerance",10^(-10), ...
               %"FunctionTolerance",10^(-10), ...
               %"StepTolerance",10^(-10), ...
               
    [x_1,resnorm] = lsqnonlin(fun1,guess,[],[],options);
    syms t;
    figure(1)
    scatter(vX,vY);
    hold on
    fplot((x_1(1) + x_1(2) / (exp(x_1(3) / t) - 1)),...
        [vX(1)-10, vX(length(vX))+10]);
    hold off

    % AICc calculation
    n = length(vX);
    sigmahat_sq = resnorm/n;
    K = length(x_1) + 1;
    AICc = n * log(sigmahat_sq) + 2*K*(n/(n-K-1));
end

% nonlinear least squares for rationally approximated Debye model

function [x_1,AICc] = DebyeFit(vX,vY,guess)
    vY = transpose(vY);
    n = size(vX);
    
    % function to be fitted The rational approximation coefficients are 
    % obtained from a paper by Henry Thacher titled "Rational 
    % Approximations for the Debye Functions". Note that the "epsilon"
    % error term for the approximation is neglected here, and it is small
    % enough to justifiably do so.
    
    % First the approximation coefficients
    c_1 = 3953.632;
    c_2 = -800.6087;
    c_3 = 85.07724;
    c_4 = -4.43282;
    c_5 = 0.0946173;
    c_6 = 682.0012;
    c_7 = 143.155337;
    c_8 = 15.121491;
    
    % The following is the Debye function approximation to be fitted.
    fun1 = @(x) x(1) + (x(2)*vX)...
    .*(((c_1*ones(n))+(c_2*(x(3)*vX.^(-1)))...
    +(c_3*(x(3)*vX.^(-1)).^(2))+(c_4*(x(3)*vX.^(-1)).^(3))...
    +(c_5*(x(3)*vX.^(-1)).^(4)))...
    .*((c_1*ones(n))+(c_6*(x(3)*vX.^(-1)))...
    +(c_7*(x(3)*vX.^(-1)).^(2))+(c_8*(x(3)*vX.^(-1)).^(3))...
    +((x(3)*vX.^(-1)).^(4))).^(-1)) - vY;
    
    options = optimoptions("lsqnonlin","Algorithm",'levenberg-marquardt');
               %"MaxFunctionEvaluations",10^(100), ...
               %"MaxIterations",10^(4));
               %"Display","iter-detailed", ...
               %"OptimalityTolerance",10^(-10), ...
               %"FunctionTolerance",10^(-10), ...
               %"StepTolerance",10^(-10), ...
    
    % x_1 denotes the coefficient vector obtained from nonlinear fit.
    [x_1,resnorm] = lsqnonlin(fun1,guess,[],[],options);
    syms t;
    figure(2)
    scatter(vX,vY);
    hold on
    fplot(x_1(1) + (x_1(2)*t)*(c_1+c_2*(x_1(3)/t)...
        +c_3*(x_1(3)/t)^(2)+c_4*(x_1(3)/t)^(3)...
        +c_5*(x_1(3)/t)^(4))/(c_1+c_6*(x_1(3)/t)...
        +c_7*(x_1(3)/t)^(2)+c_8*(x_1(3)/t)^(3)...
        +(x_1(3)/t)^(4)),[vX(1)-10, vX(length(vX))+10]);
    hold off
    
    % AICc calculation
    n = length(vX);
    sigmahat_sq = resnorm/n;
    K = length(x_1) + 1;
    AICc = n * log(sigmahat_sq) + 2*K*(n/(n-K-1));
end
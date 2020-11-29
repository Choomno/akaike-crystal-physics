function workflow(varargin)
    addpath('./');
    cInput = varargin;
    [vX,vY,guess,p] = input_data(cInput);
    [cModelsLin,cX] = ordinary_least_squares(vX,vY,p);
    cModelsNonlin = nonlinear_least_squares(vX,vY,guess);

    % p+2 due to the additional two models from nonlinear_least_squares
    cAICc_data = AICc_calc(cX,vY,cModelsLin,cModelsNonlin,p);

    % outputs results into xls files and pdf images
    output_data(cModelsLin,cModelsNonlin,cAICc_data,p);
end
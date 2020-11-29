% AICc Calculations

% figure vBeta
function cAICc_data = AICc_calc(cX,vY,cModelsLin,cModelsNonlin,p)
    for i = 1:p
        [vPoly_AICc(i),vs_sq_MLE(i)] = ...
            compute_poly_AICc(cX{i},vY,cModelsLin{1}{i});
    end
    vAICc = [vPoly_AICc,cModelsNonlin{2},cModelsNonlin{4}];
    
    % add two to number of models to account for Einstein and Debye fits.
    p = p+2;
    
    vAkaike_weights = Akaike_weights(vAICc,p);
    
    cEvidence_ratios = evidence_ratios(vAkaike_weights,p);
    
    cAICc_data = {vAICc,vAkaike_weights,cEvidence_ratios};
end

% Calculates the AICc, as well as the MLE estimate of the variance for the
% polynomial models.
function [poly_AICc,s_sq_MLE] = compute_poly_AICc(mX,vY,vBeta)

    n = length(vY);
    
    % 'i' is the number of rows in the matrix 'mX'.
    [i,~] = size(mX);
    
    % residuals
    vRes = vY - mX * vBeta;
    
    % residual sum of squares
    RSS = vRes' * vRes;
    
    % MLE estimate of the variance
    s_sq_MLE = RSS/i;
    
    % number of estimable parameters
    K = length(vBeta) + 1;
    
    % AICc computation
    poly_AICc = n * log(s_sq_MLE) + 2*K*(n/(n-K-1));
end

% Calculates Akaike weights. Takes a vector of AICc values and returns a
% vector of akaike weights corresponding, in identical order, the weights
% of each model in consideration.
function [vAkaike_weights] = Akaike_weights(vAICc,p)
    vDelta = [];
    
    vAkaike_weights = [];
    
    m = min(vAICc);
    
    % the difference in AIC value
    for i = 1:p
        vDelta(i) = vAICc(i) - m;
    end
    
    % Formula for Akaike weights as provided in Anderson's book
    denom = sum(exp(-0.5 * vDelta));
    for i = 1:p
        vAkaike_weights(i) = exp(-0.5 * vDelta(i))/denom;
    end
    
end

% Calculates evidence ratios from Akaike weights; takes a vector of Akaike
% weights and returns a vector of evidence ratios corresponding to the
% order of the models as it is in the input vector.
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
% Output data handling

function output_data(cModelsLin,cModelsNonlin,cAICc_data,p)
    cData_tables = {};
    cTemp = {};
    
    labelsModelsLin = {'Parameters', 'Std Error'};
    for i = 1:p
        cTemp{1} = transpose(cModelsLin{1}{i});
        cTemp{2} = transpose(cModelsLin{2}{i});
        cTemp = pad_cell(cTemp); % vectors need to be same size for table
        cData_tables{i} = table(cTemp{1},cTemp{2},'Variablenames',...
            labelsModelsLin);
    end
    
    labelsModelsNonlin = {'Einstein parameters', 'Debye Parameters'};
    
    cTemp = {};
    cTemp{1} = transpose(cModelsNonlin{1});
    cTemp{2} = transpose(cModelsNonlin{3});
    cTemp = pad_cell(cTemp);
    cData_tables{p+1} = table(cTemp{1},cTemp{2},'Variablenames',...
        labelsModelsNonlin);
    
    modelLabels = {};
    for i = 1:p
        modelLabels = [modelLabels, strcat('Poly_deg_',string(i))];
    end
    modelLabels{p+1} = 'Einstein';
    modelLabels{p+2} = 'Debye';
    modelLabels = transpose(modelLabels);
    
    labelsAICc = {'Model','AICc','Akaike_Weight','Std_Dev', 'R^2', ...
        'Adj_R^2', 'OLS_Var','Evidence_Ratios'};
    cTemp = {};
    cTemp{1} = transpose(cAICc_data{1});
    cTemp{2} = transpose(cAICc_data{2});
    for j = 3:6
        cTemp{j} = cModelsLin{j};
    end
    cTemp = pad_cell(cTemp);
    cTemp = [cTemp, {transpose(cAICc_data{3})}];
    cData_tables{p+2} = table(modelLabels,cTemp{1},cTemp{2},cTemp{3},cTemp{4},...
        cTemp{5},cTemp{6},cTemp{7},'VariableNames', labelsAICc);
    
    cFileNames = strings(p+2);
    for i=1:p
        cFileNames(i) = strcat('degree_',string(i), '.xls');
    end
    cFileNames = strings(p+2);
    for i=1:p
        cFileNames(i) = strcat('degree_',string(i), '.xls');
    end
    cFileNames(p+1) = 'nonlinear.xls';
    cFileNames(p+2) = 'AIC.xls';
    
    for i=1:p+2
        writetable(cData_tables{i}, cFileNames(i));
    end
    
end

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

% Outputs data as multiple .xls excel files.
% function [data_table] = output_data(cModelsLin,cModelsNonlin,cAICc_data)
%     cData_tables = {};
%     cRegression = {};
%     cAIC = {};
%     variable_names = {'Parameters','Standard_error', 'adjusted_R^2'};
%     for i = 1:p
%         cRegression{1,1} = cModelslin{1,i};
%         cRegression{1,2} = cModelslin{2,i};
%         cRegression{1,3} = vRsq_adj(i);
%         cRegression = pad_cell(cRegression); % vectors need to be same size
%         data_table = table(cRegression{1,1},cRegression{1,2}, ...
%                            cRegression{1,3},'Variablenames', ...
%                            variable_names);
%         cData_tables{1,i} = data_table;
%     end
%     
%     cAIC{1,1} = transpose(vAICc); % vectors are made to be columns
%     cAIC{1,2} = transpose(vAkaike_weights);
%     cAIC{1,3} = cEvidence_ratios{1,1};
%     cAIC = pad_cell(cAIC);
%     variable_names = {'AICc','Akaike_weights','Evidence_ratios'};
%     data_table = table(cAIC{1,1},cAIC{1,2}, cAIC{1,3},'Variablenames', ...
%                        variable_names);
%     cData_tables{1,p+1} = data_table;
%     for i = 1:p
%         name = strcat('degree_',string(i), '.xls');
%         writetable(cData_tables{1,i}, name);
%     end
%     writetable(data_table, 'polynomial_models.xls');
% end
function [X2, Os, Costdif, TotalError] = groupBrainSync(S)
    
    [numT, numV, SubNum] = size(S);

    % init for Os
    Os = zeros(numT, numT, SubNum);
    for i = 1:SubNum  %initializeing O
        R = 2 * rand(numT, numT) - 1; %define a random matrix with unity distributian from -1 to 1
        R = (R * R')^(-1/2) * R;  %orthogonal rows of matrix
        Os(:, :, i) = R';
    end
    
    Error = 1;
    PreError = 1;
    relcost = 1;
    
    alpha = 1e-6;
    var = 0;
    Costdif = zeros(10000, 1);
    
    disp('init done');
    
    while relcost > alpha
        var = var + 1;
        
        for i = 1:SubNum
            disp('subject iteration');
            disp(i);
            X = zeros(numT, numV);
            for j = 1:SubNum  %calculate X
                if j ~= i
                    X = Os(:, :, j) * S(:, :, j) + X;
                end
            end
            
            Y = X * (sqrt(1 ./ (SubNum * (SubNum-1))));
            scale = 1 ./ sqrt((SubNum-1)./SubNum); 
            
            [U, ~, V] = svd(Y * S(:, :, i)'); 
            Os(:, :, i) = scale * U * V'; % estimating trasfer matrices from svd of XY'      
        end
        
        disp('calculate error');
        Error = 0;
        X2 = (X + Os(:, :, i) * S(:, :, i)) ./ SubNum;
        
        for j = 1:SubNum
            etemp = Os(:, :, j) * S(:, :, j) - X2;
            Error = Error + trace(etemp*etemp'); %calculating error
        end
        
%         Av(var,:,:) = X2;
      
        Costdif(var) = PreError - Error;
        
        relcost = abs(Error - PreError) ./ abs(PreError);
        
        PreError = Error;
        
        var
        relcost
    end
    
    Costdif(var+1:end) = [];
    Costdif = Costdif(2:end);
    TotalError = Error;
    
end

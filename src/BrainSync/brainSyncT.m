function [Y2, R] = brainSyncT(X, Y)

    X = X'; Y = Y';
    [Y2, R] = brainSync(X, Y);
    Y2 = Y2';
    
end

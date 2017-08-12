Usage:

LBOL = getLaplaceBeltramiOperator(surfaceL, 1000);
LBOR = getLaplaceBeltramiOperator(surfaceR, 1000);

[dataSmL, dataSmR] = laplaceBeltramiSmoothFMRI(dataL, dataR, LBOL, LBOR, sigmaLB);

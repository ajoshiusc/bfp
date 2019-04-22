% Make the mexhyp2f1 mex file;

% Siyi Deng; 10-10-2013;

mex -I. -largeArrayDims -c round.c
mex -I. -largeArrayDims -c gamma.c
mex -I. -largeArrayDims -c polevl.c
mex -I. -largeArrayDims -c psi.c
mex -I. -largeArrayDims -c hyp2f1.c


if ispc
    mex -I. -largeArrayDims mexhyp2f1.c gamma.obj hyp2f1.obj polevl.obj psi.obj round.obj
    delete *.obj
else
    mex -I. -largeArrayDims mexhyp2f1.c gamma.o hyp2f1.o polevl.o psi.o round.o
    delete *.o
end
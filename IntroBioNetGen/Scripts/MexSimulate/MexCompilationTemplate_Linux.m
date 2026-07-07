% MexConstructionTemplate.m
% See the Faeder Lab document: CVODE_MEXCompilationOfBNGModels.pdf for 
% additional details

% Model file:
    modelcvode = 'Example_SimpleSTAT_cvode.c';

% Standard Linux Installation Paths:
    % specify_compiler        = 'GCC="/usr/bin/gcc-11"'; %Note: The default gcc compiler is Linux distribution dependent
    % specify_include_path    = '-I/usr/local/include';
    % specify_lib_path        = '-L/usr/local/lib';

% My HPC Cluster Installation Paths: 
    specify_compiler        = 'GCC="/usr/bin/gcc-11"'; %Note: The default gcc compiler is Linux distribution dependent
    specify_include_path    = '-I/opt/sundials-2.4.0/include';
    specify_lib_path        = '-L/opt/sundials-2.4.0/lib';

% Create mex file 
% Note: This file is computer specific and must be created with the computer on which it will be used
mex(modelcvode,specify_compiler,specify_include_path,specify_lib_path,'-lsundials_cvode','-lsundials_nvecserial','-lm')
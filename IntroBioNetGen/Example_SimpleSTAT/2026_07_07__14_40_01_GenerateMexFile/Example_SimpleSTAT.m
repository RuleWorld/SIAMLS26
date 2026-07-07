function [err, timepoints, species_out, observables_out ] = Example_SimpleSTAT( timepoints, species_init, parameters, suppress_plot )
%EXAMPLE_SIMPLESTAT Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'Example_SimpleSTAT' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the CVode library interfaced
%   to MATLAB via the MEX interface. Before running this script, the model
%   source in file Example_SimpleSTAT_cvode.c must be compiled (see that file for details).
%   EXAMPLE_SIMPLESTAT returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = Example_SimpleSTAT( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   timepoints      : column vector of time points returned by integrator.
%   species_init    : row vector of 17 initial species populations.
%   parameters      : row vector of 26 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 0.200815, 0.00000774704, 3.15924, 0.00000361185, 1003.25, 0.0000476335, 1.08116, 0.0000000197800, 9.23454, 0.00000415259, 1.77990, 3.16032, 0.046620, 0.111965, 0.00556053, 11.00256433, 1, 2.12893, 5889.76, 6.83256, 0, 1, 1, 10.2990, 48.2472, 929991 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 26  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 26].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 17  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 17].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,90,91+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'il10_il10r1_binding', 'il10_il10r1_unbinding', 'il10r1_il10r2_binding', 'il10r1_il10r2_unbinding', 'il10_complex_jak1_binding', 'il10_complex_jak1_unbinding', 'il10_jak1_med_STAT3_act', 'il10_jak1_med_STAT1_act', 'SOCS1_jak1_binding', 'SOCS1_jak1_unbinding', 'pSTAT3_rec_dissoc', 'pSTAT1_rec_dissoc', 'PTP_med_STAT3_deact', 'PTP_med_STAT1_deact', 'STAT3_SOCS1_ind', 'STAT1_SOCS1_ind', 'IL10_0', 'IL10R1_0', 'IL10R2_0', 'JAK1_0', 'SOCS1_0', 'PTP3_0', 'PTP1_0', 'SOCS1_degrad', 'STAT3_0', 'STAT1_0' };



%% Integrate Network Model
try 
    % run simulation
    [err, species_out, observables_out] = Example_SimpleSTAT_cvode( timepoints, species_init, parameters );
catch
    fprintf( 1, 'Error: some problem integrating ODE network! (CVODE exitflag %d)\n', err );
    err = 1;
    return;
end



%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'total_pSTAT3', 'total_pSTAT1', 'total_pSTAT3_explicitdefinition', 'total_pSTAT1_explicitdefinition', 'total_pSTAT3_species', 'total_pSTAT1_species' };

    % construct figure
    plot(timepoints,observables_out);
    title('Example_SimpleSTAT observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end



%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%



% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,17);
    species_init(1) = params(17);
    species_init(2) = params(18);
    species_init(3) = params(19);
    species_init(4) = params(20);
    species_init(5) = params(21);
    species_init(6) = params(22);
    species_init(7) = params(23);
    species_init(8) = params(25);
    species_init(9) = params(26);
    species_init(10) = 0;
    species_init(11) = 0;
    species_init(12) = 0;
    species_init(13) = 0;
    species_init(14) = 0;
    species_init(15) = 0;
    species_init(16) = 0;
    species_init(17) = 0;

end


end

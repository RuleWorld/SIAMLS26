function [err, timepoints, species_out, observables_out] = Example_SimpleSTAT( timepoints, species_init, parameters, suppress_plot )
%EXAMPLE_SIMPLESTAT Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'Example_SimpleSTAT' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. EXAMPLE_SIMPLESTAT returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = Example_SimpleSTAT( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 17 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
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
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...
               'AbsTol',   1e-10,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 6 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
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

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

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


% user-defined functions



% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,26);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = parameters(8);
    expressions(9) = parameters(9);
    expressions(10) = parameters(10);
    expressions(11) = parameters(11);
    expressions(12) = parameters(12);
    expressions(13) = parameters(13);
    expressions(14) = parameters(14);
    expressions(15) = parameters(15);
    expressions(16) = parameters(16);
    expressions(17) = parameters(17);
    expressions(18) = parameters(18);
    expressions(19) = parameters(19);
    expressions(20) = parameters(20);
    expressions(21) = parameters(21);
    expressions(22) = parameters(22);
    expressions(23) = parameters(23);
    expressions(24) = parameters(24);
    expressions(25) = parameters(25);
    expressions(26) = parameters(26);
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,6);
    observables(1) = species(14) +species(16);
    observables(2) = species(15) +species(17);
    observables(3) = species(14) +species(16);
    observables(4) = species(15) +species(17);
    observables(5) = species(14) +species(16);
    observables(6) = species(15) +species(17);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,6);
    ratelaws(1) = expressions(1)*species(2)*species(1);
    ratelaws(2) = expressions(24)*species(5);
    ratelaws(3) = expressions(2)*species(10);
    ratelaws(4) = expressions(3)*species(10)*species(3);
    ratelaws(5) = expressions(4)*species(11);
    ratelaws(6) = expressions(5)*species(11)*species(4);
    ratelaws(7) = expressions(6)*species(12);
    ratelaws(8) = expressions(9)*species(12)*species(5);
    ratelaws(9) = expressions(7)*species(12)*species(8);
    ratelaws(10) = expressions(8)*species(12)*species(9);
    ratelaws(11) = expressions(10)*species(13);
    ratelaws(12) = expressions(11)*species(14);
    ratelaws(13) = expressions(12)*species(15);
    ratelaws(14) = expressions(13)*species(6)*species(16);
    ratelaws(15) = expressions(14)*species(7)*species(17);
    ratelaws(16) = expressions(15)*species(16);
    ratelaws(17) = expressions(16)*species(17);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(17,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -ratelaws(1) +ratelaws(3);
    Dspecies(2) = -ratelaws(1) +ratelaws(3);
    Dspecies(3) = -ratelaws(4) +ratelaws(5);
    Dspecies(4) = -ratelaws(6) +ratelaws(7);
    Dspecies(5) = -ratelaws(2) -ratelaws(8) +ratelaws(11) +ratelaws(16) +ratelaws(17);
    Dspecies(6) = 0.0;
    Dspecies(7) = 0.0;
    Dspecies(8) = -ratelaws(9) +ratelaws(14);
    Dspecies(9) = -ratelaws(10) +ratelaws(15);
    Dspecies(10) = ratelaws(1) -ratelaws(3) -ratelaws(4) +ratelaws(5);
    Dspecies(11) = ratelaws(4) -ratelaws(5) -ratelaws(6) +ratelaws(7);
    Dspecies(12) = ratelaws(6) -ratelaws(7) -ratelaws(8) -ratelaws(9) -ratelaws(10) +ratelaws(11) +ratelaws(12) +ratelaws(13);
    Dspecies(13) = ratelaws(8) -ratelaws(11);
    Dspecies(14) = ratelaws(9) -ratelaws(12);
    Dspecies(15) = ratelaws(10) -ratelaws(13);
    Dspecies(16) = ratelaws(12) -ratelaws(14);
    Dspecies(17) = ratelaws(13) -ratelaws(15);

end


end

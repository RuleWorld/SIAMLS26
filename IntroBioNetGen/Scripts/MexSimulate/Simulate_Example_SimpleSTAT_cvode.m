% Timepoints
timepoints=0:1:90;
timepoints=timepoints';

% Parameters
il10_il10r1_binding          = 0.200815;  
il10_il10r1_unbinding        = 0.00000774704;  
il10r1_il10r2_binding        = 3.15924;  
il10r1_il10r2_unbinding      = 0.00000361185;  
il10_complex_jak1_binding    = 1003.25;  
il10_complex_jak1_unbinding  = 0.0000476335;  
il10_jak1_med_STAT3_act      = 1.08116;  
il10_jak1_med_STAT1_act      = 0.0000000197800;  
SOCS1_jak1_binding           = 9.23454;  
SOCS1_jak1_unbinding         = 0.00000415259;  
pSTAT3_rec_dissoc            = 1.77980;  
pSTAT1_rec_dissoc            = 3.16032;  
PTP_med_STAT3_deact          = 0.0466199;  
PTP_med_STAT1_deact          = 0.111965;  
STAT3_SOCS1_ind              = 0.00556053;  
STAT1_SOCS1_ind              = 11.0026;  
IL10_0                       = 10; %1;  
IL10R1_0                     = 2.12893;  
IL10R2_0                     = 5889.76;  
JAK1_0                       = 6.83256;  
SOCS1_0                      = 0;  
PTP3_0                       = 1;  
PTP1_0                       = 1;  
SOCS1_degrad                 = 10.2990;  
STAT3_0                      = 48.2472;  
STAT1_0                      = 929991.67;  

parameters(01) = il10_il10r1_binding;
parameters(02) = il10_il10r1_unbinding;
parameters(03) = il10r1_il10r2_binding;  
parameters(04) = il10r1_il10r2_unbinding;
parameters(05) = il10_complex_jak1_binding; 
parameters(06) = il10_complex_jak1_unbinding;
parameters(07) = il10_jak1_med_STAT3_act;
parameters(08) = il10_jak1_med_STAT1_act;
parameters(09) = SOCS1_jak1_binding;
parameters(10) = SOCS1_jak1_unbinding;
parameters(11) = pSTAT3_rec_dissoc;
parameters(12) = pSTAT1_rec_dissoc; 
parameters(13) = PTP_med_STAT3_deact;
parameters(14) = PTP_med_STAT1_deact;
parameters(15) = STAT3_SOCS1_ind;
parameters(16) = STAT1_SOCS1_ind;
parameters(17) = IL10_0;
parameters(18) = IL10R1_0;
parameters(19) = IL10R2_0;
parameters(20) = JAK1_0;
parameters(21) = SOCS1_0;
parameters(22) = PTP3_0; 
parameters(23) = PTP1_0;
parameters(24) = SOCS1_degrad;
parameters(25) = STAT3_0; 
parameters(26) = STAT1_0;

% Initial Conditions
species_init(01) = IL10_0;  % IL10(il10r1) IL10_0
species_init(02) = IL10R1_0;% IL10R1(il10r2,jak1,l2,stat) IL10R1_0
species_init(03) = IL10R2_0;% IL10R2(il10r1) IL10R2_0
species_init(04) = JAK1_0;  % JAK1(rec,socs1) JAK1_0
species_init(05) = SOCS1_0; % SOCS1(jak1) SOCS1_0
species_init(06) = PTP3_0;  % PTP3() PTP3_0
species_init(07) = PTP1_0;  % PTP1() PTP1_0
species_init(08) = STAT3_0; % STAT3(Y~0) STAT3_0
species_init(09) = STAT1_0; % STAT1(Y~0) STAT1_0
species_init(10) = 0;       % IL10(il10r1!1).IL10R1(il10r2,jak1,l2!1,stat) 0
species_init(11) = 0;       % IL10(il10r1!1).IL10R1(il10r2!2,jak1,l2!1,stat).IL10R2(il10r1!2) 0
species_init(12) = 0;       % IL10(il10r1!1).IL10R1(il10r2!2,jak1!3,l2!1,stat).IL10R2(il10r1!2).JAK1(rec!3,socs1) 0
species_init(13) = 0;       % IL10(il10r1!1).IL10R1(il10r2!2,jak1!3,l2!1,stat).IL10R2(il10r1!2).JAK1(rec!3,socs1!4).SOCS1(jak1!4) 0
species_init(14) = 0;       % IL10(il10r1!1).IL10R1(il10r2!2,jak1!3,l2!1,stat!4).IL10R2(il10r1!2).JAK1(rec!3,socs1).STAT3(Y~P!4) 0
species_init(15) = 0;       % IL10(il10r1!1).IL10R1(il10r2!2,jak1!3,l2!1,stat!4).IL10R2(il10r1!2).JAK1(rec!3,socs1).STAT1(Y~P!4) 0
species_init(16) = 0;       % STAT3(Y~P) 0
species_init(17) = 0;       % STAT1(Y~P) 0

% Simulate and generate plot
species_init(01) = 10;  % IL10(il10r1) IL10_0
[err, species_out, observables_out] = Example_SimpleSTAT_cvode(timepoints, species_init, parameters);

gemColors = orderedcolors("gem"); 
figure
plot(timepoints,observables_out(:,1),'Color',gemColors(1,:),'LineWidth',3)
title('pSTAT3 dynamics in response to 10 units of IL10')
xlabel('Time (min)')
ylabel('Concentration')
legend('pSTAT3')

figure
plot(timepoints,observables_out(:,2),'Color',gemColors(2,:),'LineWidth',3)
title('pSTAT1 dynamics in response to 10 units of IL10')
xlabel('Time (min)')
ylabel('Concentration')
legend('pSTAT1')

% Simulate and generate plot
species_init(01) = 1;  % IL10(il10r1) IL10_0
[err, species_out, observables_out] = Example_SimpleSTAT_cvode(timepoints, species_init, parameters);

gemColors = orderedcolors("gem"); 
figure
plot(timepoints,observables_out(:,1),'Color',gemColors(1,:),'LineWidth',3)
title('pSTAT3 dynamics in response to 1 units of IL10')
xlabel('Time (min)')
ylabel('Concentration')
legend('pSTAT3')

figure
plot(timepoints,observables_out(:,2),'Color',gemColors(2,:),'LineWidth',3)
title('pSTAT1 dynamics in response to 1 units of IL10')
xlabel('Time (min)')
ylabel('Concentration')
legend('pSTAT1')
%%
%%% ======================================================================
%%% run_hcmaq.m
%%% Created by Sina Voshtani
%%% Created on 21/09/2021
%%% =======================================================================

%%% make environment variables 
kk = num2str(i-1);
setenv('i',kk);
ll = num2str((nstep-1)*(24/step_size));
setenv('j',ll);

%%% call CMAQ for concentrations run
system('sh /home/sinavo/pkg/cmaq/5.3/CCTM/da_9/conc/fwd_hrda.sh') % path to runscript
disp('*** successfully run cmaq concentration ***')

%%% call CMAQ for error variance run
system('sh /home/sinavo/pkg/cmaq/5.3/CCTM/da_9/vari/fwd_hrda.sh')% path to runscript
disp('*** successfully run cmaq variance ***')
%%% =======================================================================
%%% END
%%% =======================================================================


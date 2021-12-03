%% 
%%% ======================================================================
%%% cv_make_3fold.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% Clear workspace and command window 
clear
clc

%%% Header
fprintf('\n ----------------------------------------\n')
fprintf(' --- Start of 3-fold cross-validation ---\n')
fprintf(' ----------------------------------------\n')

n = 1:672; % number of timesteps
% load Bias_corrected_gosat_AMF_201004.mat
% load BI_COR_GOSAT_AMF_201004.mat
load  gosat_cmaq_biascorrected_qc_on_latitude.mat % load bias-corrected data

%%% sort satellite observations data
lat_sat2  = nonzeros(lat_sat_all_bc(:,n));
lon_sat2  = nonzeros(lon_sat_all_bc(:,n));
tim_sat2  = nonzeros(tim_sat_all_bc(:,n));
ch4_sat2  = nonzeros(ch4_sat_all_bc(:,n));
sig_sat2  = nonzeros(sig_sat_all_bc(:,n));

plv_sat_w = reshape(plv_sat_all_bc,5,[]);
pwe_sat_w = reshape(pwe_sat_all_bc,4,[]);
avk_sat_w = reshape(avk_sat_all_bc,4,[]);
fgs_sat_w = reshape(fgs_sat_all_bc,4,[]);
all0     = all(plv_sat_w == 0);
plv_sat2 = plv_sat_w(:, ~all0);
pwe_sat2 = pwe_sat_w(:, ~all0);
avk_sat2 = avk_sat_w(:, ~all0);
fgs_sat2 = fgs_sat_w(:, ~all0);

%%% filtering observation data (thinning)- there is no thinning by default
% loc_r        = [lat_sat1,lon_sat1];
% loc_rr       = round(loc_r,3);   
% [~,ia]       = unique(loc_rr,'rows','stable');
% ch4_sat2  = ch4_sat1(ia);
% sig_sat2  = sig_sat1(ia);
% lat_sat2  = lat_sat1(ia);
% lon_sat2  = lon_sat1(ia);
% tim_sat2  = tim_sat1(ia);
% plv_sat2  = plv_sat1(:,ia);
% pwe_sat2  = pwe_sat1(:,ia);
% fgs_sat2  = fgs_sat1(:,ia);
% avk_sat2  = avk_sat1(:,ia);

%%% splitting into 3 fold 
kfold      = 3;
idx_all   = 1:1:numel(ch4_sat2);
if numel(ch4_sat2) >= 2
    idx_test  = 1:kfold:numel(ch4_sat2);
else
    idx_test= [];
end

%%% create train set (2/3 of data) and test set(1/3 of data) 
idx_train = setdiff(idx_all,idx_test);
ch4_sat_train  = ch4_sat2(idx_train);
sig_sat_train  = sig_sat2(idx_train);
lat_sat_train  = lat_sat2(idx_train);
lon_sat_train  = lon_sat2(idx_train);
tim_sat_train  = tim_sat2(idx_train);
plv_sat_train  = plv_sat2(:,idx_train);
pwe_sat_train  = pwe_sat2(:,idx_train);
fgs_sat_train  = fgs_sat2(:,idx_train);
avk_sat_train  = avk_sat2(:,idx_train);
ch4_sat_1 = ch4_sat_train(1:2:end);
lat_sat_1 = lat_sat_train(1:2:end);
lon_sat_1 = lon_sat_train(1:2:end);
tim_sat_1 = tim_sat_train(1:2:end);
ch4_sat_2 = ch4_sat_train(2:2:end);
lat_sat_2 = lat_sat_train(2:2:end);
lon_sat_2 = lon_sat_train(2:2:end);
tim_sat_2 = tim_sat_train(2:2:end);

ch4_sat_test  = ch4_sat2(idx_test);
sig_sat_test  = sig_sat2(idx_test);
lat_sat_test  = lat_sat2(idx_test);
lon_sat_test  = lon_sat2(idx_test);
tim_sat_test  = tim_sat2(idx_test);
plv_sat_test  = plv_sat2(:,idx_test);
pwe_sat_test  = pwe_sat2(:,idx_test);
fgs_sat_test  = fgs_sat2(:,idx_test);
avk_sat_test  = avk_sat2(:,idx_test);
t_sat_train   = datetime(tim_sat_train,'ConvertFrom','posixtime');
t_sat_test    = datetime(tim_sat_test,'ConvertFrom','posixtime');
t_sat_train1  = datetime(tim_sat_1,'ConvertFrom','posixtime');
t_sat_train2  = datetime(tim_sat_2,'ConvertFrom','posixtime');

%%% sort each set of data (train and test) in an hourly order form reqired
%%% by PvKF assimilation
step_size = 24;
for i = 1:4 % change 4 to 28 for four weeks simulation of a month
    t1       = datetime(2010,04,i+1,0,0,0);
    t_step   = dateshift(t1,'start','hour',0:1:24);
    for nstep = 1:1:24
        s = 1;
        for k=1:1:length(t_sat_test)
            if t_step(nstep)<t_sat_test(k) && t_sat_test(k)<t_step(nstep+1)
                %lat1(s)  = lat_sat_train(k);
                %lon1(s)  = lon_sat_train(k);
                %tim1(s)  = tim_sat_train(k);
                %ch41(s)  = ch4_sat_train(k);
                %sig1(s)  = sig_sat_train(k);
                %avk1(:,s)  = avk_sat_train(:,k);
                %plv1(:,s)  = plv_sat_train(:,k);
                %pwe1(:,s)  = pwe_sat_train(:,k);
                %fgs1(:,s)  = fgs_sat_train(:,k);
                lat1(s)   = lat_sat_test(k);
                lon1(s)   = lon_sat_test(k);
                tim1(s)   = tim_sat_test(k);
                ch41(s)   = ch4_sat_test(k);
                sig1(s)   = sig_sat_test(k);
                avk1(:,s) = avk_sat_test(:,k);
                plv1(:,s) = plv_sat_test(:,k);
                pwe1(:,s) = pwe_sat_test(:,k);
                fgs1(:,s) = fgs_sat_test(:,k);
                s=s+1;
            end
        end
        if (~exist('ch41','var')) || (isempty(ch41))
            lat1=[];lon1=[];tim1=[];ch41=[];sig1=[];avk1=[];plv1=[];pwe1=[];fgs1=[];
        end
         %lat_sat_all_train(1:length(lat1),step_size*(i-1)+nstep)   = lat1;
         %lon_sat_all_train(1:length(lat1),step_size*(i-1)+nstep)   = lon1;
         %tim_sat_all_train(1:length(lat1),step_size*(i-1)+nstep)   = tim1;
         %ch4_sat_all_train(1:length(lat1),step_size*(i-1)+nstep)   = ch41;
         %sig_sat_all_train(1:length(lat1),step_size*(i-1)+nstep)   = sig1;
         %plv_sat_all_train(1:5,1:length(lat1),step_size*(i-1)+nstep) = plv1;
         %pwe_sat_all_train(1:4,1:length(lat1),step_size*(i-1)+nstep) = pwe1;
         %fgs_sat_all_train(1:4,1:length(lat1),step_size*(i-1)+nstep) = fgs1;
         %avk_sat_all_train(1:4,1:length(lat1),step_size*(i-1)+nstep) = avk1;
         
        lat_sat_all_test(1:length(lat1),step_size*(i-1)+nstep)   = lat1;
        lon_sat_all_test(1:length(lat1),step_size*(i-1)+nstep)   = lon1;
        tim_sat_all_test(1:length(lat1),step_size*(i-1)+nstep)   = tim1;
        ch4_sat_all_test(1:length(lat1),step_size*(i-1)+nstep)   = ch41;
        sig_sat_all_test(1:length(lat1),step_size*(i-1)+nstep)   = sig1;
        plv_sat_all_test(1:5,1:length(lat1),step_size*(i-1)+nstep) = plv1;
        pwe_sat_all_test(1:4,1:length(lat1),step_size*(i-1)+nstep) = pwe1;
        fgs_sat_all_test(1:4,1:length(lat1),step_size*(i-1)+nstep) = fgs1;
        avk_sat_all_test(1:4,1:length(lat1),step_size*(i-1)+nstep) = avk1;
%         nstep
        lat1=[];lon1=[];tim1=[];ch41=[];sig1=[];avk1=[];plv1=[];pwe1=[];fgs1=[];
    end
    i
end
%%% save train and test file into current directory in an hourly form -
%%% uncomment for generating your own use of train and test sets
% % save('cv_input_train23_all.mat','lat_sat_all_train','lon_sat_all_train',...
% %     'ch4_sat_all_train','sig_sat_all_train','plv_sat_all_train','pwe_sat_all_train',...
% %     'fgs_sat_all_train','avk_sat_all_train','tim_sat_all_train')
% % save('cv_input_test1_all.mat','lat_sat_all_test','lon_sat_all_test',...
% %     'ch4_sat_all_test','sig_sat_all_test','plv_sat_all_test','pwe_sat_all_test',...
% %     'fgs_sat_all_test','avk_sat_all_test','tim_sat_all_test')

%%% if do not generate your own train and test sets, then two below 
load cv_input_test1_all.mat
load cv_input_train23_all.mat

%%% ploting the train and test dataset over Northern hemisphere
n=1*24:24*4;
lat_11  = nonzeros(lat_sat_all_train(:,n));
lon_11  = nonzeros(lon_sat_all_train(:,n));
ch4_11  = nonzeros(ch4_sat_all_train(:,n));
lat_22  = nonzeros(lat_sat_all_test(:,n));
lon_22  = nonzeros(lon_sat_all_test(:,n));
ch4_22  = nonzeros(ch4_sat_all_test(:,n));
% tim_22  = nonzeros(tim_sat_all_train(:,n));
% tq       = datetime(tim_11,'ConvertFrom','posixtime');

figure(4)
% worldmap('World')
load coastlines
axesm('stereo','Origin',[90 -98],'MapLatLimit',[1 90])
axis off
framem on
gridm on
mlabel on
plabel on;
setm(gca,'MLabelParallel',0)
plotm(coastlat,coastlon)

scatterm(lat_11(:),lon_11(:),10,ch4_11(:),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
scatterm(lat_22(:),lon_22(:),10,ch4_22(:),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
% % scatterm(lat_sat_test(en),lon_sat_test(en),20,ch4_sat_test(en),'Marker','s','MarkerFaceColor',gg(7,:),'MarkerEdgeColor',gg(7,:));
% % hold on
% scatterm(lat_sat_test,lon_sat_test,30,ch4_sat_test,'MarkerFaceColor',bb(6,:),'MarkerEdgeColor','k');
% title('Innovation in 1 hour $C_{v}$: SOAR ($L_{\sigma}= 1 \sigma)$',...
%     'FontSize',20,'interpreter','latex') % - $ppm$
title('3-fold , 1 hour , all data',...
    'FontSize',20,'interpreter','latex') % - $ppm$
% colormap(parula(10))
% % % % h =colorbar;%('location', 'SouthOutside','FontSize',20);
% % % % colorTitleHandle = get(h,'Title');
% % % % set(colorTitleHandle ,'String','$ppm$','interpreter','latex');
% % % % colormap(brewermap(10,'*Spectral'))
set(gca,'FontSize',20)
caxis([-0.2 0.2])
set(gca,'FontSize',40)
hold on
h = zeros(2, 1);
h(1) = plot(nan,nan,'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','b');
h(2) = plot(nan,nan,'LineStyle','none','Marker','^','MarkerSize',10,'MarkerFaceColor','r');
% h(3) = plot(nan,nan,'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor',gg(7,:));
 legend(h,'set 1','set 2');
aa = [lat_sat2,lon_sat2];
[ll,ll_idx,qq] = unique(aa,'rows','stable');
ch4_new     = ch4_sat2(ll_idx);

%%% Finished simulation
fprintf('\n --------------------------------------\n')
fprintf(' --- End of 3-fold cross-validation ---\n')
fprintf(' --------------------------------------\n\n')
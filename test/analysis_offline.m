%%
%%% ======================================================================
%%% analysis_offline.m
%%% Created by Sina Voshtani
%%% Created on 21/09/2021
%%% =======================================================================

if ~isempty(ch4_sat)
    pf_ii        = c_var ;
    %fgs_var_sat  = fgs_sat * 1e-3 * 1e-6; % fgs is on ppb   

%%% covaraince modelling
    
%%% construction of a vertical correlation matrix: C_z
%%% there are three options: SOAR, FOAR, and Gaussian

    Lz = 1;
    alpha = 1; 
    for l1  = 1:numel(hcmaq_siglvl)-1
        for l2  = 1:numel(hcmaq_siglvl)-1
            %d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
            %    ((hcmaq_siglvl(1)-hcmaq_siglvl(end))); % 2-sigma correlation lenght
	        if (l2 > Lz) && (l2 < 45-Lz)
        	    d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
        	        abs(alpha*(hcmaq_siglvl(l2-Lz)-hcmaq_siglvl(l2+Lz)));
       		 elseif (l2 <= Lz)
            		d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
                	abs(2*alpha*(hcmaq_siglvl(l2)-hcmaq_siglvl(l2+Lz)));
	        elseif (l2 >= 45-Lz)
        	  	d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
	                abs(2*alpha*(hcmaq_siglvl(l2-Lz)-hcmaq_siglvl(l2)));
        	end
	        C_z = (1+d_z/1.3494) .* exp(-d_z/1.3494); %SOAR
               % C_z = exp(-d_z/0.5005); %FOAR -d_z/0.5005
               % C_z = exp(-d_z.^2/2); %Gaussian
        end
    end
    C_z = triu(C_z.',1) + tril(C_z);
    % C_z = toeplitz([1 0.5 zeros(1,numel(hcmaq_siglvl)-3)]);%Trapizoid-2L
    % C_z4 = toeplitz([1:-1/44:0 zeros(1,numel(hcmaq_siglvl)-46)]);%Trapizoid-2L

%%% construction of a horizontal correlation matrix: C_xy
%%% there are three options: SOAR, FOAR, and Gaussian    

    L = 5 * abs(xproj_hcmaq(1,1) - xproj_hcmaq(2,1)); %  X108000 km
    d_xy_o        = pdist2(r_hcmaq(:,:),r_sat(:,:),'euclidean');
    C_xy_1d       = exp(-((d_xy_o.^2)/(2*(L^2)))); %Gaussian
    %C_xy_1d       = (1+(d_xy_o/L)) .*  exp(-(d_xy_o/L)); %SOAR
    %C_xy_1d       = exp(-(d_xy_o/L)); %FOAR
    C_xy          = reshape(C_xy_1d,187,187,[]);
    %C_xy(C_xy<1e-1)=0;  
    %C_xy_hv    = reshape(C_xy,187^2,[]);
    C_xy_hh    = reshape(C_xy,187^2,[])';

%%% initialization of background error covaraince matrix 
    pp_o_avg   = zeros([numel(xproj_hcmaq),44, numel(x_sat)]);
    pf_ii_hv   = reshape(pf_ii,187^2,44);
    pf_ii_hh   = reshape(pf_ii(:,:,1),1,[]);
    pf_ii_vec  = reshape(pf_ii,1,44*187^2);

%%% varaince on observation location - model levels
    for l1 = 1:44
         H_hor_pf        = griddedInterpolant(xproj_hcmaq,yproj_hcmaq,pf_ii(:,:,l1));
         %pf_ii_m2o(l1,:)= H_hor_pf(x_sat(:),y_sat(:));
         pf_ii_m2o(l1,:) = H_hor_pf(xproj_sat(:),yproj_sat(:));
    end

%%% pressure on observation location - model levels
    for l1 = 1:45
         H_hor_pres       = griddedInterpolant(xproj_hcmaq,yproj_hcmaq,hcmaq_p(:,:,l1));
         pres_m2o_l(l1,:) = H_hor_pres(xproj_sat(:),yproj_sat(:));
         %pres_m2o_l(l1,:) = H_hor_pres(x_sat(:),y_sat(:));
    end
    
%%% indexing matching levels betweeen model and observations
    num_l     = 1:45; %numel(pres_m2o_l(:,p));
    % avk_sat   = avk_sat*0 +1;
    for p = 1:numel(x_sat(:))   
        indx_l(:,p)    = interp1(pres_m2o_l(:,p),num_l,plv_sat(:,p),'linear','extrap');
    end
    indx_l(indx_l<1) = 1;
    indx_l(indx_l>45) = 45;
    low_l = floor(indx_l);
    % up_l = floor(indx_l)+1;
    % up_l(up_l>44)=44; up_l(up_l<1)=1; 
    low_l(low_l>44)=44; low_l(low_l<1)=1; 
    indx_l(1,:) = 1; inxd_l(numel(plv_sat(:,1)),:) = 45;
    low_l(1,:)  = 1; low_l(numel(plv_sat(:,1)),:)  = 44;
    % up_l(1,:)   = 1; up_l(numel(plv_sat(:,1)),:)  = 44;   

%%% compute pressure weight on model space 
    % pwe_m  = zeros(numel(indx_l(:,1))-1,numel(x_sat));
    for p=1:numel(x_sat(:))
        for l2 = 1:numel(indx_l(:,1))-1
            pwe_m{p,l2} = (pres_m2o_l(low_l(l2,p):low_l(l2+1,p),p)-pres_m2o_l(low_l(l2,p)+1:low_l(l2+1,p)+1,p))...
                        ./ (pres_m2o_l(low_l(l2,p),p)-pres_m2o_l(low_l(l2+1,p)+1,p));
        end
    end
    %%% uniform veritcal 
    for l2 = 1:numel(indx_l(:,1))
        low(l2) = mode(low_l(l2,:));
    end
    for l2 = 1:numel(indx_l(:,1))-1
        
        pwe_m1{l2} = (pres_m2o_l(low(l2):low(l2+1),p)-pres_m2o_l(low(l2)+1:low(l2+1)+1,p))...
                        ./ (pres_m2o_l(low(l2),p)-pres_m2o_l(low(l2+1)+1,p));
    end 

%%% Computation of model column averaged (before get an analysis)    
    xf = c_con; % ppm
    xf_hv = reshape(xf,numel(xproj_hcmaq),44); %% 3D concentrations    
    fgs_x_sat  = fgs_sat * 1e-3; % fgs is on ppb that's why multiply by 1e-3

%%% horizontal intepolation    
    for l1 = 1:44
        H_xf        = griddedInterpolant(xproj_hcmaq,yproj_hcmaq,xf(:,:,l1));
        xf_m2o(l1,:)  = H_xf(xproj_sat(:),yproj_sat(:));
    end
    xf_o = zeros([numel(indx_l(:,1))-1,numel(x_sat)]);
    fgs_x  = (( ones(numel(indx_l(:,1))-1,numel(x_sat)) - avk_sat) .* fgs_x_sat);
    %%% find a corresponding equivalent layer in model from observation
    for p = 1:numel(x_sat)
        for l2 = 1:numel(indx_l(:,1))-1
        	xf_o(l2,p) = (1*xf_m2o(low_l(l2,p):low_l(l2+1,p),p)'...
            	* pwe_m{p,l2}(1:end));
		% xf_o(l2,p) = (1*xf_m2o(low(l2):low(l2+1),p)'...
		% * pwe_m1{l2}(1:end));
        end
    
   	fgs_x_3d =  fgs_x(:,p);
    	avk_x_3d =  avk_sat(:,p);
    	x_o_avg_l  =  fgs_x_3d + avk_x_3d .* xf_o(:,p); %fgs_ph_3d +
    	xf_o_avg(p) = x_o_avg_l' * pwe_sat(:,p); %this is X-avg for each pi
    end    
    
else 
    xf_3d = c_con; % = xf_3d
    pf_3d = c_var; % = pf_3d
    xa_3d = c_con; % = xf_3d
    pa_3d = c_var; % = pf_3d
    xd_3d = 0;
    pd_3d = 0;
    innov = 0;
    chi_sq= 0;
    x_sat = [];
    xf_o_avg = [];
end
%%% store main variables for each timestep
%innov_all(1:length(innov),step_size*(i-1)+nstep) = innov;
lat_sat_all(1:length(lat_sat),step_size*(i-1)+nstep)   = lat_sat;
lon_sat_all(1:length(lon_sat),step_size*(i-1)+nstep)   = lon_sat;
tim_sat_all(1:length(lon_sat),step_size*(i-1)+nstep)   = tim_sat;
ch4_sat_all(1:length(lat_sat),step_size*(i-1)+nstep)   = ch4_sat;
sig_sat_all(1:length(lat_sat),step_size*(i-1)+nstep)   = sig_sat;
sen_sat_all(1:length(lat_sat),step_size*(i-1)+nstep)   = sen_sat;
sol_sat_all(1:length(lat_sat),step_size*(i-1)+nstep)   = sol_sat;
plv_sat_all(1:size(plv_sat,1),1:length(lat_sat),step_size*(i-1)+nstep) = plv_sat;
pwe_sat_all(1:size(pwe_sat,1),1:length(lat_sat),step_size*(i-1)+nstep) = pwe_sat;
fgs_sat_all(1:size(fgs_sat,1),1:length(lat_sat),step_size*(i-1)+nstep) = fgs_sat;
avk_sat_all(1:size(avk_sat,1),1:length(lat_sat),step_size*(i-1)+nstep) = avk_sat;

xm_o_avg_all(1:length(xf_o_avg),step_size*(i-1)+nstep) = xf_o_avg; % analysis
%xf_o_avg_all(1:length(xf_o_avg),step_size*(i-1)+nstep) = xf_o_avg; % background
%%%%var_A    = diag(hph);
%%%%var_A_all(1:length(var_A),step_size*(i-1)+nstep) = var_A; % diagonal of hph
%var_B    = diag(hph);
%var_B_all(1:length(var_B),step_size*(i-1)+nstep) = var_B; % diagonal of hph

%%% chi2 diagnostic metric
%chi_sq_all(step_size*(i-1)+nstep)       = chi_sq;
obs_num_all(step_size*(i-1)+nstep)      = numel(x_sat);
%chi_sq
obs_num = obs_num_all(end)
%norm_chi=chi_sq/obs_num
%pause(0.5);
%A1_cv_passive(1:length(xf_o_avg),step_size*(i-1)+nstep) = xf_o_avg;
%A1_cv_active(1:length(xf_o_avg),step_size*(i-1)+nstep)  = xf_o_avg;

%%% =======================================================================
%%% END
%%% =======================================================================

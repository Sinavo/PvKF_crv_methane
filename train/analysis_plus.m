%%
%%% ======================================================================
%%% analysis_plus.m
%%% Created by Sina Voshtani
%%% Created on 21/09/2021
%%% =======================================================================

%%% condition to ensure if observations exist or not. If not, then directly
%%% start the next forecast
if ~isempty(ch4_sat)
    pf_ii        = c_var ;
    
%%% covaraince modelling
    
%%% construction of a vertical correlation matrix: C_z
%%% there are three options: SOAR, FOAR, and Gaussian

    Lz = 1;
    ff = 1; %%f = 1e-5: no v-correaltion
    L_sigma=2; % Lc=1:40
    Lc =ff*(hcmaq_siglvl(1)-hcmaq_siglvl(1+L_sigma)); % Menard et al., Monthly Weather Review, 2016
    for l1  = 1:numel(hcmaq_siglvl)-1
        for l2  = 1:numel(hcmaq_siglvl)-1
            d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1)));
            C_z(l2,l1)   = (1+(d_z(l2,l1)/Lc)/1.3494) .* exp(-(d_z(l2,l1)/Lc)/1.3494); %SOAR
            %%%            %d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
            %%%            %    ((hcmaq_siglvl(1)-hcmaq_siglvl(end))); % 2-sigma correlation lenght
            %%%	        if (l2 > Lz) && (l2 < 45-Lz)
            %%%        	    d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
            %%%        	        abs(alpha*(hcmaq_siglvl(l2-Lz)-hcmaq_siglvl(l2+Lz)));
            %%%       		 elseif (l2 <= Lz)
            %%%            		d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
            %%%                	abs(2*alpha*(hcmaq_siglvl(l2)-hcmaq_siglvl(l2+Lz)));
            %%%	        elseif (l2 >= 45-Lz)
            %%%        	  	d_z(l2,l1)   =  abs((hcmaq_siglvl(l2)-hcmaq_siglvl(l1))) / ...
            %%%	                abs(2*alpha*(hcmaq_siglvl(l2-Lz)-hcmaq_siglvl(l2)));
            %%%        	end
            %%%	        C_z = (1+d_z/1.3494) .* exp(-d_z/1.3494); %SOAR
            %%%               % C_z = exp(-d_z/0.5005); %FOAR -d_z/0.5005
            %%%               % C_z = exp(-d_z.^2/2); %Gaussian
        end
    end
    
%%% construction of a horizontal correlation matrix: C_xy
%%% there are three options: SOAR, FOAR, and Gaussian    

    L = 4 * abs(xproj_hcmaq(1,1) - xproj_hcmaq(2,1)); %  X108000 km
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
    for p = 1:numel(x_sat(:))
        indx_l(:,p)    = interp1(pres_m2o_l(:,p),num_l,plv_sat(:,p),'linear','extrap');
    end
    indx_l(indx_l<1) = 1;
    indx_l(indx_l>45) = 45;
    low_l = floor(indx_l);
    low_l(low_l>44)=44; low_l(low_l<1)=1;
    indx_l(1,:) = 1; inxd_l(numel(plv_sat(:,1)),:) = 45;
    low_l(1,:)  = 1; low_l(numel(plv_sat(:,1)),:)  = 44;
    
%%% compute pressure weight on model space 
    % pwe_m  = zeros(numel(indx_l(:,1))-1,numel(x_sat));
    for p=1:numel(x_sat(:))
        for l2 = 1:numel(indx_l(:,1))-1
            pwe_m{p,l2} = (pres_m2o_l(low_l(l2,p):low_l(l2+1,p),p)-pres_m2o_l(low_l(l2,p)+1:low_l(l2+1,p)+1,p))...
                ./ (pres_m2o_l(low_l(l2,p),p)-pres_m2o_l(low_l(l2+1,p)+1,p));
        end
    end
    %%% with uniform veritcal
    for l2 = 1:numel(indx_l(:,1))
        low(l2) = mode(low_l(l2,:));
    end
    for l2 = 1:numel(indx_l(:,1))-1
        
        pwe_m1{l2} = (pres_m2o_l(low(l2):low(l2+1),p)-pres_m2o_l(low(l2)+1:low(l2+1)+1,p))...
            ./ (pres_m2o_l(low(l2),p)-pres_m2o_l(low(l2+1)+1,p));
    end
    
%%% computation of HP and HP' (outer loop)
    hp  = zeros(numel(x_sat(:)), numel(xproj_hcmaq(:))*44);
    tic
    for p=1:numel(x_sat(:))
        
        %hh = sqrt(pf_ii_hh(:,1)) .* C_xy_hh(p,:); % H^h (sigma_ii . C_xy_hh);
        C_hh_p = C_xy_hh(p,:);
        v1 =  C_z * diag(realsqrt(pf_ii_m2o(:,p)));
        
        for l2 = 1:numel(indx_l(:,1))-1
            v_o(l2,:)  =  (1*v1(:,low_l(l2,p):low_l(l2+1,p))...
                * pwe_m{p,l2}(1:end));
        end
        vv         =  pwe_sat(:,p)' * (avk_sat(:,p) .* v_o);
        C_hp       = kron(vv,C_hh_p);
        hp(p,:)    = realsqrt(pf_ii_vec) .* C_hp;
        %hp(p,:) = kron(vv,hh);
        %    p
    end
    p
    toc
    hp_3d = reshape(hp',187,187,44,numel(x_sat));
    
%%% computation of HPH' or H(HP')' (inner loop)
    tic
    hph = zeros(numel(x_sat),numel(x_sat));
    for p1 = 1:numel(x_sat)
        for l1 = 1:44
            H_hor_hp        = griddedInterpolant(xproj_hcmaq,yproj_hcmaq,squeeze(hp_3d(:,:,l1,p1)));
            hp_m2o(l1,:)  = H_hor_hp(xproj_sat(:),yproj_sat(:));
        end
        % fast approximation
        % % %     for l2 = 1:numel(indx_l(:,1))-1
        % % %         hp_o(l2,:) =  ( hp_m2o(low(l2):low(l2+1),:)'...
        % % %                        * pwe_m1{l2}(1:end));
        % % %     end
        % % %     hph_l       =  pwe_sat .* (avk_sat .* hp_o);
        % % %     hph(p1,:)   =  sum(hph_l,1);
        
        % original method
        for p = 1:numel(x_sat)
            for l2 = 1:numel(indx_l(:,1))-1
                hp_o(l2,p) = (1*hp_m2o(low_l(l2,p):low_l(l2+1,p),p)'...
                    * pwe_m{p,l2}(1:end));
            end
            hph(p1,p) =  pwe_sat(:,p)' * (avk_sat(:,p) .* hp_o(:,p));
        end
    end
    p1
    
%%% matrix inversion or Gamma = (R+HPH')^-1
    toc
    Gamma      = hph + diag((sig_sat*1e-3).^2);
    eigval = eig(Gamma);
    % Gamma_inv  = invChol_mex(Gamma);
    isposdef = all(eigval>0);
    issym    = issymmetric(Gamma);
    [~,cp] = chol(Gamma);
    if (cp == 0)
        L1_Gamma  = chol(Gamma);
        Gamma_inv = inv_chol(L1_Gamma');
        ee=1
    elseif (issym==0) && (isposdef==1)
        Gamma_inv  = inv(Gamma);
        e2=2
    else
        nspd_Gamma = nearestSPD(Gamma);
        %L2_Gamma  = chol(nspd_Gamma);
        L2_Gamma  = chol(nspd_Gamma,'lower');
        Gamma_inv = inv_chol(L2_Gamma');
        ee=3
    end
    % Gamma_inv  = inv(Gamma);
    
%%% Computation of model column averaged (before get an analysis)
    xf = c_con; % ppm
    xf_hv = reshape(xf,numel(xproj_hcmaq),44); %% 3D concentrations
    fgs_x_sat  = fgs_sat * 1e-3; % fgs is on ppb -  multiply by 1e-3
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
    
%%% Computation of the analysis
    innov         = (ch4_sat*1e-3 - xf_o_avg'); 
    chi_sq        = innov' * Gamma_inv * innov;
    hp_hv         = reshape(hp_3d,187*187,44,[]);
    tic
    %%% treat residual and sparsity
    hp_hv1 = hp_hv;
    hp_hv1(hp_hv1<1e-4*mean(hp_hv1(:))) = 0;
    % Gamma_inv_spa = sparse(Gamma_inv);
    % pf_ii_hv_spa  = sparse(pf_ii_hv);
    tic
    %%% intialization of analysis and analysis error varaince 
    xa_hv       = zeros(numel(xproj_hcmaq),44);
    pa_ii_hv    = zeros(numel(xproj_hcmaq),44);
    %%% main loop to compute the analysis and analysis error variance
    for l = 1:44
        %%% full matrix calculation for analysis
        hp_l        = squeeze(hp_hv(:,l,:));
        xa_hv_inc_l = (hp_l * Gamma_inv * innov);
        xa_hv(:,l)  = xf_hv(:,l) + xa_hv_inc_l;
        %%% matrix calculation for the error variance 
        %     hp_hv1_l         = sparse(squeeze(hp_hv1(:,l,:)));
        hp_hv1_l         = squeeze(hp_hv1(:,l,:));
        pf_ii_inc_l      = sum(hp_hv1_l*Gamma_inv.*hp_hv1_l,2); % =diag(hp_hv_spa*Gamma_inv_spa*hp_hv_spa')
        pa_ii_hv(:,l)    = abs(pf_ii_hv(:,l) - pf_ii_inc_l);
        l
    end
    toc
    
%%% convert from vector to 3D matrices for concentration and error varaince     
    pa_ii_hv = full(pa_ii_hv);
    xa_3d = reshape(xa_hv,187,187,44);
    xf_3d = c_con; % = xf = reshape(xf_hv,187,187,44);
    xd_3d = xa_3d - xf_3d;
    
    pa_3d = reshape(pa_ii_hv,187,187,44);
    pf_3d = c_var; % = pf_ii = reshape(pf_ii_hv,187,187,44);
    pd_3d = pa_3d - pf_3d;
    
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
    xf_o_avg = 0;
    hph = 0;
end

%%% store main variables for each timestep
innov_all(1:length(innov),step_size*(i-1)+nstep) = innov;
lat_sat_all(1:length(lat_sat),step_size*(i-1)+nstep) = lat_sat;
lon_sat_all(1:length(lon_sat),step_size*(i-1)+nstep) = lon_sat;
ch4_sat_all(1:length(lat_sat),step_size*(i-1)+nstep) = ch4_sat;
sig_sat_all(1:length(lat_sat),step_size*(i-1)+nstep) = sig_sat;
plv_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep) = plv_sat;
pwe_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep) = pwe_sat;
fgs_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep) = fgs_sat;
avk_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep) = avk_sat;
xf_o_avg_all(1:length(xf_o_avg),step_size*(i-1)+nstep) = xf_o_avg;
var_B    = diag(hph);
var_B_all(1:length(var_B),step_size*(i-1)+nstep) = var_B; % diagonal of hph

%%% chi2 diagnostic metric
chi_sq_all(step_size*(i-1)+nstep)       = chi_sq;
obs_num_all(step_size*(i-1)+nstep)      = numel(x_sat);
chi_sq
obs_num = obs_num_all(step_size*(i-1)+nstep)
norm_chi=chi_sq/obs_num

%%% =======================================================================
%%% END
%%% =======================================================================
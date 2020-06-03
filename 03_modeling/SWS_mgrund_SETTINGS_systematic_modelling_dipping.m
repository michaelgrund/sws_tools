function SWS_mgrund_SETTINGS_systematic_modelling_dipping()

% 1) setup dipping layer models (top)
% 2) setup 1-layer models (bottom)


%%
clear all
clc

% define incidence angle that should be used
INCvals=11;

% this setting needs ~ 12 min to be created (~485 MB disk space)
downdipdir=0:5:360;
dips=0:5:75;
thickness=0:25:500;

% dips=0:30:60;
% downdipdir=0:30:360;
% thickness=0:50:250;

% get all possible cominantions
comb_vecs=combvec(downdipdir,dips,thickness);

N = length(comb_vecs);   
MODELLstruct = repmat(struct('fast_eff',zeros(1,360), 'tlag_eff',zeros(1,360),...
    'azi',zeros(1,360), 'downdipdir',zeros(1,1), 'dip',zeros(1,1),...
    'thick',zeros(1,360), 'azi4plot',zeros(1,360), 'fast_eff4plot',zeros(1,360), 'tlag_eff4plot',zeros(1,360)), N, 1 );


start=datestr(now);

parfor ZZ=1:N

    current_model_in=comb_vecs(:,ZZ);
    
    downdipdir=current_model_in(1,:);
    dips=current_model_in(2,:);
    thickness=current_model_in(3,:);
    
    
     [fast_eff,tlag_eff,fast_eff1,tlag_eff1,azi,azi4plot]=SWS_mgrund_MSAT_split_model_4_systematic_search(INCvals,dips,downdipdir,thickness);

                MODELLstruct(ZZ).fast_eff=fast_eff;
                MODELLstruct(ZZ).tlag_eff=tlag_eff;
                MODELLstruct(ZZ).azi=azi;
                MODELLstruct(ZZ).downdipdir=downdipdir;
                MODELLstruct(ZZ).dip=dips;
                MODELLstruct(ZZ).thick=thickness;
                MODELLstruct(ZZ).azi4plot=azi4plot;
                MODELLstruct(ZZ).fast_eff4plot=fast_eff1;
                MODELLstruct(ZZ).tlag_eff4plot=tlag_eff1;
     
     
     
    

end

save(['models_out_DIPPING_layer_domper' num2str(8) 's' num2str(now) '.mat'],'MODELLstruct', '-v7.3')

endtime=datestr(now);


%% gen one layer models

stepdt=0.2;
stepphi=5;
  
fast1=-90:stepphi:90;
delaydt1=0:stepdt:4;

comb_vecs=combvec(fast1,delaydt1);
N = length(comb_vecs);  

    
for ii=1:N
    
         
    current_model_in=comb_vecs(:,ii); 

    model_phis=current_model_in(1:length(current_model_in)/2,:);
    model_dts=current_model_in(length(current_model_in)/2+1:end,:);

        model_out(ii).phi_eff=ones(1,length([0:1:360]))*model_phis;
        model_out(ii).dt_eff=ones(1,length([0:1:360]))*model_dts;
        model_out(ii).phi_in=model_phis;
        model_out(ii).dt_in=model_dts;
        model_out(ii).counter=1;

end
    
    
    
 save(['SS_models_out_layers1_dfreq' num2str(0.125) '_stepphi5_stepdt0.2_minus90_to_plus90.mat'],'model_out', '-v7.3')
    
    



%%


























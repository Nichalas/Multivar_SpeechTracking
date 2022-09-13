function hcp_15_lcmv_copula(ss,bb)
%% hcp configuration
homeDir = '/scratch/tmp/nchalas/';
tooldir = '/home/n/nchalas/fieldtrip/fieldtrip-20200513';
projDir = [homeDir 'hcp_cspeech'];

%load fieldtrip
addpath([tooldir])
addpath(genpath('home/n/nchalas/gcmi-master/'));
ft_defaults
% %%
% %fi =  find(source.avg.pow(~isnan(source.avg.pow)));
% load([homeDir '/results/dics/source0.mat']);
% %source.freq = [1:16];
% source.pow=zeros(42432,1);
% source=rmfield(source,'avg');
%% this returns the location where FieldTrip is installed

% and this is where the template source models are

 template = load(fullfile(tooldir, 'template', 'sourcemodel','standard_sourcemodel3d5mm'));
% 
% % load the atlas 
%  fname_in = fullfile(tooldir,'template','atlas','atlas_MMP1.0_4k.mat');
%  atlas = ft_read_atlas(fname_in);
% % 
%  cfg = [];
%  cfg.interpmethod = 'nearest';
%  cfg.parameter = 'parcellation';
%  atemplate = ft_sourceinterpolate(cfg, atlas, template.sourcemodel);
load([projDir '/results/atemplate.mat']);

% fname_in = fullfile(tooldir,'template','atlas','melb_subcortical','melb_sub.mat');
% melb_sub = ft_read_atlas(fname_in);
% 
% cfg                         = [];
% cfg.interpmethod            = 'nearest';
% cfg.parameter               = 'tissue';
% sub_atemplate     = ft_sourceinterpolate(cfg, melb_sub,template.sourcemodel);
% bM = logical(sub_atemplate.tissue);
% sub_atemplate.tissue = (sub_atemplate.tissue + length(atemplate.parcellationlabel)).*bM;
% 
% mixed_labels=atemplate.parcellationlabel; 
% for i=1:length(sub_atemplate.tissuelabel)
%     mixed_labels{length(atemplate.parcellationlabel)+i}=sub_atemplate.tissuelabel{i};
% end
% 
% mixed_atemplate = atemplate;
% mixed_atemplate.parcellation = atemplate.parcellation + sub_atemplate.tissue;
% mixed_atemplate.parcellationlabel = mixed_labels; 
%% 
%number of participants
pps = {'01','03','05','06','07','08','09','10','11','13','15','17','19','21','22','23','24','26','27','28','29','30','31','33'};
%load the speech envelope
load([projDir '/audio/cEnvNogap'])
cEnv_run = reshape(cEnv, 55000, 6); % 6 runs

nBlocks = 7;
delay_samples = [-30:30]; %old delay_samples = 1:30 
dels = [1:numel(delay_samples)];
delays = delay_samples*10;
Ndel = length(delay_samples);
%d=5; %fixed 50ms delay
ncomps=3;
nsurr=500;
%for bb = 2:nBlocks
fprintf('run %d \n', bb); %
load([homeDir 'hcp_cspeech/hcpmmp1/lcmv_hcpmmp1_ori_pca_ss' pps{ss} '_run' num2str(bb) '.mat'],'pca_label_ts')

allmi=[];
Env=cEnv_run(:,bb-1);
% Env=circshift(Env,d,1);
dEnv=[0 diff(Env)']';
fb = cwtfilterbank('SignalLength',numel(Env),'SamplingFrequency',100,'FrequencyLimits',[0.5 40],'WaveletParameters',[3,10]); %3 60
[Envm,fw]=wt(fb,Env);
Envm=Envm./abs(Envm);
Envm(isnan(Envm))=0;
[dEnvm,fw]=wt(fb,dEnv);
dEnvm=dEnvm./abs(dEnvm);
dEnvm(isnan(dEnvm))=0;

% stgl_idx = strmatch('L_STGa_ROI', atemplate.parcellationlabel);
% stgr_idx = strmatch('R_STGa_ROI', atemplate.parcellationlabel);
% 
% lStg=pca_label_ts{stgl_idx}(1,:);
% rStg=pca_label_ts{stgr_idx}(1,:);
% [lStgc,fw]=wt(fb,lStg);
% lStgc=lStgc./abs(lStgc);
% [rStgc,fw]=wt(fb,rStg);
% rStgc=rStgc./abs(rStgc);
% 
% z1=copnorm(real(lStgc)') ;
% z2= copnorm(imag(lStgc)');
% z3=copnorm(real(rStgc)') ;
% z4= copnorm(imag(rStgc)');
for l=1:length(atemplate.parcellationlabel)
    fprintf('parcel %d \n', l); %
    meg=pca_label_ts{l,1};
    
    if size(meg,1)==1
        
        [megc1,fw]=wt(fb,meg(1,:));
        megc1=megc1./abs(megc1);
        mi=zeros(numel(delay_samples),length(fw));
        
        for i=1:numel(delay_samples)
            
            de=delay_samples(i);
            d=dels(i);
            Envc = circshift(Envm,de,2); %delay of the envelope
            dEnvc = circshift(dEnvm,de,2);
            
            x1=copnorm(real(Envc)') ;
            x2= copnorm(imag(Envc)');
            x3=copnorm(real(dEnvc)') ;
            x4= copnorm(imag(dEnvc)');
            y1=copnorm(real(megc1)') ;
            y2= copnorm(imag(megc1)');
            
            nf=size(Envm,1);
            
            %                 for k=1:nf
            %                     mi(d,k)=mi_gg([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],[y1(:,k) y2(:,k)]);
            %                 end
            for k=1:nf
                mi(d,k)=mi_gg([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],[y1(:,k) y2(:,k)]);
            end
        end
        
        tmp=zeros(nsurr,length(fw));
        nfs=size(Envm,2);
        
        for is=1:nsurr
            %ri=randperm(nsamp);
            sh=randi([round(0.05*nfs) nfs-round(0.05*nfs)]);
            for k=1:nf
                tmp(is,k)=mi_gg(circshift([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],sh,1),[y1(:,k) y2(:,k)]);
            end
        end
        
        sm=mean(tmp,1);sd=std(tmp,0,1);
        
        miz(:,:)=(mi(:,:)-sm)./sd;
        
        allmi{l}=miz;
    elseif size(meg,1)==2
        [megc1,fw]=wt(fb,meg(1,:));
        [megc2,fw]=wt(fb,meg(2,:));
        
        megc1=megc1./abs(megc1);
        megc2=megc2./abs(megc2);
        mi=zeros(numel(delay_samples),length(fw));
        
        for i=1:numel(delay_samples)
            de=delay_samples(i);
            d=dels(i);
            
            Envc = circshift(Envm,de,2); %delay of the envelope
            dEnvc = circshift(dEnvm,de,2);
            
            x1=copnorm(real(Envc)') ;
            x2= copnorm(imag(Envc)');
            x3=copnorm(real(dEnvc)') ;
            x4= copnorm(imag(dEnvc)');
            y1=copnorm(real(megc1)') ;
            y2= copnorm(imag(megc1)');
            y3=copnorm(real(megc2)') ;
            y4= copnorm(imag(megc2)');
            
            nf=size(Envm,1);
            
            for k=1:nf
                mi(d,k)=mi_gg([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],[y1(:,k) y2(:,k)]);
            end
        end
        
        tmp=zeros(nsurr,length(fw));
        nfs=size(Envm,2);
        
        for is=1:nsurr
            %ri=randperm(nsamp);
            sh=randi([round(0.05*nfs) nfs-round(0.05*nfs)]);
            for k=1:nf
                tmp(is,k)=mi_gg(circshift([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],sh,1),[y1(:,k) y2(:,k) y3(:,k) y4(:,k)]);
            end
        end
        
        sm=mean(tmp,1);sd=std(tmp,0,1);
        miz(:,:)=(mi(:,:)-sm)./sd;
        
        allmi{l}=miz;
    else
        [megc1,fw]=wt(fb,meg(1,:));
        [megc2,fw]=wt(fb,meg(2,:));
        [megc3,fw]=wt(fb,meg(3,:));
        
        megc1=megc1./abs(megc1);
        megc2=megc2./abs(megc2);
        megc3=megc3./abs(megc3);
        mi=zeros(numel(delay_samples),length(fw));
        
        for i=1:numel(delay_samples)
            de=delay_samples(i);
            d=dels(i);
            Envc = circshift(Envm,de,2); %delay of the envelope
            dEnvc = circshift(dEnvm,de,2);
            
            x1=copnorm(real(Envc)') ;
            x2= copnorm(imag(Envc)');
            x3=copnorm(real(dEnvc)') ;
            x4= copnorm(imag(dEnvc)');
            y1=copnorm(real(megc1)') ;
            y2= copnorm(imag(megc1)');
            y3=copnorm(real(megc2)') ;
            y4= copnorm(imag(megc2)');
            y5=copnorm(real(megc3)') ;
            y6= copnorm(imag(megc3)');
            
            nf=size(Envm,1);
            
            for k=1:nf
                mi(d,k)=mi_gg([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],[y1(:,k) y2(:,k) y3(:,k) y4(:,k) y5(:,k) y6(:,k)]);
            end
        end
        tmp=zeros(nsurr,length(fw));
        nfs=size(Envc,2);
        
        for is=1:nsurr
            %ri=randperm(nsamp);
            sh=randi([round(0.05*nfs) nfs-round(0.05*nfs)]);
            for k=1:nf
                tmp(is,k)=mi_gg(circshift([x1(:,k) x2(:,k) x3(:,k) x4(:,k)],sh,1),[y1(:,k) y2(:,k) y3(:,k) y4(:,k) y5(:,k) y6(:,k)]);
            end
        end
        
        sm=mean(tmp,1);sd=std(tmp,0,1);
        miz(:,:)=(mi(:,:)-sm)./sd;
        
        allmi{l}=miz;
    end
end
save([projDir '/results/lcmv_hcpmmp1_copula_negdel_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
%end
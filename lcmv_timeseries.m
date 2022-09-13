function hcp_10_lcmv(ss)
%% hcp configuration
homeDir = '/scratch/tmp/nchalas/';
tooldir = '/home/n/nchalas/fieldtrip/fieldtrip-20200513';
projDir = [homeDir 'hcp_cspeech'];

%load fieldtrip
addpath([tooldir])
ft_defaults
% %%
% %fi =  find(source.avg.pow(~isnan(source.avg.pow)));
% load([homeDir '/results/dics/source0.mat']);
% %source.freq = [1:16];
% source.pow=zeros(42432,1);
% source=rmfield(source,'avg');
%% this returns the location where FieldTrip is installed

% and this is where the template source models are

%templatedir = fullfile(tooldir, 'external', 'spm8', 'templates');
%template_mri = ft_read_mri(fullfile(templatedir, 'T1.nii'));

template = load(fullfile(tooldir, 'template', 'sourcemodel','standard_sourcemodel3d5mm'));

% load the atlas 
%fname_in = fullfile(tooldir,'template','atlas','aal','ROI_MNI_V4.nii');
fname_in = fullfile(tooldir,'template','atlas','atlas_MMP1.0_4k.mat');
atlas = ft_read_atlas(fname_in);

% and call ft_sourceinterpolate:
% cfg = [];
% cfg.interpmethod = 'nearest';
% cfg.parameter = 'parcellation';
% atemplate = ft_sourceinterpolate(cfg, atlas, template.sourcemodel);
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

% and call ft_sourceinterpolate:
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'parcellation';
atemplate = ft_sourceinterpolate(cfg, atlas, template.sourcemodel);

%% extract pca time series within label from 3 orientations 
%number of participants
pps = {'01','03','05','06','07','08','09','10','11','13','15','17','19','21','22','23','24','26','27','28','29','30','31','33'};
%load the speech envelope
load([projDir '/audio/cEnvNogap'])
cEnv_run = reshape(cEnv, 55000, 6); % 6 runs

nBlocks = 7;
delay_samples = 1:30;
delays = delay_samples*10;
Ndel = length(delay_samples);
ff=[0.5:0.5:8];
nfoi = length(ff);


load([projDir '/vol/vol_ss',pps{ss},'.mat'],'vol')
for bb = 2:nBlocks
    fprintf('run %d \n', bb);
    
    load([projDir '/leadfields/leGrid_ss',pps{ss},'_bb',num2str(bb),'_5mm.mat'],'leGrid')
    load([projDir '/data/data_ss' pps{ss} '_run' num2str(bb) '.mat'],'data');
    
%     cfg = [];
%     cfg.resamplefs = 200; % was 125 for ICA
%     cfg.detrend = 'no';
%     cfg.demean = 'no';

%     data = ft_resampledata(cfg,data);
    
    grids = find(leGrid.inside);

    cfg = [];
    cfg.covariance = 'yes';
    data = ft_timelockanalysis(cfg,data);

    pca_label_ts={};%zeros(length(atemplate.tissuelabel),55000);
    %pca_label_ori_ts=zeros(3,length(atemplate.tissuelabel),55000);

    %all_ts = lcmv.avg.mom(fi);

    for l=1:length(atemplate.parcellationlabel)        
        tmpindx = find(atemplate.parcellation(grids)==l);

        cfg = [];
        cfg.covariance = 'yes';
        data = ft_timelockanalysis(cfg,data);

        %constrain the source analysis within the label 
        points = grids(tmpindx);
        leGrid.inside = false(size(leGrid.inside,1),1);
        leGrid.inside(points)=true;

        % build spatial filter, project data in source space
        cfg = [];
        cfg.method = 'lcmv';
        cfg.grid = leGrid;
        cfg.headmodel = vol;
        cfg.sourcemodel.pos = leGrid;
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.projectnoise = 'yes';
        cfg.lcmv.fixedori = 'no';
        cfg.lcmv.projectmom = 'no';
      % cfg.lcmv.lambda = '5%'; %'0%' when used with oasCovariance
        lcmv = ft_sourceanalysis(cfg,data);

        tmp=cat(1,lcmv.avg.mom{points})*1e12;

                         % #1 pca for orientatiotion

        [coeff,score,latent,tsquared,explained,mu] = pca(tmp');

        tmp_var=0;
        for i=1:length(explained)
            tmp_var=tmp_var+explained(i);
            if tmp_var>90    % extract >99% of variance        
                sig_pc=score(:,1:i); 
                break
            else                    
                continue                   
            end
        end

        pca_label_ts{l,1}=sig_pc';
        pca_label_ts{l,2}=explained(1:size(sig_pc,2)); %insert to the container
    end
    save([projDir '/results/lcmv200_hcpmmp1_ori_pca_ss' pps{ss} '_run' num2str(bb) '.mat'],'pca_label_ts');
end



















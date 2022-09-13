homeDir = '/remotedata/AgGross/cspeech/';
tooldir = '/usr/local1/ctf';
projDir = [homeDir 'hcp_cspeech'];
sourceDir = [homeDir 'data/MEG/'];

cd '/remotedata/AgGross/cspeech/tmp_figures'
addpath(genpath('/home/chalas/Documents/MATLAB/meth_2019_04/meth/'));
addpath('/home/chalas/Documents/MATLAB/BrewerMap/')
addpath('/home/chalas/Documents/MATLAB/Colormaps/')
addpath('/home/chalas/Documents/MATLAB/BCT/2019_03_03_BCT')
addpath('/home/chalas/Documents/MATLAB/circularGraph/')
addpath('/remotedata/AgGross/Toolboxes/');
addpath('/remotedata/AgGross/cspeech/functions/');
addpath('~/Documents/Codebase_VisualGamma_Elife/External/');

get(groot, 'defaultAxesFontName')
set(groot, 'defaultAxesFontName', 'Arial')

% load the atlas 
fname_in = fullfile('/home','chalas','Documents','MATLAB','my_fieldtrip','fieldtrip-20200513','template','atlas','atlas_MMP1.0_4k.mat');
atlas_8k = ft_read_atlas(fname_in);

atlas = ft_read_cifti('/home/chalas/Documents/MATLAB/my_fieldtrip/fieldtrip-20200513/template/atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

%load the sourcemodel
filename = '/home/chalas/Documents/MATLAB/fieldtrip/template/sourcemodel/S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii';
sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

filename = '/home/chalas/Documents/MATLAB/fieldtrip/template/sourcemodel/S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii';
sourcemodel_left = ft_read_headshape({filename, strrep(filename, '.R.', '.L.')});
% weird but otherwise we face problems
sourcemodel_left.pos(length(sourcemodel.pos)/2:end,:)=NaN; 
sourcemodel_left.tri(length(sourcemodel.pos):end,:)=NaN;

filename = '/home/chalas/Documents/MATLAB/fieldtrip/template/sourcemodel/S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii';
sourcemodel_right=ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

%% study configuration
%number of participants
pps = {'01','03','05','06','07','08','09','10','11','13','15','17','19','21','22','23','24','26','27','28','29','30','31','33'};
%load the speech envelope
load([projDir '/audio/cEnvNogap'])
cEnv_run = reshape(cEnv, 55000, 6); % 6 runs

nBlocks = 7;
delay_samples = 1:30;
%delay_samples = [1 10 20 30];
delays = delay_samples*10;
Ndel = length(delay_samples);
nlab=362;
load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/fw.mat');
tmp = fw<21;
cfw = nonzeros(tmp.*fw);
ffw= flip(fw);

%% figure #1
% gather the data in one matrix and average (we are using only one cond)
block_dat = zeros(2,6,nlab,Ndel,length(fw)); % conditions | runs | parcels | delays | frequencies
in_dat = zeros(2,24,nlab,Ndel,length(fw)); % conditions | subjs | parcels | delays | frequencies

for ss=1:numel(pps)
    fprintf('subject %d \n', ss);
    for bb=2:7
            load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_denv_tpca_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
             for ilab=1:nlab
                block_dat(1,bb-1,ilab,:,:)=allmi{ilab};
            end
            load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
            for ilab=1:nlab
                block_dat(2,bb-1,ilab,:,:)=allmi{ilab};
        end
in_dat(:,ss,:,:,:) = mean(block_dat,2);
    end
end

% cluster permutation stats
freq=[];
freq.label=atlas_8k.parcellationlabel;
freq.freq=[1:64];
freq.powspctrm=[];
freq.dimord='chan_freq';

dat1={};
dat2={};
for ss=1:24
    dat1{ss}=freq;
    dat1{ss}.powspctrm=squeeze(max(in_dat(2,ss,:,:,:),[],4)); 
    dat2{ss}=freq;
    dat2{ss}.label=atlas_8k.parcellationlabel;
    dat2{ss}.powspctrm=1.6*ones(362,64);%compare with 95th percentile
end
nsubj=24;
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster'; %cluster
cfg.parameter = 'powspctrm';
cfg.alpha            = 0.05;
cfg.clusterstatistic = 'wcm';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 0;          % minimum number of neighborhood frechannels that is
                                   % required for a selected sample to be included   % in the clustering algorithm (default=0)
cfg.tail             =  1;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      =  1;
cfg.neighbours=[];
cfg.numrandomization =5000;  %set this to 5000 for final statistics (correcting for multiple comparisons)
% set up design matrix
cfg.design  = [1:nsubj 1:nsubj; ones(1, nsubj) 2*ones(1, nsubj)];
cfg.uvar    = 1;
cfg.ivar    = 2;
% run t-test
stat = ft_freqstatistics(cfg, dat1{:}, dat2{:});

tstat=stat.mask.*stat.stat;
tstat(tstat==0) = NaN;

% find the total number of significant clusters
p_clusters = [stat(:).posclusters.prob];
n_sig_cluster = length(find(p_clusters<0.05));

% XXXXXXXX Spectra
% find left and right significant parcels
fi =[];
for f=1:length(fw)
    tmp = find(stat.mask(:,f));
    fi = [fi; tmp];
end
fi = unique(fi);
fi = fi';

lh=[];
rh=[];
for l=1:length(fi)
    ilab=fi(l);
    tmp_tstat=nansum(tstat(ilab,:))/length(fw);    
    tmp_roiidx=find(atlas_8k.parcellation==ilab);
    tmp_pos=atlas_8k.pos(tmp_roiidx,:);
    if tmp_pos(1,1)<0
        lh=[lh ilab];
    else 
        rh=[rh ilab];
    end
end

cmap = brewermap(nlab,'PuBuGn');

% plot spctrm for the masked seperately for lh and rh
labs={lh,rh}; 
% for the first conditionvas
for i=1:2 % left and right hemisphere
    ilab = labs{i};
    y = squeeze(mean(max(in_dat(1,:,ilab,:,1:length(fw)),[],3),2));
    tmp = squeeze(max(in_dat(1,:,ilab,:,1:length(fw)),[],3));
    subplot(2,2,i)
    h=shadedErrorBar(fw,mean(y,1),std(y),{'black','LineWidth',3});
    xlim([fw(end) 20]);
    ylim([0 27]);   
    ax = gca;
    set(gca,'box','off')
    set(gca,'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
    if i==1
        ylabel('MI')
        title('Left','FontWeight','normal');
    else
        title('Right','FontWeight','normal');
    end
    ax.FontSize = 18;
    subplot(2,2,i+2)
    plot(fw,tstat(ilab,1:length(fw)),'lineWidth',3,'color',[cmap(300,:) 0.05]); %change to opacity of the lines with the 4th number
    hold on
    plot(fw,nanmean(tstat(ilab,1:length(fw)),1),'lineWidth',3,'color',[cmap(300,:) 1]);
    ax = gca;
    set(gca,'box','off')
    set(gca,'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
    xlabel('Frequency (Hz)') 
    ylim([min(min(tstat(ilab,1:length(fw)))) 11.5]);
    xlim([fw(end) 20]);
    hold on;
    if i==1
        ylabel('t-values')
    else              
    end    
    ax.FontSize = 17;
    hold off;   
end
print('95th_tstat_mi','-djpeg','-r300');

cmap = brewermap(200,'BuPu');
% for the first conditionvas
d=10;
y = squeeze(trimmean(in_dat(1,:,fi,d,1:length(fw)),95,2));
colormap(cmap)
imagesc(fw,1:length(ilab),y)
xlim([fw(end) 20])
ax=gca;
ax.FontSize = 17;
xlabel('Frequency (Hz)')
ylabel('Parcel #')
set(ax,'yticklabel',[])
cb=colorbar;
caxis([0 25])
cb.FontSize = 20;

print('parcels_freq','-djpeg','-r300');

% XXXXXXXXXXXXXXXXXXX Surface plots
atlas.left = zeros(1,length(atlas.indexmax));
atlas.right = zeros(1,length(atlas.indexmax));
for ilab=1:nlab
    str_label = atlas_8k.parcellationlabel{ilab};
    index = find(contains(atlas.indexmaxlabel,str_label));
    tmp_tstat=nansum(tstat(ilab,:))/length(fw);
    try
    tmp_roiidx=find(atlas.indexmax==index);
    tmp_pos = atlas.brainstructure(tmp_roiidx);
    if tmp_pos(1) == 1

        atlas.left(tmp_roiidx)=tmp_tstat;
        atlas.right(tmp_roiidx)=NaN;
    else 
        atlas.left(tmp_roiidx)=NaN;
        atlas.right(tmp_roiidx)=tmp_tstat;
    end
    end
end

tmp = brewermap(362,'PuBuGn');
cmap = [0.9 0.9 0.9; tmp];
hcp_allviews(atlas,sourcemodel_left,sourcemodel_right,cmap,'t-values')

print('95th_surface_32k','-djpeg','-r300');
%% Figure #2 (3vs1 components)
% gather the data in one matrix and average (we are using only one cond)
block_dat = zeros(2,5,nlab,Ndel,length(fw)); % conditions | runs | parcels | delays | frequencies
in_dat = zeros(2,24,nlab,Ndel,length(fw)); % conditions | subjs | parcels | delays | frequencies
 
for ss=1:numel(pps)
    fprintf('subject %d \n', ss);
    for bb=2:7
        if bb==7
            continue
        else
            load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
            for ilab=1:nlab
                block_dat(1,bb-1,ilab,:,:)=allmi{ilab};
            end
            load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_fpca_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
            for ilab=1:nlab
                block_dat(2,bb-1,ilab,:,:)=allmi{ilab};
            end
        end
    end
in_dat(:,ss,:,:,:) = mean(block_dat,2);
end

% cluster permutation stats
freq=[];
freq.label=atlas_8k.parcellationlabel;
freq.freq=[1:64];
freq.powspctrm=[];
freq.dimord='chan_freq';

for ss=1:24
    dat1{ss}=freq;
    dat1{ss}.powspctrm=squeeze(max(in_dat(1,ss,:,:,:),[],4)); 
    %dat1{ss}.powspctrm=squeeze((in_dat3(ss,:,[1 10 20 30],:)));
    dat2{ss}=freq;
 %   dat2{ss}.label=atlas.indexmaxlabel;
    dat2{ss}.powspctrm=squeeze(max(in_dat(2,ss,:,:,:),[],4));%compare with 95th percentile
end
nsubj=24;
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
%            cfg.statistic        = 'ft_statfun_pooledT';
cfg.correctm         = 'cluster'; %cluster
cfg.parameter = 'powspctrm';
cfg.alpha            = 0.05;
cfg.clusterstatistic = 'wcm';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 0;          % minimum number of neighborhood frechannels that is
                                   % required for a selected sample to be included   % in the clustering algorithm (default=0)
cfg.tail             =  1;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      =  1;
cfg.neighbours=[];
cfg.numrandomization =5000;  %set this to 5000 for final statistics (correcting for multiple comparisons)
cfg.design  = [1:nsubj 1:nsubj; ones(1, nsubj) 2*ones(1, nsubj)];
cfg.uvar    = 1;
cfg.ivar    = 2;
% run t-test
stat = ft_freqstatistics(cfg, dat1{:}, dat2{:});

tstat=stat.mask.*stat.stat;
tstat(tstat==0) = NaN;

% find the total number of significant clusters
p_clusters = [stat(:).posclusters.prob];
n_sig_cluster = length(find(p_clusters<0.05));

fi =[];
for f=1:length(fw)
    tmp = find(stat.mask(:,f));
    fi = [fi; tmp];
end
fi = unique(fi);
fi = fi';

% XXXXXXXX Spectra
% find left and right significant parcels
lh=[];
rh=[];
for l=1:length(fi)
    ilab=fi(l);
    tmp_tstat=nansum(tstat(ilab,:))/length(fw);    
    tmp_roiidx=find(atlas_8k.parcellation==ilab);
    tmp_pos=atlas_8k.pos(tmp_roiidx,:);
    if tmp_pos(1,1)<0
        lh=[lh ilab];
    else 
        rh=[rh ilab];
    end
end

cmap=magma(nlab);
labs={lh,rh};  
for i=1:2 % left and right hemisphere
    ilab = labs{i};
    y = squeeze(mean(max(in_dat(1,:,ilab,:,1:length(fw)),[],3),2));
    subplot(2,2,i)
    tmpcols={};
    tmp=cmap(180,:);
    tmpcols{1} = tmp;
    y_mean = mean(y,1); y_std = std(y,[],1)/sqrt(size(y,1)).*1.96;
    lineProps.col = tmpcols;
    h=mseb(fw',y_mean,y_std,lineProps,1);
    hold on
    y = squeeze(mean(max(in_dat(2,:,ilab,:,1:length(fw)),[],3),2));
    tmpcols={};
    tmp=brewermap(6,'Greys');
    tmpcols{1} = tmp(6,:);
    y_mean = mean(y,1); y_std = std(y,[],1)/sqrt(size(y,1)).*1.96;
    lineProps.col = tmpcols;
    h=mseb(fw',y_mean,y_std,lineProps,1);
    ax = gca;
    set(gca,'box','off')
    set(gca,'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
    if i==1
        ylabel('MI')
        title('Left','FontWeight','normal');
    else
        title('Right','FontWeight','normal');
    end 
    ax.FontSize = 18;
    lgd = legend({'multivariate','univariate'},'FontSize',12,'TextColor','black','Location','north');
    lgd.NumColumns = 1;
    legend boxoff 
    xlim([fw(end) fw(1)/2]);
    hold off
    subplot(2,2,i+2)
    plot(fw,tstat(ilab,1:length(fw)),'lineWidth',3,'color',[cmap(180,:) 0.05]); %change to opacity of the lines with the 4th number
    hold on
    plot(fw,nanmean(tstat(ilab,1:length(fw)),1),'lineWidth',3,'color',[cmap(180,:) 1]);
    ax = gca;
    set(gca,'box','off')
    set(gca,'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
    %ylim([min(min(tstat(ilab,1:length(fw)))) max(max(stat.stat(ilab,1:length(fw))))]);
    ylim([min(min(tstat(ilab,1:length(fw)))) 10]);

    xlim([fw(end) fw(1)/2]);
    xlabel('Frequency (Hz)') 
    hold on;
    if i==1
        ylabel('t-values')
    else   
    end
    ax.FontSize = 17;
    hold off;   
end
print('3vs1_tstat_mi','-djpeg','-r300');

% XXXXXXXXXXXXXXXXXXX Surface plots
% k-means

% find the significant areas
sig = [lh rh];

% do the kmeans clustering 
k=3;
[c1,c2]=kmeans(stat.stat(sig,:),k);
% 
cmap=[];
tmp = magma(k+2);

cmap = [0.9686 0.9843 1; tmp([1 2 3],:)];

atlas.left = zeros(1,length(atlas.indexmax));
atlas.right = zeros(1,length(atlas.indexmax));
for l=1:numel(sig)
    str_label = atlas_8k.parcellationlabel{ilab};
    index = find(contains(atlas.indexmaxlabel,str_label));
    ilab=sig(l);
    tmp_tstat=c1(l);
    try
        tmp_roiidx=find(atlas.indexmax==index);
        tmp_pos = atlas.brainstructure(tmp_roiidx);
        if tmp_pos(1)==1
            atlas.left(tmp_roiidx)=tmp_tstat;
            atlas.right(tmp_roiidx)=NaN;
        else
            atlas.left(tmp_roiidx)=NaN;
            atlas.right(tmp_roiidx)=tmp_tstat;

        end
    end
end
hcp_allviews(atlas,sourcemodel_left,sourcemodel_right,cmap,'t-values')

print('3vs1_surface_kmeans_32k','-djpeg','-r300');

%cmap = brewermap(nlab,'PuBuGn');
cmap(round(nlab/k)*i-2)
figure;
for i=1:k
    plot(fw, c2(i,:),'lineWidth',7,'Color',cmap(i+1,:));
    xlim([fw(end) fw(1)/2]);
    lgd = legend({'1','2','3'},'FontSize',26,'TextColor','black','Location','northoutside');
    lgd.NumColumns = 4;
    legend boxoff 
    %set(gca,'box','off')
    hold on
    ax = gca;
    set(gca,'box','off')
    xlabel('Frequency (Hz)') 
    ylabel('t-values') 
    ax.FontSize = 28;
end
print('kmean_spectra','-djpeg','-r300');

%% Figure #3 (env vs -only- derivat) 
% gather the data in one matrix and average 
block_dat = zeros(2,5,nlab,Ndel,length(fw)); % conditions | runs | parcels | delays | frequencies
in_dat = zeros(2,24,nlab,Ndel,length(fw)); % conditions | subjs | parcels | delays | frequencies

for ss=1:numel(pps)
    fprintf('subject %d \n', ss);
    for bb=2:7
        if bb==7
            continue
        else
            load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
            for ilab=1:nlab
                block_dat(1,bb-1,ilab,:,:)=allmi{ilab};
            end
            load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_odenv_tpca_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
            for ilab=1:nlab
                block_dat(2,bb-1,ilab,:,:)=allmi{ilab};
            end
        end
    end
in_dat(:,ss,:,:,:) = mean(block_dat,2);
end

% cluster permutation stats
freq=[];
freq.label=atlas_8k.parcellationlabel;
freq.freq=[1:64];
freq.powspctrm=[];
freq.dimord='chan_freq';

for ss=1:24
    dat1{ss}=freq;
    dat1{ss}.powspctrm=squeeze(max(in_dat(1,ss,:,:,:),[],4)); 
    %dat1{ss}.powspctrm=squeeze((in_dat3(ss,:,[1 10 20 30],:)));
    dat2{ss}=freq;
    dat2{ss}.label=atlas_8k.parcellationlabel;
    dat2{ss}.powspctrm=squeeze(max(in_dat(2,ss,:,:,:),[],4));
end
nsubj=24;
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
%            cfg.statistic        = 'ft_statfun_pooledT';
cfg.correctm         = 'cluster'; %cluster
cfg.parameter = 'powspctrm';
cfg.alpha            = 0.05;
cfg.clusterstatistic = 'wcm';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 0;          % minimum number of neighborhood frechannels that is
                                   % required for a selected sample to be included   % in the clustering algorithm (default=0)
cfg.tail             =  0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      =  0;
cfg.neighbours=[];
%cfg.neighbours=neighbours;
cfg.numrandomization =5000;  %set this to 5000 for final statistics (correcting for multiple comparisons)
% set up design matrix
cfg.design  = [1:nsubj 1:nsubj; ones(1, nsubj) 2*ones(1, nsubj)];
cfg.uvar    = 1;
cfg.ivar    = 2;
% run t-test
stat = ft_freqstatistics(cfg, dat2{:}, dat1{:}); % envelope > only derivative

tstat=stat.mask.*stat.stat;
tstat(tstat==0) = NaN;

% find the total number of significant clusters
p_clusters = [stat(:).posclusters.prob];
n_sig_cluster = length(find(p_clusters<0.05));

% XXXXXXXX Spectra
% find left and right significant parcels
fi =[];
for f=1:length(fw)
    tmp = find(stat.mask(:,f));
    fi = [fi; tmp];
end
fi = unique(fi);
fi = fi';

lh=[];
rh=[];
for l=1:length(fi)
    ilab=fi(l);
    tmp_tstat=sum(tstat(ilab,:))/length(fw);  
    tmp_roiidx=find(atlas_8k.parcellation==ilab);
    tmp_pos=atlas_8k.pos(tmp_roiidx,:);
    if tmp_pos(1,1)<0
        lh=[lh ilab];
    else 
        rh=[rh ilab];
    end
end


%%
cmap = brewermap(nlab,'PuBuGn');
labs={lh,rh};  
for i=1:2 % left and right hemisphere
    ilab = labs{i};
    y = squeeze(mean(max(in_dat(:,:,ilab,:,1:length(fw)),[],3),2));
    subplot(2,2,i)
    tmpcols={};
    tmp=cmap(300,:);
    tmpcols{1} = tmp;
    y_mean = mean(y,1); y_std = std(y,[],1)/sqrt(size(y,1)).*1.96;
    lineProps.col = tmpcols;
    h=mseb(fw',y_mean,y_std,lineProps,1);
    hold on
    y = squeeze(mean(max(in_dat(2,:,ilab,:,1:length(fw)),[],3),2));
    tmpcols={};
    tmp=brewermap(6,'Greys');
    tmpcols{1} = tmp(6,:);
    y_mean = mean(y,1); y_std = std(y,[],1)/sqrt(size(y,1)).*1.96;
    lineProps.col = tmpcols;
    h=mseb(fw',y_mean,y_std,lineProps,1);
    ax = gca;
    set(gca,'box','off')
    if i==1
        ylabel('MI') 
        title('Left','FontWeight','normal');
    else   
        title('Right','FontWeight','normal');
    end 
    ax.FontSize = 18;
    lgd = legend({'Env','EnvD'},'FontSize',12,'TextColor','black','Location','north');
    lgd.NumColumns = 1;
    legend boxoff 
    xlim([fw(end) fw(1)/2]);
    hold off
    subplot(2,2,i+2)
    plot(fw,tstat(ilab,1:length(fw)),'lineWidth',3,'color',[cmap(300,:) 0.05]); %change to opacity of the lines with the 4th number
    hold on
    plot(fw,nanmean(tstat(ilab,1:length(fw)),1),'lineWidth',3,'color',[cmap(300,:) 1]);
    ax = gca;
    set(gca,'box','off')
    ylim([min(min(tstat(ilab,1:length(fw)))) max(max(stat.stat(ilab,1:length(fw))))]);
    xlim([fw(end) fw(1)/2]);
    xlabel('Frequency (Hz)') 
    hold on;
    if i==1
        ylabel('t-values')     
    end
    ax.FontSize = 17;
    hold off;   
end
print('envVsder_spectra','-djpeg','-r300');


% XXXXXXXXXXXXXXXXXXX Surface plots
av_tstat = nanmean(tstat,1); 
peaks = [62,21];% for f=1.7 Hz & f = 10 Hz 
cmap=[];
tmp = brewermap(nlab-1,'PuBuGn');
cmap = [0.9686 0.9843 1; tmp];

for i = 1:1 % for paper only the first peak (derivative > envelope)
    figure;
    p=peaks(i);
    atlas.left = zeros(1,length(atlas.indexmax));
    atlas.right = zeros(1,length(atlas.indexmax));
    for ilab=1:nlab
        str_label = atlas_8k.parcellationlabel{ilab};
        index = find(contains(atlas.indexmaxlabel,str_label));
        tmp_tstat=nansum(tstat(ilab,p));
        try
        tmp_roiidx=find(atlas.indexmax==index);  
        end
        tmp_pos = atlas.brainstructure(tmp_roiidx);
        if tmp_pos(1) == 1           
            atlas.left(tmp_roiidx)=tmp_tstat;
            atlas.right(tmp_roiidx)=NaN;
        else
            atlas.left(tmp_roiidx)=NaN;
            atlas.right(tmp_roiidx)=tmp_tstat;
        end
end

    hcp_allviews(atlas,sourcemodel_left,sourcemodel_right,cmap,'t-values')
    end
    print(sprintf('envder_surface'),'-djpeg','-r300');
end
%% envelope / acoustic landmarks
cmap=inferno(362);

dEnv_run = diff(cEnv_run);
for i=5000:500:10000
figure;
i
plot(zscore(dEnv_run(i:i+100,2)));hold on;plot(zscore(cEnv_run(i:i+100,2)));
end

plot(zscore(dEnv_run(10000:10100,2)),'lineWidth',3,'color',cmap(120,:));hold on;plot(zscore(cEnv_run(10000:10100,2)),'lineWidth',3,'color',[0.1451    0.1451    0.1451]);
lgd = legend({'EnvD','Env'},'FontSize',14,'TextColor','black','Location','northOutside');
lgd.NumColumns = 1;
legend boxoff 
set(gca,'XColor', 'none','YColor','none')
possub  = get(lgd,'Position');
lgd.Position = [possub(1)/2 possub(2)/2 possub(3) possub(4)];
print('envelope_envD','-djpeg','-r300');


tmp_col=brewermap(6,'Blues');
%  tmpcols={};
%  for d=1:2
%     tmp=brewermap(6,'Greens');
%     tmpcols{d} = tmp(2+d,:);
%  end
%   y_mean = mean(y,1); y_std = std(y,[],1)/sqrt(size(y,1)).*1.96;
%  lineProps.col = tmpcols;
%  h=mseb(fw,y_mean,y_std,lineProps,1);
for z=1:4
    d=del(z);
    plot(fw, tstat(d,:,i),'lineWidth',3,'Color',tmp_col(z+2,:));
    hold on
end 
lgd = legend({'Envelope','Acoustic landmarks},'FontSize',9,'TextColor','black','Location','northoutside');
lgd.NumColumns = 1;
xlim([fw(end) 20]);
legend boxoff 

%% NMF
% generate the data structure
dat = zeros(24,nlab,Ndel,length(fw));

block_dat = zeros(6,nlab,Ndel,length(fw));
in_dat = zeros(24,nlab,Ndel,length(fw));
sds=[];
Ndel = 61;
% gather the data in one matrix and average
for ss=1:numel(pps)
    fprintf('subject %d \n', ss);
    for bb=2:7
        if bb==8
            continue
        else
        %load([homeDir '/results/lcmv/aal/copula_mi/lcmv_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
        load([homeDir '/results/lcmv/hcpmmp1/copula_mi/lcmv_denv_tpca_hcpmmp1_mi_copula_del_s' pps{ss} '_run' num2str(bb) '.mat'],'allmi');
         for ilab=1:nlab
%            indptf(ilab,:,:)=normalize(squeeze(indptf(ilab,:,:)),1,'zscore');         
        block_dat(bb-1,ilab,:,:)=allmi{ilab};
         end
         end
    end
in_dat(ss,:,:,:) = mean(block_dat,1);
 tmp=reshape(in_dat(ss,:,:,:),[1 nlab*Ndel*length(fw)]);
% % % % XXX try normalisatiions 
%   sd=std(tmp);
%   sds=[sds; sd];
%   in_dat(ss,:,:,:)=in_dat(ss,:,:,:)/sd;
end

if replicate==1
    dat=permute(in_dat,[1 3 4 2]);
    
    for isub = 1:numel(pps)
        tmp = dat(isub,:,:,:);
        [u s v p] = ChoosingR(squeeze(tmp(1,10,:,:))); % for a fixed 100ms delay / was checked for all delays - same result
        allp(isub,d) = p;
    end
    
    dat=reshape(dat,[length(pps)*Ndel*length(fw) length(atemplate.parcellationlabel)]);
    dat(dat<0)=0;
    
    c = median(allp);
    c=15; % fixed number of components
    
    w={};
    h={};
    parfor k=1:300
        k
        opt =statset('maxiter',10,'display','final');
        [w0,h0] = nnmf(dat,c,'rep',10,'opt',opt,'alg','mult');
        opt = statset('maxiter',1000,'display','final');
        [w{k},h{k},D(k)] = nnmf(dat,c,'w0',w0,'h0',h0,'opt',opt,'alg','als');
    end
    [o,j]=sort(D);
    W=w{j(1)};
    H=h{j(1)};
else
%     load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/nmf/www.mat');
%     load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/nmf/hhh.mat');
%     load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/nmf/DDD.mat');
    load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/nmf/h_nmf15_negdel.mat');
    load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/nmf/w_nmf15_negdel.mat');
    load('/remotedata/AgGross/cspeech/results/lcmv/hcpmmp1/copula_mi/nmf/D_nmf15_negdel.mat');
    [o,j]=sort(D);
    W=w{j(1)};
    H=h{j(1)};

end
c=15;
W=reshape(W,[length(pps) Ndel length(fw) c]);

%% cluster stats
freq=[];
freq.label=atlas_8k.parcellationlabel;
freq.freq=[1:64];
freq.powspctrm=[];
freq.dimord='chan_freq';

for ss=1:24
    dat1{ss}=freq;
    dat1{ss}.powspctrm=squeeze(W(ss,:,:,:));
    dat2{ss}=freq;
    dat2{ss}.powspctrm=zeros(61,64,c);%
end
nsubj=24;
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
%            cfg.statistic        = 'ft_statfun_pooledT';
cfg.correctm         = 'fdr'; %cluster
cfg.parameter = 'powspctrm';
cfg.alpha            = 0.05;
cfg.clusterstatistic = 'wcm';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 0;          % minimum number of neighborhood frechannels that is
                                   % required for a selected sample to be included                                   % in the clustering algorithm (default=0)
cfg.tail             = 1;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 1;
cfg.neighbours=[];
%cfg.neighbours=neighbours;
cfg.numrandomization =5000;  %set this to 5000 for final statistics (correcting for multiple comparisons)
% set up design matrix
cfg.design  = [1:nsubj 1:nsubj; ones(1, nsubj) 2*ones(1, nsubj)];
cfg.uvar    = 1;
cfg.ivar    = 2;
% run t-test
stat = ft_freqstatistics(cfg, dat1{:}, dat2{:});
stat.stat = stat.stat.*stat.mask;
tstat = stat.mask.*stat.stat;
tstat(tstat==0)=NaN;


cmaps={'Blues','Purples','Greens','Greys','Oranges','RdPu','PuBu','Reds','YlGn','OrRd','YlGnBu','BuPu','YlOrBr','YlOrRd','PuBuGn','BrBG'};
%cmaps={'BrBG','PRGn','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','BrBG','PRGn','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral'};
% plot spectra / delays per component 
 

delay_samples = [-30:30];
%delay_samples = [1, 10, 20, 30];
delays = delay_samples*10;

del=[32 41 51 61];
pos_del=[1:10:length(delay_samples)];

%W=reshape(W,[length(pps) Ndel length(fw) c]);
for icom = 1:c
%plot nmf spectra
    fig(icom)=figure;
    sgtitle(sprintf('component %d',icom));
    sp1=subplot(1,2,1);
    tmp_col=brewermap(6,cmaps{icom});
    %tmp_col=magma(10);
    dels = [1:numel(delay_samples)];
    z=1;
    for d=del
        plot(fw, tstat(d,:,icom),'lineWidth',3,'Color',tmp_col(z+2,:));
        hold on
        z=z+1;
    end 
    lgd = legend({'10 ms','100 ms','200 ms','300 ms'},'FontSize',9,'TextColor','black','Location','northoutside');
    lgd.NumColumns = 2;
    xlim([fw(end) 15]);
    ylim([1.5 10]);
    legend boxoff 
    %set(gca,'box','off')
    ax = gca;
    xlabel('Frequency (Hz)') 
    %%ylabel('t-values') 
    ax.FontSize = 11;
    possub  = get(ax,'Position');

    % plot component frequency/delays
    sp2=subplot(1,2,2);
    tmpsum = nansum(tstat(:,:,icom));
    [v,id]=find(tmpsum);
    h=pcolor(delays,fw,tstat(:,:,icom)');
    col=colorbar;
    ylabel(col, 't-values')
    col.FontSize = 11;
    h.EdgeColor = 'none';
    col.Limits = [1.5 8];
    colormap(brewermap(nlab,cmaps{icom}));
    %xlim([fw(end) 10]);
    ylim([1 15]);
    ax = gca;
    ax.FontSize = 11;
    xlabel('Delays (ms)'); 
    ylabel('Frequency (Hz)');
    possub  = get(ax,'Position');
    %ax.Position = [possub(1) possub(2) possub(3) possub(4)-0.5];
    set(gca,'box','off')
    %set(gcf, 'Units', 'Normalized', 'Position',[0.2000 0.6833 0.1873 0.1042]);
  %  set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.03, 0.132, 0.26]);
    possub1  = get(sp1,'Position');
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.7263 0.4323 0.1710 0.1552]);
    print(sprintf('%d_component_spctra',icom),'-djpeg','-r300');

    % for plotting all components in one cortical space run with c=15
%     for c=15:15
    Hth=zeros(size(H));
    for icomp=icom:icom % icom:icom 
        sigind = find(H(icomp,:) > prctile(H(icomp,:),95));
        Hth(icomp,sigind)=H(icomp,sigind);
    end

    atlas.left=zeros(1,length(atlas.indexmax));
    atlas.right=zeros(1,length(atlas.indexmax));
    for ilab=1:nlab
        [m1,m2]=max(Hth(:,ilab));
        str_label = atlas_8k.parcellationlabel{ilab};
        index = find(contains(atlas.indexmaxlabel,str_label));
        if ilab == 1 || ilab == 182
            continue
        else
        tmp_roiidx=find(atlas.indexmax==index);
        tmp_pos = atlas.brainstructure(tmp_roiidx);
        end
        if m1>0
            if tmp_pos(1) == 2
                
                atlas.right(tmp_roiidx)=m2;
            elseif tmp_pos(1) == 1
                
                atlas.left(tmp_roiidx)=m2;
            end
        else
            continue
        end
    end

    %plot right hemisphere (both)
    %figure;
    cmap=[];
    cmap = [[1 1 1]; tmp_col(5,:)];
    tmp_col=brewermap(6,cmaps{icom});

    figure;
    hcp_allviews_nmf(atlas,sourcemodel,sourcemodel_left,sourcemodel_right,cmap)

    clear cmap;
    clear tmp;
    print(sprintf('%d_component_surfaces',icom),'-djpeg','-r300');
   % print(fig,'-dpsc2','-append' ,['sparseNMF_comps.ps']);
   %pause;
   close all;
end

% plot the nmf coeeficient spectra (for supplemmentary)
del=[32 41 51 61];
for icom=1:15
    subplot(4,4,icom)
    tmp_col=brewermap(6,cmaps{icom});
    z=1;
    for d=del
    hold on
    plot(fw,squeeze(mean(W(:,d,1:length(fw),icom),1)),'lineWidth',3,'color',[tmp_col(z+2,:)]);
    title(sprintf('component %d',icom));
    ax = gca;
    set(gca,'box','off')
    xlabel('Frequency (Hz)')
    ylabel('NMF coefficients')
    xlim([fw(end) 15])
    ax.FontSize = 11;
    z=z+1;
    end
    lgd = legend({'10 ms','100 ms','200 ms','300 ms'},'FontSize',9,'TextColor','black','Location','northoutside');
    lgd.NumColumns = 2;
    legend boxoff 

end
print(sprintf('MI_nmf_supplementary1'),'-djpeg','-r300');

%% all nmf surfaces plot
Hth=zeros(size(H));
for icomp=1:c % icom:icom 
    sigind = find(H(icomp,:) > prctile(H(icomp,:),95));
    Hth(icomp,sigind)=H(icomp,sigind);
end

atlas.left=zeros(1,length(atlas.indexmax));
atlas.right=zeros(1,length(atlas.indexmax));
for ilab=1:nlab
    [m1,m2]=max(Hth(:,ilab));
    str_label = atlas_8k.parcellationlabel{ilab};
    index = find(contains(atlas.indexmaxlabel,str_label));
    if ilab == 1 || ilab == 182
        continue
    else
    tmp_roiidx=find(atlas.indexmax==index);
    tmp_pos = atlas.brainstructure(tmp_roiidx);
    end
    if m1>0
        if tmp_pos(1) == 2
            atlas.right(tmp_roiidx)=m2;
        elseif tmp_pos(1) == 1         
            atlas.left(tmp_roiidx)=m2;
        end
    else
        continue
    end
end

%figure;

cmap=[];
for i=1:(c+1)
    if i==1
        tmp=brewermap(6,cmaps{i});
         cmap=[cmap; tmp(1,:)];
    else
    tmp=brewermap(6,cmaps{i-1});
    cmap=[cmap; tmp(5,:)];
    end
end

figure;
hcp_allviews_nmf(atlas,sourcemodel,sourcemodel_left,sourcemodel_right,cmap)

    clear cmap;
    clear tmp;
    print(sprintf('all_components_surfaces',icom),'-djpeg','-r300');
   % print(fig,'-dpsc2','-append' ,['sparseNMF_comps.ps']);
   %pause;
   close all;

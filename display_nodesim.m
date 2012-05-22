%==========================%
% FUNCTION DISPLAY_NODESIM %
%==========================%
%
%DISPLAY_NODESIM displays figures related to statistics on similarity between different
% categories of pairs of probe sets calculated on a particular network.
%
%INPUT PARAMETERS
% 1 ProbeNbLimit: minimal number of probes targeting a gene
%                 (used to make a bipartition of genes)
% 2     TestFlag: 
% 3      ProbeNb: number of probes in a probe set
% 4    ChipRank: chip ChipRank
% 5      NetRank: network rank
% 6    FigRanks: list of figures to be displayed
% 7       SimNb:  number of similarity types to be displayed in figure 4 (from 1 to SimNb)
%                 order is (Sim, OutSim, HighSim, LowHighSim, LowSim)
% 8  AceviewFlag: indicates if Aceview is used
%
%FIGURES 19 to 29

%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

function display_nodesim(ProbeNbLimit,TestFlag,Species,ChipRank,NetRank,FigRanks,AceviewFlag,PairCorrLimit)
global K
DataDir=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(DataDir)
ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
if isempty(ChipPos)
    h=errordlg(sprintf('chip m%u does not exist',ChipRank));
    waitfor(h)
    error('process canceled')
end
ProbeNb=K.chip.probeNb(ChipPos);

%% Load SIM structures
%            Sim.corr: positive correlations values for each pair of probesets
%            Sim.anti: negative correlations values for each pair of probesets
%     Sim.firstPsRank: probeset ranks of the first probeset in each pair
%       Sim.sndPsRank: probeset ranks of the second probeset in each pair
% Sim.commonNodeNb{i}: common number of neighbouhrs
%  Sim.firstNodeNb{i}: number of neighbouhrs for the first probeset
%    Sim.sndNodeNb{i}: number of neighbouhrs for the second probeset
%      Sim.overlap{i}: percentage overlap (common number * 100 / min(firstNodeNb,sndNodeNb))
%           Sim.pv{i}: p-value calculated with hypergeometric distribution;

if TestFlag
    % Pairs of probe sets targeting multiple genes
    SIM=cell(2,5);
    for SimL1=1:2
        for SimL2=1:5
            SIM{SimL1,SimL2}=[];
        end
    end

    % Pairs of probe sets targeting a single gene
    % Pairs of probe sets inside exons of the same gene(s)
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_single.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{1,1}=Sim;
    end
    % Pairs of probe sets outside exons of the same gene(s)
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_single_testout.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{1,2}=Sim;
    end
    % Pairs of randomly matched probe sets present in ~single and targeting genes with more
    % or = than max(1,ProbeNbLimit-2 probes)
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_single_testhigh.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{1,3}=Sim;
    end
    % Pairs of randomly matched probe sets, one present ~testhigh and the other in ~testLow
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_single_testlowhigh.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{1,4}=Sim;
    end
    % Pairs of randomly matched probe sets present in ~single targeting with genes with less
    % than 3 probes
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_single_testlow.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{1,5}=Sim;
    end

    % Pairs of probe sets inside exons of the same gene(s)
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_multiple.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{2,1}=Sim;
    end
    % Pairs of probe sets outside exons of the same gene(s)
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testout.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{2,2}=Sim;
    end
    % Pairs of randomly matched probe sets present in ~mutliple and targeting genes with more
    % or = than max(1,ProbeNbLimit-2 probes)
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testhigh.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{2,3}=Sim;
    end
    % Pairs of randomly matched probe sets, one present ~testhigh and the other in ~testLow
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testlowhigh.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{2,4}=Sim;
    end
    % Pairs of randomly matched probe sets present in ~mutliple targeting with genes with
    % less than 3 probes
    FileName=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testlow.mat',ChipRank,NetRank,ProbeNbLimit);
    if exist(FileName,'file')
        load(FileName)
        SIM{2,5}=Sim;
    end

else
    % Pairs of probe sets inside exons of the same gene(s)with at leat ProbeNbLimit probes
    FileName=sprintf('m%u_n%u_nodesim_probenb%u.mat',ChipRank,NetRank,ProbeNbLimit);
    load(FileName)
    SIM{1,1}=Sim;
end
clear Sim


%% load gene information of probe sets
cd(DataDir)
eval(sprintf('load m%u_ensembl_psinfo',ChipRank))
PsInfo{1}=EPsInfo;
clear EPsInfo
if AceviewFlag
    eval(sprintf('load m%u_aceview_psinfo',ChipRank))
    PsInfo{2}=APsInfo;
    clear APsInfo
end

%% load or fill CompInfo files
cd(DataDir)
if TestFlag
    FileName=sprintf('m%u_compinfo_probenb%u_single.mat',ChipRank,ProbeNbLimit);
else
    FileName=sprintf('m%u_compinfo_probenb%u.mat',ChipRank,ProbeNbLimit);
end
if exist(FileName,'file')
    if TestFlag
        FileName=sprintf('m%u_compinfo_probenb%u_single.mat',ChipRank,ProbeNbLimit);
        load(FileName)
        CurrCompInfo=CompInfo;
        FileName=sprintf('m%u_compinfo_probenb%u_multiple.mat',ChipRank,ProbeNbLimit);
        load(FileName)
        CompInfo={};
        CompInfo{1}=CurrCompInfo;
        CompInfo{2}=MCompInfo;
        clear MCompInfo CurrCompInfo
    else
        FileName=sprintf('m%u_compinfo_probenb%u.mat',ChipRank,ProbeNbLimit);
        load(FileName)
        CurrCompInfo=CompInfo;
        CompInfo={};
        CompInfo{1}=CurrCompInfo;
        clear CurrCompInfo
    end
else
    if TestFlag
        GopFlag=1;
        [CompInfo,Ps]=fill_psinfo(PsInfo,ProbeNbLimit,SIM{1,1},SIM{1,2},SIM{1,3},SIM{1,4},SIM{1,5},AceviewFlag,GopFlag);
        cd(DataDir)
        FileName=sprintf('m%u_compinfo_probenb%u_single.mat',ChipRank,ProbeNbLimit);
        eval(sprintf('save %s CompInfo Ps',FileName))

        GopFlag=0;
        [MCompInfo,MPs]=fill_psinfo(PsInfo,ProbeNbLimit,SIM{2,1},SIM{2,2},SIM{2,3},SIM{2,4},SIM{2,5},AceviewFlag,GopFlag);
        cd(DataDir)
        FileName=sprintf('m%u_compinfo_probenb%u_multiple.mat',ChipRank,ProbeNbLimit);
        eval(sprintf('save %s MCompInfo MPs',FileName))

        CurrCompInfo=CompInfo;
        CurrMCompInfo=MCompInfo;
        CompInfo={};
        CompInfo{1}=CurrCompInfo;
        CompInfo{2}=CurrMCompInfo;
        clear CurrCompInfo CurrMCompInfo
    else
        GopFlag=0;
        [CompInfo,Ps]=fill_psinfo(PsInfo,ProbeNbLimit,SIM{1,1},[],[],[],[],AceviewFlag,GopFlag);
        cd(DataDir)
        FileName=sprintf('m%u_compinfo_probenb%u.mat',ChipRank,ProbeNbLimit);
        eval(sprintf('save %s CompInfo Ps',FileName))

        CurrCompInfo=CompInfo;
        CompInfo={};
        CompInfo{1}=CurrCompInfo;
        clear CurrCompInfo
    end
end

PvCorrLimit=SIM{1,1}.pvCorrLimit;
Colors='rmbcg';
Titles{1}={'InSim single','OutSim single','HSim single','LHSim single','LSim single'};
Titles{2}={'InSim multiple','OutSim multiple','HSim multiple','LHSim multiple','LSim multiple'};
PvNb=length(PvCorrLimit);

if TestFlag
    %% FIG 19  pv(sim) for pv corr limit in PvCorrLimit > 0 vs pv(sim) for pv corr limit = 0
    %  if PvCorrLimit=[0,40:60] for example
    Letters='abcdefgh';
    if ~isempty(find(FigRanks==19))
        FigRank=19;
        for PvL=1:length(PvCorrLimit)-1
            h=figure;
            set(h,'name',sprintf('FIG%u%c- Effect of pv corr limit on pv(similarity) - Corr%u vs Corr%u',FigRank,Letters(PvL),PvCorrLimit(1+PvL),PvCorrLimit(1)))
            set(gcf,'color',[1,1,1])
            for PlotL=1:2
                for SimL=1:5
                    subplot(2,5,SimL+(PlotL-1)*5)
                    if ~isempty(SIM{PlotL,SimL})
                    plot(SIM{PlotL,SimL}.pv{1},SIM{PlotL,SimL}.pv{1+PvL},sprintf('%c+',Colors(SimL)))

                    xlabel(sprintf('Pv at corr > %u',PvCorrLimit(1)))
                    ylabel(sprintf('Pv at corr > %u',PvCorrLimit(1+PvL)))
                    set(gca,'xlim',[-350,0])
                    set(gca,'ylim',[-350,0])
                    end
                    set(gca,'box','on')
                    title(Titles{PlotL}{SimL})
                end
            end
            Position=get(gcf,'position');
            Position(3:4)=[1000,350];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(PvL)),'png')
            close(h)
        end
    end

    %% FIG 20  pv(sim) vs corr(probe set pair) for different values of pv corr limit
    if ~isempty(find(FigRanks==20))
        FigRank=20;
        for PvL=1:length(PvCorrLimit)
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('FIG%u%c - p-value vs corr of pairs of probe sets for pv corr limit > %u ',FigRank,Letters(PvL),PvCorrLimit(PvL)))
            for PlotL=1:2
                for SimL=1:5
                    subplot(2,5,SimL+(PlotL-1)*5)
                    if ~isempty(SIM{PlotL,SimL})
                    plot(SIM{PlotL,SimL}.corr,SIM{PlotL,SimL}.pv{PvL},sprintf('%c+',Colors(SimL)))                
                    xlabel('corr')
                    ylabel('pv')
                    set(gca,'ylim',[-350,0])
                    end
                    set(gca,'box','on')
                    title(sprintf('pv(%s) vs CORR',Titles{PlotL}{SimL}));
                end
            end
            Position=get(gcf,'position');
            Position(3:4)=[1000,350];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(PvL)),'png')    
            close(h)
        end
    end

    %% FIG 21  pv(sim) distribution for different values of pair corr limit
    if ~isempty(find(FigRanks==21))
        FigRank=21;
        Lines='-:';
        Legend=cell(length(PairCorrLimit)+1,1);
        for i=1:length(Legend)
            Legend{i}=cell(1,10);
            for j=1:10
                Legend{i}{j}='-1';
            end
        end
        for PvL=1:length(PvCorrLimit)
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('FIG%u%c - Pv distribution at pv corr limit> %u for  different Corr values between probe set of the same pair',FigRank,Letters(PvL),PvCorrLimit(PvL)))
            for PlotL=1:2
                subplot(2,ceil((length(PairCorrLimit)+2)/2),1)
                hold on
                Pos{PlotL}=[];
                for SimL=1:5
                    if ~isempty(SIM{PlotL,SimL})
                        Pos{PlotL}=[Pos{PlotL},SimL];
                        plot(sort(SIM{PlotL,SimL}.pv{PvL}(SIM{PlotL,SimL}.corr==0)),0:1/(length(find(SIM{PlotL,SimL}.corr==0))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))                    
                    end
                end
                if PlotL==2
                    title('CORR=0');
                    xlabel('pv')
                    set(gca,'ylim',[0,1])
                    set(gca,'box','on')
                    legend([Titles{1}(Pos{1}),Titles{2}(Pos{2})],'location','NorthWest')               
                end

                subplot(2,ceil((length(PairCorrLimit)+2)/2),2)
                hold on
                for SimL=1:5
                    if ~isempty(SIM{PlotL,SimL})               
                    plot(sort(SIM{PlotL,SimL}.pv{PvL}(SIM{PlotL,SimL}.corr>=0)),0:1/(length(find(SIM{PlotL,SimL}.corr>=0))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))
                    Legend{1}{(PlotL-1)*5+SimL}=sprintf('%u',length(find(SIM{PlotL,SimL}.corr>=0)));
                    end
                end
                if PlotL==2
                    title('CORR>=0');
                    xlabel('pv')                
                    set(gca,'ylim',[0,1])
                    set(gca,'box','on')
                    CurrLegend=Legend{1};
                    CurrLegend(strmatch('-1',CurrLegend))=[];
                    legend(CurrLegend,'location','NorthWest')
                end

                for CorrL=1:length(PairCorrLimit)
                    %for CorrL=1:2
                    subplot(2,ceil((length(PairCorrLimit)+2)/2),CorrL+2)
                    hold on
                    for SimL=1:5
                        if ~isempty(SIM{PlotL,SimL})                       
                        plot(sort(SIM{PlotL,SimL}.pv{PvL}(SIM{PlotL,SimL}.corr>PairCorrLimit(CorrL))),0:1/(length(find(SIM{PlotL,SimL}.corr>PairCorrLimit(CorrL)))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))
                        Legend{CorrL+1}{(PlotL-1)*5+SimL}=sprintf('%u',length(find(SIM{PlotL,SimL}.corr>PairCorrLimit(CorrL))));
                        end
                    end
                    if PlotL==2
                        title(sprintf('CORR>%u',PairCorrLimit(CorrL)));
                        xlabel('pv')                  
                        set(gca,'ylim',[0,1])
                        set(gca,'box','on')
                        CurrLegend=Legend{CorrL+1};
                        CurrLegend(strmatch('-1',CurrLegend))=[];
                        legend(CurrLegend,'location','NorthWest')
                    end
                end
            end
            Position=get(gcf,'position');
            Position(3:4)=[1000,800];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(PvL)),'png')                    
            close(h)
        end
    end

    %% FIG 22  pv(sim) distribution of pairs with corr>0 for different values of pv corr limit
    if ~isempty(find(FigRanks==22))
        FigRank=22;
        h=figure;
        set(gcf,'color',[1,1,1])    
        set(h,'name',sprintf('FIG%u - Pv distributions',FigRank))
        for PvL=1:length(PvCorrLimit)
            subplot(2,ceil(length(PvCorrLimit)/2),PvL)
            hold on
            for PlotL=1:2
                for SimL=1:5
                    if ~isempty(SIM{PlotL,SimL})
                        plot(sort(SIM{PlotL,SimL}.pv{PvL}(SIM{PlotL,SimL}.corr>0)),0:1/(length(find(SIM{PlotL,SimL}.corr>0))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))
                    end
                end
                if PlotL==2
                    xlabel('pv(sim)')
                    title(sprintf('pv corr limit=%u',PvCorrLimit(PvL)))
                    set(gca,'box','on')
                end
            end
        end 
        Position=get(gcf,'position');
        Position(3:4)=[500,500];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
        close(h)
    end
    %% FIG 23 corr(pair), anti(pair) distribution  for pairs with corr>0
    if ~isempty(find(FigRanks==23))
        FigRank=23;
        h=figure;
        set(gcf,'color',[1,1,1])
        subplot(1,3,1)
        set(h,'name',sprintf('FIG%u - Corr & Anti distributions',FigRank))
        subplot(1,3,1)
        hold on
        for PlotL=1:2
            Pos{PlotL}=[];
            for SimL=1:5
                if ~isempty(SIM{PlotL,SimL})
                    Pos{PlotL}=[Pos{PlotL},SimL];
                    plot(sort(SIM{PlotL,SimL}.corr(SIM{PlotL,SimL}.corr>0)),0:1/(length(find(SIM{PlotL,SimL}.corr>0))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))
                end
            end
            if PlotL==2
                xlabel('CORR')
                title('Corr distribution')
                set(gca,'box','on')
                legend([Titles{1}(Pos{1}),Titles{2}(Pos{2})],'location','NorthWest')
            end
        end

        subplot(1,3,2)
        hold on
        for PlotL=1:2
            for SimL=1:5
                if ~isempty(SIM{PlotL,SimL})
                plot(sort(SIM{PlotL,SimL}.anti(SIM{PlotL,SimL}.corr>0)),0:1/(length(find(SIM{PlotL,SimL}.corr>0))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))
                end
            end
            if PlotL==2
                xlabel('ANTI')
                title('Anti distribution')
                set(gca,'xlim',[0,20])
                set(gca,'box','on')
            end
        end

        subplot(1,3,3)
        hold on
        for PlotL=1:2
            for SimL=1:5
                if ~isempty(SIM{PlotL,SimL})
                plot(sort(SIM{PlotL,SimL}.corr(SIM{PlotL,SimL}.corr>0)-SIM{PlotL,SimL}.anti(SIM{PlotL,SimL}.corr>0)),0:1/(length(find(SIM{PlotL,SimL}.corr>0))-1):1,sprintf('%c%c',Colors(SimL),Lines(PlotL)))
                end
            end
            if PlotL==2
                xlabel('CORR-ANTI')
                title('Corr - Anti distribution')
                set(gca,'box','on')
            end
        end
        Position=get(gcf,'position');
        Position(3:4)=[1000,350];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
        close(h)
    end

    %% FIG 24 properties of Sim pairs

    if ~isempty(find(FigRanks==24))
        FigRank=24;
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('FIG%u - Distributions of all properties for different Sim',FigRank))
        subplot(3,4,1)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);            
                plot(sort(CompInfo{PlotL}{1}{1}.comGeneIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.comGeneIn(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('Common gene nb')
            end
            if PlotL==2
                legend('InSim single','OutSim single','InSim multiple','OutSim multiple')
                set(gca,'xlim',[0,10])
            end

        end


        subplot(3,4,2)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.meanComGeneIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.meanComGeneIn(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('GMean of common genes')
            end
        end

        subplot(3,4,3)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.uncomGeneIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.uncomGeneIn(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('Uncommon gene nb')                      
                set(gca,'xlim',[0,10])
            end
        end

        subplot(3,4,4)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.comTscriptIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if PlotL==2
                xlabel('Common transcript nb')
                set(gca,'xlim',[0,10])
            end
        end

        subplot(3,4,5)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.meanTscriptIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if PlotL==2
                xlabel('Gmean of common transcripts')            
            end
        end

        subplot(3,4,6)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.uncomTscriptIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if PlotL==2
                xlabel('uncommon Transcript Nb')
                set(gca,'xlim',[0,5])
            end
        end

        subplot(3,4,7)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.cMeanGroupProbeIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
                hold on
            end
            if PlotL==2            
                xlabel('GMean of probe nb in common exons')           
            end
        end

        subplot(3,4,8)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.tMeanGroupProbeIn(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
                hold on
            end
            if PlotL==2
                xlabel('GMean of probe nb in all exons')      
            end
        end

        subplot(3,4,9)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.maxProbe1In(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.maxProbe1In(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('max probe 1 Nb')
                set(gca,'xlim',[0,ProbeNb+1])
            end
        end

        subplot(3,4,10)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.maxProbe2In(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
                hold on
            end
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.maxProbe2In(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('max probe 2 Nb')
                %legend('Sim single','OutSim single','HSim single','LSim single','Sim multiple','OutSim multiple','HSim multiple','LSim multiple','location','SouthEast')
            end
        end

        subplot(3,4,11)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.minProbe1In(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.minProbe1In(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('min probe 1 Nb')
                %legend('Sim single','OutSim single','HSim single','LSim single','Sim multiple','OutSim multiple','HSim multiple','LSim multiple','location','SouthEast')
            end
        end

        subplot(3,4,12)
        for PlotL=1:2
            if ~isempty(SIM{PlotL,1})
                Pos=find(SIM{PlotL,1}.corr>0);
                plot(sort(CompInfo{PlotL}{1}{1}.minProbe2In(Pos)),0:1/(length(Pos)-1):1,sprintf('r%c',Lines(PlotL)))
            end
            hold on
            if ~isempty(SIM{PlotL,2})
                Pos=find(SIM{PlotL,2}.corr>0);
                plot(sort(CompInfo{PlotL}{2}{1}.minProbe2In(Pos)),0:1/(length(Pos)-1):1,sprintf('m%c',Lines(PlotL)))
            end
            if PlotL==2
                xlabel('min probe 2 Nb')
                %legend('Sim single','OutSim single','HSim single','LSim single','Sim multiple','OutSim multiple','HSim multiple','LSim multiple','location','SouthEast')
            end
        end
        Position=get(gcf,'position');
        Position(3:4)=[1000,700];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
        close(h)
    end



    %RELATION BETWEEN CORR ANTI and PV
    %shows
    %   that the correlation is a powerful indicator of neighboorh similarity:
    %      when CORR>0 => pv <=-10
    %       ANTI values are inversely correlated to CORR values
    %   that there exist probeset with O correlation (and 0 anti) but a similar neighboorh

    %% FIG 25 pv(sim), corr-anti (pair), anti (pair) ordered on corr or on anti
    if ~isempty(find(FigRanks==25))
        FigRank=25;
        SimL=1;
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('FIG%u - Distributions of all properties for different Sim',FigRank))
        for PlotL=1:2
            for SortL=1:2
                set(gcf,'color',[1,1,1])
                Pos1=1;
                Pos2=2;
                if SortL==1
                    if PlotL==1
                        PlotCol=1;
                    else
                        PlotCol=2;                    
                    end
                else
                    PlotCol=3;
                    if PlotL==2
                       Pos1=3;
                       Pos2=4;
                    end
                end
                if SortL==1
                    %sort on CORR
                    txt1='ordered on corr';
                    txt2='ordered on corr';
                    switch SimL
                        case {1,2}
                            Index=find(CompInfo{PlotL}{SimL}{1}.minProbe1In>=ProbeNbLimit&CompInfo{PlotL}{SimL}{1}.minProbe2In>=ProbeNbLimit);
                            Corr=SIM{PlotL,SimL}.corr(Index);
                            Anti=SIM{PlotL,SimL}.anti(Index);
                            Pv=SIM{PlotL,SimL}.pv{1}(Index);
                            %MinGroupNb=CompInfo{PlotL}{SimL}{1}.mingroupprobe_in(Index);
                            [Corr SortIndex]=sort(Corr);
                            Anti=Anti(SortIndex);
                            Pv=Pv(SortIndex);
                        case {3,4,5}
                            Corr=SIM{PlotL,SimL}.corr;
                            Anti=SIM{PlotL,SimL}.anti;
                            Pv=SIM{PlotL,SimL}.pv{1};
                            [Corr SortIndex]=sort(Corr);
                            Anti=Anti(SortIndex);
                            Pv=Pv(SortIndex);
                    end
                else

                    %sort on ANTI
                    txt1='ordered on anti';
                    txt2='ordered on anti';
                    switch SimL
                        case {1,2}
                            Index=find(CompInfo{PlotL}{SimL}{1}.minProbe1In>=ProbeNbLimit&CompInfo{PlotL}{SimL}{1}.minProbe2In>=ProbeNbLimit);
                            Corr=SIM{PlotL,SimL}.corr(Index);
                            Anti=SIM{PlotL,SimL}.anti(Index);
                            Pv=SIM{PlotL,SimL}.pv{1}(Index);
                            [Anti SortIndex]=sort(Anti);
                            Corr=Corr(SortIndex);
                            Pv=Pv(SortIndex);
                        case {3,4,5}
                            Corr=SIM{PlotL,SimL}.corr;
                            Anti=SIM{PlotL,SimL}.anti;
                            Pv=SIM{PlotL,SimL}.pv{1};
                            [Anti SortIndex]=sort(Anti);
                            Corr=Corr(SortIndex);
                            Pv=Pv(SortIndex);
                    end
                end

                subplot(3,4,(PlotCol-1)*4+Pos1)
                [ax,h1,h2]=plotyy(1:length(Pv),Pv,1:length(Pv),Corr-Anti);
                set(h1,'linestyle','.','color','g','markersize',3)
                set(h2,'linestyle','.','color','r','markersize',3)
                set(get(ax(1),'Ylabel'),'String','pv','color','g')
                set(get(ax(2),'Ylabel'),'String','Corr-Anti','color','r')
                set(ax(2),'YLim',[0,100])
                set(ax(1),'XTick',[])
                set(ax(2),'XTick',[])
                set(ax(1),'XLim',[0,length(Pv)])
                set(ax(2),'XLim',[0,length(Pv)])
                set(ax(1),'ycolor','g')
                set(ax(2),'ycolor','r')
                xlabel(txt1)
                title(Titles{PlotL}{SimL})

                subplot(3,4,(PlotCol-1)*4+Pos2)
                [ax,h1,h2]=plotyy(1:length(Pv),Pv,1:length(Pv),Anti);
                set(h1,'linestyle','.','color','g','markersize',3)
                set(h2,'linestyle','.','color','b','markersize',3)
                set(get(ax(1),'Ylabel'),'String','pv','color','g')
                set(get(ax(2),'Ylabel'),'String','Anti','color','b')
                set(ax(2),'YLim',[0,100])
                set(ax(1),'XTick',[])
                set(ax(2),'XTick',[])
                set(ax(1),'XLim',[0,length(Pv)])
                set(ax(2),'XLim',[0,length(Pv)])
                set(ax(1),'ycolor','g')
                set(ax(2),'ycolor','b')

                xlabel(txt2)


                switch SimL
                    case {1,2}
                        Index=find(CompInfo{PlotL}{SimL}{1}.minProbe1In>=ProbeNbLimit&CompInfo{PlotL}{SimL}{1}.minProbe2In>=ProbeNbLimit);
                        Corr=SIM{PlotL,SimL}.corr(Index);
                        Anti=SIM{PlotL,SimL}.anti(Index);
                        Pv=SIM{PlotL,SimL}.pv{1}(Index);
                        [Pv SortIndex]=sort(Pv);
                        Corr=Corr(SortIndex);
                        Anti=Anti(SortIndex);
                    case {3,4,5}
                        [Pv SortIndex]=sort(SIM{PlotL,1}.pv{1});
                        Corr=SIM{PlotL,SimL}.corr(SortIndex);
                        Anti=SIM{PlotL,SimL}.anti(SortIndex);
                end

                if PlotCol<3
                    subplot(3,4,(PlotCol-1)*4+3)
                    [ax,h1,h2]=plotyy(1:length(Pv),Pv,1:length(Pv),Corr-Anti);
                    set(h1,'linestyle','.','color','g','markersize',3)
                    set(h2,'linestyle','.','color','r','markersize',3)
                    set(get(ax(1),'Ylabel'),'String','pv','color','g')
                    set(get(ax(2),'Ylabel'),'String','Corr-Anti','color','r')
                    set(ax(2),'YLim',[0,100])
                    set(ax(1),'XTick',[])
                    set(ax(2),'XTick',[])
                    set(ax(1),'XLim',[0,length(Pv)])
                    set(ax(2),'XLim',[0,length(Pv)])
                    set(ax(1),'ycolor','g')
                    set(ax(2),'ycolor','r')
                    xlabel('ordered on pv')

                    subplot(3,4,(PlotCol-1)*4+4)
                    [ax,h1,h2]=plotyy(1:length(Pv),Pv,1:length(Pv),Anti);
                    set(h1,'linestyle','.','color','g','markersize',3)
                    set(h2,'linestyle','.','color','b','markersize',3)
                    set(get(ax(1),'Ylabel'),'String','pv','color','g')
                    set(get(ax(2),'Ylabel'),'String','Anti','color','b')
                    set(ax(2),'YLim',[0,100])
                    set(ax(1),'XTick',[])
                    set(ax(2),'XTick',[])
                    set(ax(1),'XLim',[0,length(Pv)])
                    set(ax(2),'XLim',[0,length(Pv)])
                    set(ax(1),'ycolor','g')
                    set(ax(2),'ycolor','b')
                    xlabel('ordered on pv')
                end
            end
        end
        Position=get(gcf,'position');
        Position(3:4)=[1000,700];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
        close(h)

    end

    %% FIG 26 Histogram of node number and of overlap between probe sets of a pair at different
    %        values of pv corr limit
    if ~isempty(find(FigRanks==26))
        FigRank=26;
        for PvL=1:length(PvCorrLimit)
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('FIG%u%c - Histogram node nb and overlap for pv corr limit = %u',FigRank,Letters(PvL),PvCorrLimit(PvL)))
            for PlotL=1:2
                GoodIndex=find(SIM{PlotL,1}.corr>=50);
                BadIndex=find(SIM{PlotL,1}.corr==0);

                Corr=SIM{PlotL,1}.corr;
                NodeNb1=SIM{PlotL,1}.firstNodeNb{PvL};
                NodeNb2=SIM{PlotL,1}.sndNodeNb{PvL};
                NodeNb=min([NodeNb1,NodeNb2],[],2);
                Overlap=SIM{PlotL,1}.overlap{PvL};
                Overlap(find(Overlap<0))=0;
                Overlap(isinf(Overlap))=0;

                subplot(2,4,(PlotL-1)*4+1)
                hist(NodeNb(GoodIndex),50)
                title(sprintf('corr>=50 (%u%% of %u))',round(length(GoodIndex)*100/length(Corr)),length(Corr)))
                xlabel(sprintf('MIN NODE NB (med=%u)',round(median(NodeNb(GoodIndex)))))
                if PlotL==1
                    ylabel('frequency in SINGLE InSim')
                else
                    ylabel('frequency in MULTIPLE InSim')
                end

                subplot(2,4,(PlotL-1)*4+2)
                hist(NodeNb(BadIndex),50)
                title('')
                title(sprintf('corr==0 (%u%% of %u))',round(length(BadIndex)*100/length(Corr)),length(Corr)))
                xlabel(sprintf('MIN NODE NB (med=%u)',round(median(NodeNb(BadIndex)))))

                subplot(2,4,(PlotL-1)*4+3)
                hist(Overlap(GoodIndex),50)
                title(sprintf('corr>=50 (%u%% of %u))',round(length(GoodIndex)*100/length(Corr)),length(Corr)))
                xlabel(sprintf('OVERLAP (med=%.1f)',median(Overlap(GoodIndex))))

                subplot(2,4,(PlotL-1)*4+4)
                hist(Overlap(BadIndex),50)
                title(sprintf('corr==0 (%u%% of %u))',round(length(BadIndex)*100/length(Corr)),length(Corr)))
                xlabel(sprintf('OVERLAP (med=%.1f)',median(Overlap(BadIndex))))
            end
            Position=get(gcf,'position');
            Position(3:4)=[1000,450];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(PvL)),'png')
            close(h)
        end
    end

    %% FIG 27 Plots of node number and of overlap between probe sets of a pair
    %        at different values of pv corr limit
    if ~isempty(find(FigRanks==27))
        FigRank=27;
        for PlotL=1:2
            h=figure;
            set(gcf,'color',[1,1,1])
            if PlotL==1
                set(h,'name',sprintf('FIG%ua -NODE NB and OVERLAP vs PV in SINGLE',FigRank))
            else
                set(h,'name',sprintf('FIG%ub -NODE NB and OVERLAP vs PV in MULTIPLE',FigRank))
            end
            for PvL=1:PvNb
                BadBindex=SIM{PlotL,1}.corr==0;
                GoodBindex=SIM{PlotL,1}.corr>=50;
                InterBindex=SIM{PlotL,1}.corr>0&SIM{PlotL,1}.corr<50;
                NodeNb1=SIM{PlotL,1}.firstNodeNb{PvL};
                NodeNb2=SIM{PlotL,1}.sndNodeNb{PvL};
                NodeNb=min([NodeNb1,NodeNb2],[],2);
                Overlap=SIM{PlotL,1}.overlap{PvL};
                Pv=SIM{PlotL,1}.pv{PvL};

                %sort data on pv
                [Pv SortIndex]=sort(Pv);
                NodeNb=NodeNb(SortIndex);
                Overlap=Overlap(SortIndex);
                BadBindex=BadBindex(SortIndex);
                GoodBindex=GoodBindex(SortIndex);
                InterBindex=InterBindex(SortIndex);
                BadIndex=find(BadBindex);
                GoodIndex=find(GoodBindex);
                InterIndex=find(InterBindex);            

                subplot(2,PvNb,PvL)
                hold on
                %plot(1:length(Pv),NodeNb,'k.','markersize',3)
                %plot(GoodIndex,NodeNb(GoodIndex),'g.','markersize',3)
                plot(BadIndex,NodeNb(BadIndex),'r.','markersize',3)
                plot(GoodIndex,NodeNb(GoodIndex),'g.','markersize',3)            
                plot(InterIndex,NodeNb(InterIndex),'k.','markersize',3)
                title(sprintf('min node nb at pv corr limit = %u',PvCorrLimit(PvL)))
                xlabel('ordered on pv')
                set(gca,'xlim',[0,length(Pv)])
                ylabel('node nb')
                set(gca,'box','on')
                if PvL==PvNb
                    legend({'corr(pair)>0&<50','corr(pair)>50','corr(pair)=0'});
                end

                subplot(2,PvNb,PvNb+PvL)
                hold on
                %plot(1:length(Pv),abs(Overlap*100),'k.','markersize',3)
                %plot(GoodIndex,abs(Overlap(GoodIndex)*100),'g.','markersize',3)
                plot(BadIndex,abs(Overlap(BadIndex)*100),'r.','markersize',3)
                plot(GoodIndex,abs(Overlap(GoodIndex)*100),'g.','markersize',3)            
                plot(InterIndex,abs(Overlap(InterIndex)*100),'k.','markersize',3)
                plot(1:length(Pv),Pv,'g-')
                title(sprintf('overlap at pv corr limit = %u',PvCorrLimit(PvL)))
                xlabel('ordered on pv')
                set(gca,'xlim',[0,length(Pv)])
                set(gca,'ylim',[-100,100])
                ylabel('overlap and pv(sim)')
                set(gca,'box','on')
            end
            Position=get(gcf,'position');
            Position(3:4)=[1000,500];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            if PlotL==1
                saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
                close(h)
            else
               saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
               close(h)
            end

        end
    end



    %% FIG 28 Distribution of probe nb in different pv, corr and anti value combinations

    if ~isempty(find(FigRanks==28))
        FigRank=28;
        for PlotL=1:2
            h=figure;
            set(gcf,'color',[1,1,1])
            if PlotL==1
                set(h,'name',sprintf('FIG%ua - MIN(Probe Nb) in SINGLE',FigRank))
            else
                set(h,'name',sprintf('FIG%ub - MIN(Probe Nb) in MULTIPLE',FigRank))
            end

            Plot=0;
            for PvL=[0,-50,-100]
                for CorrL=50:10:90
                    %for AntiL=30:-10:0
                    for AntiL=10
                        Plot=Plot+1;
                        Pos=find(SIM{PlotL,1}.pv{1}<PvL&SIM{PlotL,1}.corr>=CorrL&SIM{PlotL,1}.anti<=AntiL);
                        MinProbeNb=min([CompInfo{PlotL}{1}{1}.minProbe1In(Pos),CompInfo{PlotL}{1}{1}.minProbe2In(Pos)],[],2);
                        Val=histc(MinProbeNb,1:ProbeNb);
                        subplot(3,5,Plot)
                        %values with ProbeNb are not displayed (high value masks
                        %the detail ot other probeNb
                        if PlotL==1
                            bar(Val(1:ProbeNb-1)*100./length(Pos))
                        else
                            %values with ProbeNb and 1 are not displayed (high and low value masks
                            %the detail ot other probeNb
                            Val(1)=0;
                            bar(Val(1:ProbeNb-1)*100./length(Pos))
                        end
                        set(gca,'xlim',[0,ProbeNb])
                        set(gca,'xtick',[1,5,10])
                        set(gca,'xticklabel',{'1','5','10'})
                        %set(gca,'ylim',[0,6])
                        title(sprintf('pv<%d c>=%u a<=%u #%u(%u%%)',PvL,CorrL,AntiL,length(Pos),round(Val(ProbeNb)*100/length(Pos))))
                    end
                end
            end

            Position=get(gcf,'position');
            Position(3:4)=[1000,700];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            if PlotL==1
                saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
                close(h)
            else
                saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
                close(h)
            end
        end
    end


    %% FIG 29 Enrichment of probe nb in different pv, corr and anti value combinations
    if ~isempty(find(FigRanks==29))
        FigRank=29;
        for PlotL=1:2
            CurrProbeNb=min([CompInfo{PlotL}{1}{1}.minProbe1In,CompInfo{PlotL}{1}{1}.minProbe2In],[],2);
            ProbeCat=histc(CurrProbeNb,1:ProbeNb);
            ProbeCat=ProbeCat*100/sum(ProbeCat);    
            h=figure;
            set(gcf,'color',[1,1,1])
            if PlotL==1
                set(h,'name',sprintf('FIG%ua - Enrichment of MIN(Probe Nb) in SINGLE',FigRank))
            else
                set(h,'name',sprintf('FIG%ub - Enrichment of MIN(Probe Nb) in MULTIPLE',FigRank))
            end

            Plot=0;
            for PvL=[0,-50,-100]
                for CorrL=50:10:90
                    %for AntiL=30:-10:0
                    for AntiL=10
                        Plot=Plot+1;
                        Pos=find(SIM{PlotL,1}.pv{1}<PvL&SIM{PlotL,1}.corr>=CorrL&SIM{PlotL,1}.anti<=AntiL);
                        MinProbeNb=min([CompInfo{PlotL}{1}{1}.minProbe1In(Pos),CompInfo{PlotL}{1}{1}.minProbe2In(Pos)],[],2);
                        Val=histc(MinProbeNb,1:ProbeNb);
                        try
                        Enrichment=(Val*100/sum(Val))./ProbeCat;
                        catch
                        Enrichment=(Val*100/sum(Val))./ProbeCat';
                        end
                        NegPos=find(Enrichment<1);
                        Enrichment(NegPos)=-1./Enrichment(NegPos);
                        subplot(3,5,Plot)
                        bar(Enrichment)
                        set(gca,'xlim',[1,ProbeNb+1])
                        set(gca,'ylim',[-3,3])
                        set(gca,'ytick',[-3:3])
                        set(gca,'ygrid','on')
                        title(sprintf('pv<%d c>=%u a<=%u #%u',PvL,CorrL,AntiL,length(Pos)))
                    end
                end
            end

            Position=get(gcf,'position');
            Position(3:4)=[1000,700];
            set(gcf,'position',Position)
            try
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            catch
                mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            if PlotL==1
                saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
                close(h)
            else
                saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
                close(h)
            end
        end
    end
end
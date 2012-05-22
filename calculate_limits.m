%============================%
% FUNCTION CALCULATE_LIMITS  %
%============================%

% CALCULATE_LIMITS uses several small networks (1024 comparisons : 32x32 biol cond)
% to find the distributions of p-value of overlaping between neighbourhood and of
% positive and negative correlation of probe sets that target the same groups of transcripts
%
%INPUT PARAMETERS
%
% 1        Species: species
% 2       ChipRank: chip model rank
% 3   ProbeNbLimit: minimal number of probes targeting a gene
% 4       FigRanks: indicates figures that are to be displayed (set it to [] in order
% 5  FirstNetRanks: first list of networks
% 6     PvCorrRank: pv(overlap) is calculated for corr limit >[0,40,50,60]. PvCorrRank
%                   indicates the corr limit to be used, by giving its index in
%                   the corr list(,[0,40,50,60])
% 7        ValFlag: indicates if all values of corr, anti and pv of each pair are used, or
%                   only a single derived value (mean - std)
% 8       MeanFlag: indicates if limit is calculated from mean of four values (single and
%                   multiple targeted genes, and InSim and OutSim) or only from single
%                   targeted genes and InSim.
% varargin:
% 9     NoQlimitFlag: indicates if the second group of network is a set of no qlimit network
% 10   SndNetRankList: second list of networks
%
%FIGURES 30 to 39

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

function [Limit]=calculate_limits(Species,ChipRank,ProbeNbLimit,FigRanks,FirstNetRanks,PvCorrRank,ValFlag,MeanFlag,varargin)
global K

% PS ASSIGNATION
NetNb=[];
NetRanks={};
DataDir=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(DataDir)

NetRanks{1}=FirstNetRanks;
NetNb(1)=length(NetRanks{1});
NetRankListNb=1;

ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
PsNb=K.chip.probesetNb(ChipPos);
ProbeNb=K.chip.probeNb(ChipPos);

if nargin==10
    NetRankListNb=2;
    NoQlimitFlag=varargin{1};
    if NoQlimitFlag
        if length(varargin{2})~=length(FirstNetRanks)
            errodlg('NetRankss do not have the same length')
            error('process canceled')
        end
    end
    NetRanks{2}=varargin{2};
    NetNb(2)=length(NetRanks{2});
    %     NetPosList{2}=zeros(NetRankListNb,NetNb(2));
    %     for j=1:NetNb(2)
    %         NetPosList{2}(j)=find(K.net{ChipRank}.rank==NetRanks{2}(j));
    %     end
end

%load Sim
SIM=cell(NetRankListNb);
for ListL=1:NetRankListNb
    SIM{ListL}=cell(1,2);
    for NetL=1:NetNb(ListL)
        eval(sprintf('load m%u_n%u_nodesim_probenb%u_single;',ChipRank,NetRanks{ListL}(NetL),ProbeNbLimit));
        SIM{ListL}{1,1}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb%u_single_testout;',ChipRank,NetRanks{ListL}(NetL),ProbeNbLimit));
        SIM{ListL}{1,2}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb%u_single_testhigh;',ChipRank,NetRanks{ListL}(NetL),ProbeNbLimit));
        SIM{ListL}{1,3}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb%u_multiple;',ChipRank,NetRanks{ListL}(NetL),ProbeNbLimit));
        SIM{ListL}{2,1}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb%u_multiple_testout;',ChipRank,NetRanks{ListL}(NetL),ProbeNbLimit));
        SIM{ListL}{2,2}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb%u_multiple_testhigh;',ChipRank,NetRanks{ListL}(NetL),ProbeNbLimit));
        SIM{ListL}{2,3}{NetL}=Sim;
    end
end

if ~isempty(FigRanks)
    eval(sprintf('load m%u_pearson_probenb%u_single;',ChipRank,ProbeNbLimit));
    p{1}=Pearson;
    eval(sprintf('load m%u_pearson_probenb%u_multiple;',ChipRank,ProbeNbLimit));
    p{2}=Pearson;
    Pearson={};
    Pearson=p;
    clear p;
    Letters='abcdefghijklmnop';
end

%load CompInfo
load(sprintf('m%u_compinfo_probenb%u_single.mat',ChipRank,ProbeNbLimit));
CoupleInfo{1}=CompInfo;
load(sprintf('m%u_compinfo_probenb%u_multiple.mat',ChipRank,ProbeNbLimit));
CoupleInfo{2}=MCompInfo;
clear CompInfo MCompInfo


%calculate statistics on SIM
%TypeL unique/multiple
for TypeL=1:2
    ComGene{TypeL}=[];
    MeanComGene{TypeL}=[];
    UncomGene{TypeL}=[];
    ComTscript{TypeL}=[];
    MeanComTscript{TypeL}=[];
    UncomTscript{TypeL}=[];
    CMeanProbe{TypeL}=[];
    TMeanProbe{TypeL}=[];
    MaxProbeIn1{TypeL}=[];
    MaxProbeIn2{TypeL}=[];
    MinProbeIn1{TypeL}=[];
    MinProbeIn2{TypeL}=[];
    %SimL in exon/out exon/random couples
    for SimL=1:3
        ComGene{TypeL}{SimL}=[];
        MeanComGene{TypeL}{SimL}=[];
        UncomGene{TypeL}{SimL}=[];
        ComTscript{TypeL}{SimL}=[];
        MeanComTscript{TypeL}{SimL}=[];
        UncomTscript{TypeL}{SimL}=[];
        CMeanProbe{TypeL}{SimL}=[];
        TMeanProbe{TypeL}{SimL}=[];
        MaxProbeIn1{TypeL}{SimL}=[];
        MaxProbeIn2{TypeL}{SimL}=[];
        MinProbeIn1{TypeL}{SimL}=[];
        MinProbeIn2{TypeL}{SimL}=[];

        Corr{1}=[];
        Anti{1}=[];
        Pv{1}=[];
        %keep all significative values for qlimit nets
        Corrs{TypeL}{SimL}={};
        Antis{TypeL}{SimL}={};
        Pvs{TypeL}{SimL}={};
        %recover Corr,Anti and Pv values in qlimit networks
        for NetL=1:NetNb(1)
            Corr{1}=[Corr{1},SIM{1}{TypeL,SimL}{NetL}.corr];
            Anti{1}=[Anti{1},SIM{1}{TypeL,SimL}{NetL}.anti];
            Pv{1}=[Pv{1},SIM{1}{TypeL,SimL}{NetL}.pv{PvCorrRank}];
        end
        InfPos=isinf(Pv{1});
        NegPos=Pv{1}<0;
        Pv{1}(InfPos)=350;
        Pv{1}(InfPos&NegPos)=-350;
        %recover eventually Corr,Anti and Pv values in the second set of networks
        if NetRankListNb==2
            Corr{2}=[];
            Anti{2}=[];
            Pv{2}=[];
            for NetL=1:NetNb(2)
                Corr{2}=[Corr{2},SIM{2}{TypeL,SimL}{NetL}.corr];
                Anti{2}=[Anti{2},SIM{2}{TypeL,SimL}{NetL}.anti];
                Pv{2}=[Pv{2},SIM{2}{TypeL,SimL}{NetL}.pv{PvCorrRank}];
            end
            InfPos=isinf(Pv{2});
            NegPos=Pv{2}<0;
            Pv{2}(InfPos)=350;
            Pv{2}(InfPos&NegPos)=-350;
        end

        if NetRankListNb==2 & NoQlimitFlag
            RoundNb=2;
        else
            RoundNb=1;
        end

        for RoundL=1:RoundNb
            % STATISTICS WHERE DOES EXIST AT LEAST ONE SIGNIFICATIVE CORR VALUE
            % Pos index of probe set pairs that have at least one significative corr value for
            % the current SIM
            PosIndex=find(sum(Corr{RoundL},2));
            % empty variables for first set of (qlimit)  networks
            MeanCorr{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            MeanAnti{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            MeanPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            StdCorr{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            StdAnti{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            StdPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            % for each couple random selection of PosNb networks for calculating Pv on
            % an identical number of networks that have significative Corr values (PosNb)
            RandMeanPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            RandStdPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            OutMeanPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            OutStdPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);

            OutMeanCorr{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            OutMeanAnti{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            OutStdCorr{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            OutStdAnti{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);

            % complement network set of the PosNb randomly selected networks
            RandOutMeanPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            RandOutStdPv{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);


            %nb of networks with significative corr values
            PosNb=zeros(length(PosIndex),1);
            % NB OF REPLICATED VALUES (REPRODUCIBILITY OF PRESENCE
            RepNb{RoundL}{TypeL}{SimL}=zeros(length(PosIndex),1);
            %process each couple and  recover significatives values
            for i=1:size(PosIndex);
                Keep=find(Corr{RoundL}(PosIndex(i),:));
                RepNb{RoundL}{TypeL}{SimL}(i)=length(Keep);
                if RoundL==1
                    Corrs{TypeL}{SimL}{i,1}=Corr{1}(PosIndex(i),Keep);
                    Antis{TypeL}{SimL}{i,1}=Anti{1}(PosIndex(i),Keep);
                    Pvs{TypeL}{SimL}{i,1}=Pv{1}(PosIndex(i),Keep);
                end

                PosNb(i)=length(Keep);
                RandKeep=ceil(rand(PosNb(i),1)*NetNb(1));
                MeanCorr{RoundL}{TypeL}{SimL}(i,1)=mean(Corr{RoundL}(PosIndex(i),Keep));
                MeanAnti{RoundL}{TypeL}{SimL}(i,1)=mean(Anti{RoundL}(PosIndex(i),Keep));
                MeanPv{RoundL}{TypeL}{SimL}(i,1)=mean(Pv{RoundL}(PosIndex(i),Keep));
                RandMeanPv{RoundL}{TypeL}{SimL}(i,1)=mean(Pv{RoundL}(PosIndex(i),RandKeep));

                StdCorr{RoundL}{TypeL}{SimL}(i,1)=std(single(Corr{RoundL}(PosIndex(i),Keep)),1);
                StdAnti{RoundL}{TypeL}{SimL}(i,1)=std(single(Anti{RoundL}(PosIndex(i),Keep)),1);
                StdPv{RoundL}{TypeL}{SimL}(i,1)=std(single(Pv{RoundL}(PosIndex(i),Keep)),1);
                RandStdPv{RoundL}{TypeL}{SimL}(i,1)=std(single(Pv{RoundL}(PosIndex(i),RandKeep)),1);

                OutKeep=find(Corr{RoundL}(PosIndex(i),:)==0);
                RandOutKeep=setdiff(1:NetNb(1),RandKeep);
                OutMeanPv{RoundL}{TypeL}{SimL}(i,1)=mean(Pv{RoundL}(PosIndex(i),OutKeep));
                OutStdPv{RoundL}{TypeL}{SimL}(i,1)=std(single(Pv{RoundL}(PosIndex(i),OutKeep)),1);
                OutMeanCorr{RoundL}{TypeL}{SimL}(i,1)=mean(Corr{RoundL}(PosIndex(i),OutKeep));
                OutMeanAnti{RoundL}{TypeL}{SimL}(i,1)=mean(Anti{RoundL}(PosIndex(i),OutKeep));
                OutStdCorr{RoundL}{TypeL}{SimL}(i,1)=std(single(Corr{RoundL}(PosIndex(i),OutKeep)),1);
                OutStdAnti{RoundL}{TypeL}{SimL}(i,1)=std(single(Anti{RoundL}(PosIndex(i),OutKeep)),1);
                RandOutMeanPv{RoundL}{TypeL}{SimL}(i,1)=mean(Pv{RoundL}(PosIndex(i),RandOutKeep));
                RandOutStdPv{RoundL}{TypeL}{SimL}(i,1)=std(single(Pv{RoundL}(PosIndex(i),RandOutKeep)),1);

            end

            %STATISTICS WHERE DOES EXIST A SIGNIFICATIVE CORR VALUE FOR ALL NETWORKS
            AllPos=find(RepNb{RoundL}{TypeL}{SimL}==NetNb(1));
            AllPos=PosIndex(AllPos);
            %recover Ps Pair informations
            if RoundL==1
                ComGene{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.comGeneIn(AllPos);
                MeanComGene{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.meanComGeneIn(AllPos);
                UncomGene{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.uncomGeneIn(AllPos);
                ComTscript{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.comTscriptIn(AllPos);
                MeanComTscript{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.meanTscriptIn(AllPos);
                UncomTscript{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.uncomTscriptIn(AllPos);
                CMeanProbe{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.cMeanGroupProbeIn(AllPos);
                TMeanProbe{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.tMeanGroupProbeIn(AllPos);
                MaxProbeIn1{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.maxProbe1In(AllPos);
                MaxProbeIn2{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.maxProbe2In(AllPos);
                MinProbeIn1{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.minProbe1In(AllPos);
                MinProbeIn2{TypeL}{SimL}{1}=CoupleInfo{TypeL}{SimL}{1}.minProbe2In(AllPos);
            end

            if ~isempty(FigRanks)
                %STATISTICS WHERE DOES NOT EXIST AT LEAST ONE SIGNIFICATIVE CORR VALUE IN
                %QLIMIT NETWORKS
                AllNullPos=find(sum(Corr{RoundL},2)==0);
                %recover Ps Pair informations
                if SimL~=3 & RoundL==1
                    ComGene{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.comGeneIn(AllNullPos);
                    MeanComGene{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.meanComGeneIn(AllNullPos);
                    UncomGene{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.uncomGeneIn(AllNullPos);
                    ComTscript{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.comTscriptIn(AllNullPos);
                    MeanComTscript{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.meanTscriptIn(AllNullPos);
                    UncomTscript{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.uncomTscriptIn(AllNullPos);
                    CMeanProbe{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.cMeanGroupProbeIn(AllNullPos);
                    TMeanProbe{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.tMeanGroupProbeIn(AllNullPos);
                    MaxProbeIn1{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.maxProbe1In(AllNullPos);
                    MaxProbeIn2{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.maxProbe2In(AllNullPos);
                    MinProbeIn1{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.minProbe1In(AllNullPos);
                    MinProbeIn2{TypeL}{SimL}{3}=CoupleInfo{TypeL}{SimL}{1}.minProbe2In(AllNullPos);
                end

                %recover Pv values in qlimit networks
                AllNullMPv{RoundL}{TypeL}{SimL}=mean(Pv{RoundL}(AllNullPos,:),2);
                AllNullSPv{RoundL}{TypeL}{SimL}=std(single(Pv{RoundL}(AllNullPos,:)),1,2);

                %recover Corr & Anti values in no qlimit networks
                AllNullMCorr{RoundL}{TypeL}{SimL}=mean(Corr{RoundL}(AllNullPos,:),2);
                AllNullSCorr{RoundL}{TypeL}{SimL}=std(single(Corr{RoundL}(AllNullPos,:)),0,2);
                AllNullMAnti{RoundL}{TypeL}{SimL}=mean(Anti{RoundL}(AllNullPos,:),2);
                AllNullSAnti{RoundL}{TypeL}{SimL}=std(single(Anti{RoundL}(AllNullPos,:)),0,2);

            end

        end

    end
end

%% CALCULATE LIMITS
% recover either the list of intermediate limits for each pair of probe sets (Lim,
% calculated with mean +/- std),
% or a single list of all the values (Val)
% finally the limits are drawn from each of these lists ( 5th percentile for corr,
% and 95th percentile for anti and pv(overlap) and are written in LimitLim and
% LimitVal, respectively

AllVal.corr=cell(2,3);
AlltVal.anti=cell(2,3);
AllVal.pv=cell(2,3);

LimitLim.corr=cell(2,3);
LimitLim.anti=cell(2,3);
LimitLim.pv=cell(2,3);

LimitVal.corr=cell(2,3);
LimitVal.anti=cell(2,3);
LimitVal.pv=cell(2,3);

LimitsNb=cell(2,3);
%TypeL= {single,multiple targeted gene)
for TypeL=1:2
    %SimL = (In,Out,Random High)
    for SimL=1:3
        LimitLim.corr{TypeL,SimL}=ones(NetNb,1)*-1;
        LimitLim.anti{TypeL,SimL}=ones(NetNb,1)*-1;
        LimitLim.pv{TypeL,SimL}=ones(NetNb,1)*400;
        LimitVal.corr{TypeL,SimL}=ones(NetNb,1)*-1;
        LimitVal.anti{TypeL,SimL}=ones(NetNb,1)*-1;
        LimitVal.pv{TypeL,SimL}=ones(NetNb,1)*400;

        for RepL=1:NetNb
            %index of values positive in RepL networks
            PosIndex=find(RepNb{1}{TypeL}{SimL}==RepL);
            if ~isempty(PosIndex)
                for ValL=1:3
                    switch ValL
                        case 1
                            CurrVal=Corrs{TypeL}{SimL};
                        case 2
                            CurrVal=Antis{TypeL}{SimL};
                        case 3
                            CurrVal=Pvs{TypeL}{SimL};
                    end
                    CurrLim=zeros(length(PosIndex),1);
                    CurrValues=zeros(length(PosIndex)*RepL,1);
                    for PosL=1:length(PosIndex)                      
                        [CurrVal{PosIndex(PosL)},SortIndex]=sort(CurrVal{PosIndex(PosL)});
                        if ValL==1
                            CurrLim(PosL)=mean(CurrVal{PosIndex(PosL)})-std(single(CurrVal{PosIndex(PosL)}));
                            CurrValues((PosL-1)*RepL+1:PosL*RepL)=CurrVal{PosIndex(PosL)};
                        else
                            CurrLim(PosL)=mean(CurrVal{PosIndex(PosL)})+std(single(CurrVal{PosIndex(PosL)}));
                            CurrValues((PosL-1)*RepL+1:PosL*RepL)=CurrVal{PosIndex(PosL)};
                        end
                    end
                    [CurrValues,Temp]=sort(CurrValues);
                    [CurrLim,Temp]=sort(CurrLim);
                    switch ValL
                        case 1
                            LimitsNb{TypeL,SimL}(RepL)=length(PosIndex);
                            LimitLim.corr{TypeL,SimL}(RepL)=CurrLim(ceil(length(PosIndex)*0.05));
                            AllVal.corr{TypeL,SimL}{RepL}=CurrValues;
                            if ~isempty(CurrValues)
                                CurrPairNb=length(CurrValues);
                                LimitVal.corr{TypeL,SimL}(RepL)=CurrValues(max(1,round(CurrPairNb*0.05)));
                            end
                        case 2
                            LimitLim.anti{TypeL,SimL}(RepL)=CurrLim(ceil(length(PosIndex)*0.95));
                            AllVal.anti{TypeL,SimL}{RepL}=CurrValues;
                            if ~isempty(CurrValues)
                                CurrPairNb=length(CurrValues);
                                LimitVal.anti{TypeL,SimL}(RepL)=CurrValues(max(1,round(CurrPairNb*0.95)));
                            end
                        case 3
                            LimitLim.pv{TypeL,SimL}(RepL)=CurrLim(ceil(length(PosIndex)*0.95));
                            AllVal.pv{TypeL,SimL}{RepL}=CurrValues;
                            if ~isempty(CurrValues)
                                CurrPairNb=length(CurrValues);
                                LimitVal.pv{TypeL,SimL}(RepL)=CurrValues(max(1,round(CurrPairNb*0.95)));
                            end
                    end
                end
            else
                LimitsNb{TypeL,SimL}(RepL)=0;
                LimitLim.corr{TypeL,SimL}(RepL)=0;
                AllVal.corr{TypeL,SimL}{RepL}=0;
                AllVal.corr{TypeL,SimL}{RepL}=0;
                LimitLim.anti{TypeL,SimL}(RepL)=0;
                AllVal.anti{TypeL,SimL}{RepL}=0;
                AllVal.anti{TypeL,SimL}{RepL}=0;
                LimitLim.pv{TypeL,SimL}(RepL)=0;
                AllVal.pv{TypeL,SimL}{RepL}=0;
                AllVal.pv{TypeL,SimL}{RepL}=0;               
            end            
        end
    end
end


%% FIG30 - Distibution of corr, anti and pv(overlap) in LimitVal for pairs of probe sets that
%        target a single gene (FIG1) or multiple genes (FIG2)
if ~isempty(find(FigRanks==30))
    FigRank=30;
    CountFlag=1;
    SimVal=[1,2,3];
    SimLabel={'InSim','OutSim','RandSim'};
    Colors=colors(colormap,NetNb);
    for TypeL=1:2
        h=figure;
        set(gcf,'color',[1,1,1])
        if TypeL==1
            set(h,'name',sprintf('FIG%ua - m%u: Limits for Single',FigRank,ChipRank))
        else
            set(h,'name',sprintf('FIG%ub - m%u: Limits for Multiple',FigRank,ChipRank))
        end
        for SimL=1:length(SimVal)
            if CountFlag
                subplot(length(SimVal),4,(SimL-1)*4+1)
                hold on
                for NetL=1:NetNb               
                    h=bar(NetL,LimitsNb{TypeL,SimL}(NetL));
                    set(h,'facecolor',Colors(NetL,:))
                end
                if NetNb>5
                    set(gca,'xtick',[1,[5:5:NetNb]])
                else
                    set(gca,'xtick',[1:NetNb])
                end
                set(gca,'xlim',[0,NetNb+1])
                set(gca,'tickdir','out')
                set(gca,'box','on')
                if SimL==1
                    title('Nb of positive networks')
                end
            end
            for ValL=1:3
                if CountFlag
                    subplot(length(SimVal),4,(SimL-1)*4+1+ValL)
                else
                    subplot(length(SimVal),3,(SimL-1)*3+ValL)
                end
                hold on
                switch ValL
                    case 1
                        for NetL=1:NetNb
                            if ~isempty(AllVal.corr{TypeL,SimVal(SimL)}{NetL});
                                CurrValues=AllVal.corr{TypeL,SimVal(SimL)}{NetL};
                                CurrPairNb=length(CurrValues);
                                plot(CurrValues,1/CurrPairNb:1/CurrPairNb:1,'color',Colors(NetL,:),'linewidth',2)
                                XVal=CurrValues(max(1,round(CurrPairNb*0.05)));
                                plot(XVal,0.05,'o','color',Colors(NetL,:))
                                line([XVal,XVal],[0,0.05],'color',Colors(NetL,:),'linestyle','-.','linewidth',2)
                            end
                            set(gca,'xlim',[0,100])
                        end
                        if SimVal(SimL)==1
                            title('CORR distribution')
                        end
                    case 2
                        for NetL=1:NetNb
                            if ~isempty(AllVal.corr{TypeL,SimVal(SimL)}{NetL});
                                CurrValues=AllVal.anti{TypeL,SimVal(SimL)}{NetL};
                                CurrPairNb=length(CurrValues);
                                plot(CurrValues,1/CurrPairNb:1/CurrPairNb:1,'color',Colors(NetL,:),'linewidth',2)
                                XVal=CurrValues(round(CurrPairNb*0.95));
                                plot(XVal,0.95,'o','color',Colors(NetL,:))
                                line([XVal,XVal],[0,0.95],'color',Colors(NetL,:),'linestyle','-.','linewidth',2)
                            end
                            set(gca,'xlim',[0,100])
                        end
                        if SimVal(SimL)==1
                            title('ANTI distribution')
                        end

                    case 3
                        for NetL=1:NetNb
                            if ~isempty(AllVal.corr{TypeL,SimVal(SimL)}{NetL});
                                CurrValues=AllVal.pv{TypeL,SimVal(SimL)}{NetL};
                                CurrPairNb=length(CurrValues);
                                plot(CurrValues,1/CurrPairNb:1/CurrPairNb:1,'color',Colors(NetL,:),'linewidth',2)
                                XVal=CurrValues(round(CurrPairNb*0.95));
                                plot(XVal,0.95,'o','color',Colors(NetL,:))
                                line([XVal,XVal],[0,0.95],'color',Colors(NetL,:),'linestyle','-.','linewidth',2)
                            end
                            set(gca,'xlim',[-360,10])
                        end
                        if SimVal(SimL)==1
                            title('PV distribution')
                        end
                        set(gca,'xlim',[-360,10])
                end
                hold on
                set(gca,'box','on')
            end
            if CountFlag
                subplot(length(SimVal),4,(SimL-1)*4+1)
            else
                subplot(length(SimVal),3,(SimL-1)*3+1)
            end
            ylabel(sprintf('Nb of pairs in %s',SimLabel{SimVal(SimL)}))
        end
        Position=get(gcf,'position');
        Position(3:4)=[900,800];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        if TypeL==1
            saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
            close(h)
        end
    end
end

%% FIG31 - Limits versus number of networks


if ~isempty(find(FigRanks==31))
    FigRank=31;
    for LimitL=1:2
        if LimitL==1
            CurrLimit=LimitLim;
        else
            CurrLimit=LimitVal;
        end

        Lines={'-',':'};
        Markers='+o';

        h=figure;
        set(gcf,'color',[1,1,1])
        if LimitL==1
            set(h,'name',sprintf('FIG%ua - m%u: LimitLim: Limit vs network nb',FigRank,ChipRank))
        else
            set(h,'name',sprintf('FIG%ub - m%u: LimitVal: Limit vs network nb',FigRank,ChipRank))
        end
        PosIndex=1:NetNb;
        for TypeL=1:2
            for SimL=1:2
                subplot(2,4,(SimL-1)*4+1)
                hold on
                plot(LimitsNb{TypeL,SimL},sprintf('k%c',Lines{TypeL}))
                plot(LimitsNb{TypeL,SimL},sprintf('k%c',Markers(TypeL)))
                set(gca,'box','on')
                if TypeL==2&SimL==1
                    title('Nb of positive networks')
                end
                for ValL=1:3
                    subplot(2,4,(SimL-1)*4+1+ValL)
                    hold on
                    switch ValL
                        case 1
                            %PosIndex=find(CurrLimit.corr{TypeL,SimL}>-1);
                            plot(CurrLimit.corr{TypeL,SimL}(PosIndex),sprintf('r%c',Lines{TypeL}))
                            plot(CurrLimit.corr{TypeL,SimL}(PosIndex),sprintf('r%c',Markers(TypeL)))
                            if TypeL==1
                                Limit1=sprintf('%.0f',CurrLimit.corr{TypeL,SimL}(NetNb));
                            else
                                title(sprintf('CORR limits (%s,%.0f)',Limit1,CurrLimit.corr{TypeL,SimL}(NetNb)))
                            end
                        case 2
                            %PosIndex=find(CurrLimit.anti{TypeL,SimL}>-1);
                            plot(CurrLimit.anti{TypeL,SimL}(PosIndex),sprintf('b%c',Lines{TypeL}))
                            plot(CurrLimit.anti{TypeL,SimL}(PosIndex),sprintf('b%c',Markers(TypeL)))
                            if TypeL==2&SimL==1
                                title('ANTI limits')
                            end
                            if TypeL==1
                                Limit2=sprintf('%.0f',CurrLimit.anti{TypeL,SimL}(NetNb));
                            else
                                title(sprintf('ANTI limits (%s,%.0f)',Limit2,CurrLimit.anti{TypeL,SimL}(NetNb)))
                            end

                        case 3
                            %PosIndex=find(CurrLimit.pv{TypeL,SimL}<400);
                            plot(CurrLimit.pv{TypeL,SimL}(PosIndex),sprintf('g%c',Lines{TypeL}))
                            plot(CurrLimit.pv{TypeL,SimL}(PosIndex),sprintf('g%c',Markers(TypeL)))
                            if TypeL==2&SimL==1
                                title('PV limits')
                            end
                            if TypeL==1
                                Limit3=sprintf('%.0f',CurrLimit.pv{TypeL,SimL}(NetNb));
                            else
                                title(sprintf('PV limits (%s,%.0f)',Limit3,CurrLimit.pv{TypeL,SimL}(NetNb)))
                            end

                    end
                    set(gca,'box','on')
                end
                if TypeL==2
                    subplot(2,4,(SimL-1)*4+1)
                    if SimL==1
                        ylabel('InSim')
                    elseif SimL==2
                        ylabel('OutSim')
                    else
                        ylabel('RandSim')
                    end
                end
            end
        end
        Position=get(gcf,'position');
        Position(3:4)=[700,400];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        if LimitL==1
            saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
            close(h)
        end
    end
end


%% Recover limits
for TypeL=1:2
    for SimL=1:2
        LimitsLim.corr(TypeL,SimL)=LimitLim.corr{TypeL,SimL}(NetNb);
        LimitsLim.anti(TypeL,SimL)=LimitLim.anti{TypeL,SimL}(NetNb);
        LimitsLim.pv(TypeL,SimL)=LimitLim.pv{TypeL,SimL}(NetNb);
        LimitsVal.corr(TypeL,SimL)=LimitVal.corr{TypeL,SimL}(NetNb);
        LimitsVal.anti(TypeL,SimL)=LimitVal.anti{TypeL,SimL}(NetNb);
        LimitsVal.pv(TypeL,SimL)=LimitVal.pv{TypeL,SimL}(NetNb);
    end

end


'limlim corr'
LimitsLim.corr
mean(mean(LimitsLim.corr))
'limlim anti'
LimitsLim.anti
mean(mean(LimitsLim.anti))
'limlim pv'
LimitsLim.pv
mean(mean(LimitsLim.pv))

'limval corr'
LimitsVal.corr
mean(mean(LimitsVal.corr))
'limval anti'
LimitsVal.anti
mean(mean(LimitsVal.anti))
'limval pv'
LimitsVal.pv
mean(mean(LimitsVal.pv))

if ValFlag
    if MeanFlag
        %Use mean of four values
        Limit.corr=mean(LimitsVal.corr);
        Limit.anti=mean(LimitsVal.anti);
        Limit.pv=mean(LimitsVal.pv);
    else
        %Use only InSim values
        Limit.corr=LimitsVal.corr(1,1);
        Limit.anti=LimitsVal.anti(1,1);
        Limit.pv=LimitsVal.pv(1,1);
    end
else
    if MeanFlag
        %Use mean of four values
        Limit.corr=mean(LimitsLim.corr);
        Limit.anti=mean(LimitsLim.anti);
        Limit.pv=mean(LimitsLim.pv);
    else
        %Use only InSim values
        Limit.corr=LimitsLim.corr(1,1);
        Limit.anti=LimitsLim.anti(1,1);
        Limit.pv=LimitsLim.pv(1,1);
    end

end
DirName=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(DirName)
SaveFile=sprintf('limits_m%u_n%u_netnb%u_rank%u_v%u_m%u',ChipRank,FirstNetRanks(1),length(FirstNetRanks),PvCorrRank,ValFlag,MeanFlag);
eval(sprintf('save %s Limit',SaveFile))



if ~isempty(FigRanks)

    %for NbL=1:NetRankListNb
    for TypeL=1:2
        for SimL=1:2
            Corr{1}=[];
            for NetL=1:NetNb(1)
                Corr{1}=[Corr{1},SIM{1}{TypeL,SimL}{NetL}.corr];
            end
            PosIndex=find(sum(Corr{1},2));
            %STATISTICS ON PROBE SET COUPLE NOT IN SELECTION FILTER
            OutPos=find(MeanCorr{1}{TypeL}{SimL}+StdCorr{1}{TypeL}{SimL}<Limit.corr|...
                MeanAnti{1}{TypeL}{SimL}-StdAnti{1}{TypeL}{SimL}>Limit.anti|...
                MeanPv{1}{TypeL}{SimL}-StdPv{1}{TypeL}{SimL}>Limit.pv);

            %STATISTICS ON PROBE SET COUPLE NOT IN SELECTION FILTER
            %recover Couple informations
            ComGene{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.comGeneIn(PosIndex(OutPos));
            MeanComGene{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.meanComGeneIn(PosIndex(OutPos));
            UncomGene{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.uncomGeneIn(PosIndex(OutPos));
            ComTscript{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.comTscriptIn(PosIndex(OutPos));
            MeanComTscript{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.meanTscriptIn(PosIndex(OutPos));
            UncomTscript{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.uncomTscriptIn(PosIndex(OutPos));
            CMeanProbe{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.cMeanGroupProbeIn(PosIndex(OutPos));
            TMeanProbe{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.tMeanGroupProbeIn(PosIndex(OutPos));
            MaxProbeIn1{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.maxProbe1In(PosIndex(OutPos));
            MaxProbeIn2{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.maxProbe2In(PosIndex(OutPos));
            MinProbeIn1{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.minProbe1In(PosIndex(OutPos));
            MinProbeIn2{TypeL}{SimL}{2}=CoupleInfo{TypeL}{SimL}{1}.minProbe2In(PosIndex(OutPos));
        end
    end

end



%CALCUL PAR ETAPES SUCCESSIVES DES PARAMETRES STATISTIQUES
PairNb=zeros(2,3);
for TypeL=1:2
    for SimL=1:3
        PairNb(TypeL,SimL)=length(MeanPv{1}{TypeL}{SimL});
    end
end




%% FIG32

if ~isempty(find(FigRanks==32))
    FigRank=32;
    h=figure;
    set(gcf,'color',[1,1,1])
    if NetRankListNb==1
        SubNb1=1;
    else
        SubNb1=2;
    end
    SubNb2=3;
    set(h,'name',sprintf('FIG%u - m%u: Pv & Corr',FigRank,ChipRank))
    Colors='rmg';
    Lines={'-',':'};
    Symbols='o+';
    HLine=[];
    for TypeL=1:2
        for SimL=1:3
            %pv(overlap) distribution
            subplot(SubNb1,SubNb2,1)
            hold on
            Val=sort(MeanPv{1}{TypeL}{SimL});
            if isempty(Val)               
                l=plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));          
            else
                if length(Val)==1
                    l=plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
                else
                    l=plot(Val,0:1/(length(Val)-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}));
                    PvPos=find(Val<=Limit.pv);
                    if ~isempty(PvPos)
                        if length(PvPos)==1
                            plot(Val(PvPos(end)),1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        else
                            plot(Val(PvPos(end)),PvPos(end-1)*1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        end
                    end
                end
                HLine=[HLine,l];
            end
            % corr distribution
            subplot(SubNb1,SubNb2,2)
            hold on
            Val=sort(MeanCorr{1}{TypeL}{SimL});
            if isempty(Val)
                plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
            else
                if length(Val)==1
                    plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                else
                    plot(Val,0:1/(length(Val)-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}))
                    CorrPos=find(Val>=Limit.corr);
                    if ~isempty(CorrPos)
                        if CorrPos==1
                            plot(Val(CorrPos(1)),1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        else
                            plot(Val(CorrPos(1)),(CorrPos(1)-1)*1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        end
                    end
                end
            end

            % anti distribution
            subplot(SubNb1,SubNb2,3)
            hold on
            Val=sort(MeanAnti{1}{TypeL}{SimL});
            if isempty(Val)
                plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
            else
                if length(Val)==1
                    plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                else
                    plot(Val,0:1/(length(Val)-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}))
                    AntiPos=find(Val<=Limit.anti);
                    if ~isempty(AntiPos)
                        if length(AntiPos)==1
                            plot(Val(AntiPos(end)),1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        else
                            plot(Val(AntiPos(end)),AntiPos(end-1)*1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        end
                    end

                end
            end

            %no-qlimit distribution
            if NetRankListNb==2
                subplot(SubNb1,SubNb2,4)
                hold on
                Val=sort(OutMeanCorr{2}{TypeL}{SimL});
                NanPos=find(isnan(Val));
                Val(NanPos)=[];
                if isempty(Val)
                    plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
                else
                    if length(Val)==1
                        plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                    else
                        plot(Val,0:1/(length(Val)-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}))
                    end
                end
                HLine=[HLine,l];
                subplot(SubNb1,SubNb2,5)
                hold on
                plot(MeanCorr{1}{TypeL}{SimL}+(SimL-1)*100,OutMeanCorr{2}{TypeL}{SimL}+(SimL-1)*100,sprintf('%c%c',Colors(SimL),Symbols(TypeL)),'markersize',3)
                subplot(SubNb1,SubNb2,6)
                hold on
                plot(MeanAnti{1}{TypeL}{SimL}+(SimL-1)*100,OutMeanAnti{2}{TypeL}{SimL}+(SimL-1)*100,sprintf('%c%c',Colors(SimL),Symbols(TypeL)),'markersize',3)
            end
        end
    end


    Colors='kbc';
    for TypeL=1:2
        for SimL=1:3
            subplot(SubNb1,SubNb2,1)
            hold on
            Val=sort(OutMeanPv{1}{TypeL}{SimL});
            PosIndex=find(~isnan(Val));
            Val=Val(PosIndex);
            if isempty(Val)
                l=plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
            else
                if length(Val)==1
                    l=plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
                else
                    l=plot(Val,0:1/(length(Val)-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}));
                    PvPos=find(Val<=Limit.pv);
                    if ~isempty(PvPos)
                        if length(PvPos)==1
                            plot(Val(PvPos(end)),1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        else
                            plot(Val(PvPos(end)),PvPos(end-1)*1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        end
                    end
                end
            end
            HLine=[HLine,l];
        end
    end

    Colors='rmg';
    Lines={'--','-.'};
    %Symbols='sp';
    for TypeL=1:2
        for SimL=1:3
            subplot(SubNb1,SubNb2,1)
            hold on
            Val=sort(AllNullMPv{1}{TypeL}{SimL});
            if isempty(Val)
                l=plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
            else
                if length(Val)==1
                    l=plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
                else
                    l=plot(Val,0:1/(length(AllNullMPv{1}{TypeL}{SimL})-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}));
                    PvPos=find(Val<=Limit.pv);
                    if ~isempty(PvPos)
                        if length(PvPos)==1
                            plot(Val(PvPos(end)),1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        else
                            plot(Val(PvPos(end)),PvPos(end-1)*1/(length(Val)-1),sprintf('%c%c',Colors(SimL),Symbols(TypeL)))
                        end
                    end
                end
            end
            HLine=[HLine,l];
            if NetRankListNb==2
                subplot(SubNb1,SubNb2,4)
                hold on
                Val=sort(AllNullMCorr{2}{TypeL}{SimL});
                if isempty(Val)
                    l=plot(0,0,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
                else
                    if length(Val)==1
                        l=plot(Val,1,sprintf('%c%c',Colors(SimL),Symbols(TypeL)));
                    else
                        l=plot(Val,0:1/(length(Val)-1):1,sprintf('%c%c',Colors(SimL),Lines{TypeL}));
                    end
                end
                %HLine=[HLine,l];
            end
        end
    end

    legend(HLine,'sIn IsCorr Corr>0','sOut IsCorr Corr>0','sRand IsCorr Corr>0','mIn IsCorr Corr>0','mOut IsCorr Corr>0','mRand IsCorr Corr>0',...
        'sIn IsCorr Corr=0','sOut IsCorr Corr=0','sRand IsCorr Corr=0','mIn IsCorr Corr=0','mOut IsCorr Corr=0','mRand IsCorr Corr=0',...
        'sIn NoCorr','sOut NoCorr','sRand NoCorr','mIn NoCorr','mOut NoCorr','mRand NoCorr',...
        'location','NorthWest')


    subplot(SubNb1,SubNb2,1)
    set(gca,'box','on')
    title('PV distribution')
    xlabel('PV')
    ylabel('cdf')


    subplot(SubNb1,SubNb2,2)
    set(gca,'box','on')
    title('CORR distribution')
    xlabel('CORR')
    ylabel('cdf')

    subplot(SubNb1,SubNb2,3)
    set(gca,'box','on')
    title('ANTI distribution')
    xlabel('ANTI')
    ylabel('cdf')



    if NetRankListNb==2
        subplot(SubNb1,SubNb2,4)
        set(gca,'box','on')
        title('CORR distribution - no qlimit')
        xlabel('CORR')
        ylabel('cdf')

        subplot(SubNb1,SubNb2,5)
        line([0,300],[0,300])
        set(gca,'box','on')
        title('CORR no qlimit vs qlimit plot')
        xlabel('CORR - qlimit')
        ylabel('CORR - no qlimit')

        subplot(SubNb1,SubNb2,6)
        line([0,300],[0,300])
        set(gca,'box','on')
        title('ANTI no qlimit vs qlimit plot')
        xlabel('ANTI - qlimit')
        ylabel('ANTI - no qlimit')

    end
    set(gcf,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[1000,700];
    set(gcf,'position',Position)
    try
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    catch
        mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    end
    saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
    close(h)
end



%% FIG33
if ~isempty(find(FigRanks==33))
    FigRank=33;
    %Without 0 class
    Legend={'sIn IsCorr Corr>0','sOut IsCorr Corr>0','sRand IsCorr Corr>0','mIn IsCorr Corr>0',...
        'mOut IsCorr Corr>0','mRand IsCorr Corr>0'};
    h=figure;
    set(h,'name',sprintf('FIG%u - m%u: histogram number of positive networks',FigRank,ChipRank))
    set(gcf,'color',[1,1,1])
    PosIndex=0;
    for TypeL=1:2
        for SimL=1:3
            PosIndex=PosIndex+1;
            subplot(2,3,PosIndex)
            hold on
            hist(RepNb{1}{TypeL}{SimL},NetNb(1));
            set(gca,'xlim',[1,NetNb(1)])
            set(gca,'box','on')
            xlabel('rep nb')
            ylabel('freq')
            title(Legend{PosIndex})
        end
    end
    Position=get(gcf,'position');
    Position(3:4)=[600,400];
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


%% FIG34
if ~isempty(find(FigRanks==34))
    FigRank=34;
    %voisinage des couples
    %charger toutes les valeurs du réseau (CORR et ANTI)
    Corr=[];
    for NetL=1:NetNb(1)
        Corr=[Corr,SIM{1}{1,1}{NetL}.corr];
    end
    %utiliser des couples dont corr est significatif dans la moitié des réseaux
    PosIndex=find(sum(Corr,2));
    RepPos=find(RepNb{1}{1}{1}==round(NetNb(1)/2));
    RepPos=PosIndex(RepPos);
    %index des réseaux pour lesquels corr est significatif
    CoupleNb=min(10,length(RepPos));
    PosIndex=zeros(CoupleNb,round(NetNb(1)/2));
    for PosL=1:CoupleNb
        PosIndex(PosL,:)=find(Corr(RepPos(PosL),:));
    end

    for ListL=1:NetRankListNb
        for CoupleL=1:CoupleNb
            AllCorr{ListL}{CoupleL}{1}=uint8(zeros(PsNb,NetNb(1)));
            AllCorr{ListL}{CoupleL}{2}=uint8(zeros(PsNb,NetNb(1)));
            AllAnti{ListL}{CoupleL}{1}=uint8(zeros(PsNb,NetNb(1)));
            AllAnti{ListL}{CoupleL}{2}=uint8(zeros(PsNb,NetNb(1)));
        end
    end

    for ListL=1:NetRankListNb
        for NetL=1:NetNb(1)
            %load current net values
            CurrNetRank=NetRanks{ListL}(NetL);
            CNetFile=sprintf('c_m%u_n%u.4mat',ChipRank,CurrNetRank);
            ANetFile=sprintf('a_m%u_n%u.4mat',ChipRank,CurrNetRank);
            NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',CurrNetRank));
            cd(NetDir);
            for CoupleL=1:CoupleNb
                %recover corr and anti values for each probe set of the pair
                PsRank1=CoupleInfo{1}{1}{1}.psRank1(RepPos(CoupleL));
                PsRank2=CoupleInfo{1}{1}{1}.psRank2(RepPos(CoupleL));              
                AllCorr{ListL}{CoupleL}{1}(:,NetL)=load_data(CNetFile,NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsRank1);
                AllCorr{ListL}{CoupleL}{1}(PsRank1,NetL)=0;
                AllCorr{ListL}{CoupleL}{2}(:,NetL)=load_data(CNetFile,NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsRank2);
                AllCorr{ListL}{CoupleL}{1}(PsRank2,NetL)=0;
                AllAnti{ListL}{CoupleL}{1}(:,NetL)=load_data(ANetFile,NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsRank1);
                AllAnti{ListL}{CoupleL}{2}(:,NetL)=load_data(ANetFile,NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsRank2);
            end
        end
    end


    %FIG34a - FIG34b
    %CORR
    % reproducibility of correlation values of a probe set in different networks
    Colors=colors(colormap,round(NetNb(1)/2)*(round(NetNb(1)/2)-1)/2);
    CorrCoeff=cell(2,1);
    Common=cell(2,1);
    for FigL=1:2

        h=figure;
        if FigL==1
            set(h,'name',sprintf('FIG%ua - m%u: Corr between Corr>0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        else
            set(h,'name',sprintf('FIG%ub - m%u: Corr between Corr==0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        end

        for CoupleL=1:CoupleNb
            PsRank1=CoupleInfo{1}{1}{1}.psRank1(RepPos(CoupleL));
            PsRank2=CoupleInfo{1}{1}{1}.psRank2(RepPos(CoupleL));
            Pc{CoupleL}{FigL}=[];
            Comm{CoupleL}{FigL}=[];
            CurrIndex=PosIndex(CoupleL,:);
            if FigL==2
                CurrIndex=setdiff([1:NetNb(1)],CurrIndex);
            end
            if ~isempty(FigRanks)
                subplot(2,5,CoupleL)
                hold on
                CompRank=0;
                for NetL1=1:length(CurrIndex)-1
                    for NetL2=NetL1+1:length(CurrIndex)
                        CompRank=CompRank+1;
                        plot(AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL1)),AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL2)),'+','color',Colors(CompRank,:),'markersize',3)
                    end
                end
                set(gca,'box','on')
                title(sprintf('ps %u of couple %u',PsRank1,CoupleL))
                xlabel('corr in first network')
                ylabel('corr in second network')
            end
            for NetL1=1:length(CurrIndex)-1
                for NetL2=NetL1+1:length(CurrIndex);
                    ComPos=find(AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL1))&AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL2)));
                    CurrCorrCoeff=corrcoef(single(AllCorr{1}{CoupleL}{1}(ComPos,CurrIndex(NetL1))),...
                        single(AllCorr{1}{CoupleL}{1}(ComPos,CurrIndex(NetL2))));
                    CurrCommon=length(ComPos)*100/sqrt(length(find(AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL1))))*length(find(AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL2)))));
                    Comm{CoupleL}{FigL}=[Comm{CoupleL}{FigL},CurrCommon];
                    if length(CurrCorrCoeff)==2
                        Pc{CoupleL}{FigL}=[Pc{CoupleL}{FigL};CurrCorrCoeff(2)];
                    else
                        Pc{CoupleL}{FigL}=[Pc{CoupleL}{FigL};CurrCorrCoeff];
                    end
                end
            end
            CorrCoeff{1}=zeros(length(CurrIndex),2);
            Common{1}=zeros(length(CurrIndex),2);
            CorrCoeff{2}=zeros(length(CurrIndex),2);
            Common{2}=zeros(length(CurrIndex),2);

            CorrCoeff{FigL}(CoupleL,1)=mean(Pc{CoupleL}{FigL});
            CorrCoeff{FigL}(CoupleL,2)=std(Pc{CoupleL}{FigL});
            Common{FigL}(CoupleL,1)=mean(Comm{CoupleL}{FigL});
            Common{FigL}(CoupleL,2)=std(Comm{CoupleL}{FigL});
        end
        Position=get(gcf,'position');
        Position(3:4)=[1100,450];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        if FigL==1
            saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
            close(h)
        end
    end


    %FIG34c - FIG34d
    %reproductibilité entre network pour un probeset donné
    %ANTI
    for FigL=1:2
        h=figure;
        if FigL==1
            set(h,'name',sprintf('FIG%uc - m%u: Anti between Corr>0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        else
            set(h,'name',sprintf('FIG%ud - m%u: Anti between Corr==0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        end
        for CoupleL=1:CoupleNb
            PsRank1=CoupleInfo{1}{1}{1}.psRank1(RepPos(CoupleL));
            PsRank2=CoupleInfo{1}{1}{1}.psRank2(RepPos(CoupleL));
            CurrIndex=PosIndex(CoupleL,:);
            if FigL==2
                CurrIndex=setdiff([1:NetNb(1)],CurrIndex);
            end
            if ~isempty(FigRanks)
                subplot(2,5,CoupleL)
                hold on
                CompRank=0;
                for NetL1=1:length(CurrIndex)-1
                    for NetL2=NetL1+1:length(CurrIndex)
                        CompRank=CompRank+1;
                        plot(AllAnti{1}{CoupleL}{1}(:,CurrIndex(NetL1)),AllAnti{1}{CoupleL}{1}(:,CurrIndex(NetL2)),'+','color',Colors(CompRank,:),'markersize',3)
                    end
                end
                set(gca,'box','on')
                title(sprintf('ps %u of couple %u',PsRank1,CoupleL))
                xlabel('anti in first network')
                ylabel('anti in second network')
            end
        end
        Position=get(gcf,'position');
        Position(3:4)=[1100,450];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        if FigL==1
            saveas(h,sprintf('m%u_fig%uc.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%ud.png',ChipRank,FigRank),'png')
            close(h)
        end
    end


    %FIG34e - FIG34f
    %reproductibilité des CORR entre probe sets d'un même couple à l'intérieur d'un même
    %réseau
    Colors=colors(colormap,round(NetNb/2));
    for FigL=1:2
        h=figure;
        if FigL==1
            set(h,'name',sprintf('FIG%ue - m%u: Corr in Corr>0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        else
            set(h,'name',sprintf('FIG%uf - m%u: Corr in Corr==0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        end
        for CoupleL=1:CoupleNb
            PsRank1=CoupleInfo{1}{1}{1}.psRank1(RepPos(CoupleL));
            PsRank2=CoupleInfo{1}{1}{1}.psRank2(RepPos(CoupleL));
            CurrIndex=PosIndex(CoupleL,:);
            if FigL==2
                CurrIndex=setdiff([1:NetNb(1)],CurrIndex);
            end

            subplot(2,5,CoupleL)
            hold on
            CompRank=0;
            for NetL=1:length(CurrIndex)
                CompRank=CompRank+1;
                plot(AllCorr{1}{CoupleL}{1}(:,CurrIndex(NetL)),AllCorr{1}{CoupleL}{2}(:,CurrIndex(NetL)),'+','color',Colors(CompRank,:),'markersize',3)
            end
            set(gca,'box','on')
            title(sprintf('couple %u (Ps %u and %u)',CoupleL,PsRank1,PsRank2))
            xlabel('corr of first ps')
            ylabel('corr of second ps')
        end
        Position=get(gcf,'position');
        Position(3:4)=[1100,450];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        if FigL==1
            saveas(h,sprintf('m%u_fig%ue.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%uf.png',ChipRank,FigRank),'png')
            close(h)
        end
    end



    %FIG34g - FIG34h
    %reproductibilité des ANTI entre probe sets d'un même couple à l'intérieur d'un même
    %réseau
    Colors=colors(colormap,round(NetNb/2));
    for FigL=1:2
        h=figure;
        if FigL==1
            set(h,'name',sprintf('FIG%ug - m%u: Anti in Corr>0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        else
            set(h,'name',sprintf('FIG%uh - m%u: Anti in Corr==0 networks',FigRank,ChipRank))
            set(gcf,'color',[1,1,1])
        end
        for CoupleL=1:CoupleNb
            PsRank1=CoupleInfo{1}{1}{1}.psRank1(RepPos(CoupleL));
            PsRank2=CoupleInfo{1}{1}{1}.psRank2(RepPos(CoupleL));
            CurrIndex=PosIndex(CoupleL,:);
            if FigL==2
                CurrIndex=setdiff([1:NetNb(1)],CurrIndex);
            end

            subplot(2,5,CoupleL)
            hold on
            CompRank=0;
            for NetL=1:length(CurrIndex)
                CompRank=CompRank+1;
                plot(AllAnti{1}{CoupleL}{1}(:,CurrIndex(NetL)),AllAnti{1}{CoupleL}{2}(:,CurrIndex(NetL)),'+','color',Colors(CompRank,:),'markersize',3)
            end
            set(gca,'box','on')
            title(sprintf('couple %u (Ps %u and %u)',CoupleL,PsRank1,PsRank2))
            xlabel('anti of first ps')
            ylabel('anti of second ps')
        end
        Position=get(gcf,'position');
        Position(3:4)=[1100,450];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        if FigL==1
            saveas(h,sprintf('m%u_fig%ug.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%uh.png',ChipRank,FigRank),'png')
            close(h)
        end
    end




    %FIG34i
    % distribution du voisinage
    Conn{1}=zeros(CoupleNb,2,round(NetNb(1)/2));
    Conn{2}=zeros(CoupleNb,2,NetNb(1)-round(NetNb(1)/2));
    for CoupleL=1:CoupleNb
        NetIndex=PosIndex(CoupleL,:);
        for NetL=1:length(NetIndex)
            Conn{1}(CoupleL,1,NetL)=length(find(AllCorr{1}{CoupleL}{1}(:,NetIndex(NetL))));
            Conn{1}(CoupleL,2,NetL)=length(find(AllCorr{1}{CoupleL}{2}(:,NetIndex(NetL))));
        end
        NetIndex=setdiff([1:NetNb(1)],NetIndex);
        for NetL=1:length(NetIndex)
            Conn{2}(CoupleL,1,NetL)=length(find(AllCorr{1}{CoupleL}{1}(:,NetIndex(NetL))));
            Conn{2}(CoupleL,2,NetL)=length(find(AllCorr{1}{CoupleL}{2}(:,NetIndex(NetL))));
        end
    end
    MeanConn=zeros(CoupleNb,4);
    StdConn=zeros(CoupleNb,4);
    for CoupleL=1:CoupleNb
        MeanConn(CoupleL,1)=mean(Conn{1}(CoupleL,1,:));
        MeanConn(CoupleL,2)=mean(Conn{1}(CoupleL,2,:));
        MeanConn(CoupleL,3)=mean(Conn{2}(CoupleL,1,:));
        MeanConn(CoupleL,4)=mean(Conn{2}(CoupleL,2,:));
        StdConn(CoupleL,1)=std(Conn{1}(CoupleL,1,:));
        StdConn(CoupleL,2)=std(Conn{1}(CoupleL,2,:));
        StdConn(CoupleL,3)=std(Conn{2}(CoupleL,1,:));
        StdConn(CoupleL,4)=std(Conn{2}(CoupleL,2,:));

    end

    h=figure;
    set(h,'name',sprintf('FIG%ui - m%u: CONNECTIVITY',FigRank,ChipRank))
    subplot(1,2,1)
    Colors=colors(colormap,CoupleNb);
    plot(MeanConn(:,1),MeanConn(:,3),'r+')
    hold on
    plot(MeanConn(:,2),MeanConn(:,4),'ro')
    for CoupleL=1:CoupleNb
        line([MeanConn(CoupleL,1),MeanConn(CoupleL,2)],[MeanConn(CoupleL,3),MeanConn(CoupleL,4)],'color',Colors(CoupleL,:))
    end
    x=get(gca,'xlim');
    y=get(gca,'ylim');
    line([x(1),x(2)],[y(1),y(2)])
    xlabel(sprintf('mean connectivity in networks\nwhere corr>0'))
    ylabel(sprintf('mean connectivity in networks\nwhere corr=0'))
    set(gcf,'color',[1,1,1])

    subplot(1,2,2)
    plot(StdConn(:,1),StdConn(:,3),'r+')
    hold on
    plot(StdConn(:,2),StdConn(:,4),'ro')
    for CoupleL=1:CoupleNb
        line([StdConn(CoupleL,1),StdConn(CoupleL,2)],[StdConn(CoupleL,3),StdConn(CoupleL,4)],'color',Colors(CoupleL,:))
    end
    x=get(gca,'xlim');
    y=get(gca,'ylim');
    line([x(1),x(2)],[y(1),y(2)])
    xlabel(sprintf('std connectivity in networks\nwhere corr>0'))
    ylabel(sprintf('std connectivity in networks\nwhere corr=0'))
    set(gcf,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[1000,700];
    set(gcf,'position',Position)
    try
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    catch
        mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    end
    saveas(h,sprintf('m%u_fig%ui.png',ChipRank,FigRank),'png')
    close(h)



    %FIG34j
    %la diminution du nombre de liens est dûe à une disparition des liens présents dans les autres réseaux
    % Moyenne des connectiviés - chaque couple de probe sets est relié par un trait.
    % Pour chaque probe set ont calcule la connectivité dans chaque réseau et l'on fait les moyennes des valeurs
    % dans les réseaux pour lesquels les corr sont >0 et dans les réseaux pour lesquels corr=0.
    % La connectivité est systématiquement plus faible dans les réseaux pour lesquels corr=0.
    % Conforte l'hypothèse que dans ces réseaux où la corrélation du couple a été déclarée non significatives, une grande partie des corrélations avec
    % les autres probe sets qui sont reliés normalement avec les probe sets du couples ont été également déclarées non significatives.
    % Il s'agit sans doute d'un effet du choix des conditions biologiques qui ne permettent pas de mettre en évidence de nombreuses corrélations entre les probe
    % sets du couples et ses partenaires 'naturels'.
    % On le confirme ainsi :
    % Si on construit des réseaux union dans chaque catégorie (corr=0 et corr>0), on observe que le % de
    % liens spécifique au réseau union-corr=0 est faible  alors que le % de liens spécifique à union-corr>0 est plus fort.
    InterSpec1=zeros(CoupleNb,2);
    UnionSpec1=zeros(CoupleNb,2);
    InterSpec2=zeros(CoupleNb,2);
    UnionSpec2=zeros(CoupleNb,2);
    for CoupleL=1:CoupleNb
        % Values in networks where corr>0
        NetIndex=PosIndex(CoupleL,:);
        InterNet1=ones(PsNb,1);
        UnionNet1=zeros(PsNb,1);
        InterNet2=ones(PsNb,1);
        UnionNet2=zeros(PsNb,1);
        for NetL=1:length(NetIndex)
            UnionNet1(find(AllCorr{1}{CoupleL}{1}(:,NetIndex(NetL))))=1;
            InterNet1=InterNet1&AllCorr{1}{CoupleL}{1}(:,NetIndex(NetL))>0;
            UnionNet2(find(AllCorr{1}{CoupleL}{2}(:,NetIndex(NetL))))=1;
            InterNet2=InterNet2&AllCorr{1}{CoupleL}{2}(:,NetIndex(NetL))>0;
        end
        % Values in networks where corr=0
        NetIndex=setdiff([1:NetNb(1)],NetIndex);
        InterNet3=ones(PsNb,1);
        UnionNet3=zeros(PsNb,1);
        InterNet4=ones(PsNb,1);
        UnionNet4=zeros(PsNb,1);
        for NetL=1:length(NetIndex)
            UnionNet3(find(AllCorr{1}{CoupleL}{1}(:,NetIndex(NetL))))=1;
            InterNet3=InterNet3&AllCorr{1}{CoupleL}{1}(:,NetIndex(NetL))>0;
            UnionNet4(find(AllCorr{1}{CoupleL}{2}(:,NetIndex(NetL))))=1;
            InterNet4=InterNet4&AllCorr{1}{CoupleL}{2}(:,NetIndex(NetL))>0;
        end
        % % of neighbors that are specific to corr=0
        %prevent NaN
        Spec=setdiff(find(InterNet3),find(InterNet1));
        InterSpec1(CoupleL,1)=length(Spec)*100/(max(1,length(find(InterNet3))));
        Spec=setdiff(find(InterNet4),find(InterNet2));
        InterSpec1(CoupleL,2)=length(Spec)*100/max(1,length(find(InterNet4)));
        Spec=setdiff(find(UnionNet3),find(UnionNet1));
        UnionSpec1(CoupleL,1)=length(Spec)*100/max(1,length(find(UnionNet3)));
        Spec=setdiff(find(UnionNet4),find(UnionNet2));
        UnionSpec1(CoupleL,2)=length(Spec)*100/max(1,length(find(UnionNet4)));

        % % of neighbors that are specific to corr>0
        Spec=setdiff(find(InterNet1),find(InterNet3));
        InterSpec2(CoupleL,1)=length(Spec)*100/max(1,length(find(InterNet1)));
        Spec=setdiff(find(InterNet2),find(InterNet4));
        InterSpec2(CoupleL,2)=length(Spec)*100/max(1,length(find(InterNet2)));
        Spec=setdiff(find(UnionNet1),find(UnionNet3));
        UnionSpec2(CoupleL,1)=length(Spec)*100/max(1,length(find(UnionNet1)));
        Spec=setdiff(find(UnionNet2),find(UnionNet4));
        UnionSpec2(CoupleL,2)=length(Spec)*100/max(1,length(find(UnionNet2)));
    end


    h=figure;
    set(h,'name',sprintf('FIG%uj - m%u: SPECIFICITY',FigRank,ChipRank))
    set(gcf,'color',[1,1,1])
    Colors=colors(colormap,CoupleNb);
    subplot(1,2,1)
    % Intersection - Specific to corr==0
    plot(InterSpec1(:),'b')
    hold on
    plot(InterSpec1(:),'b+')
    title('intersection of networks')
    xlabel(' ps pairs')
    ylabel('% of specificity')


    % Union - Specific to corr==0
    subplot(1,2,2)
    plot(UnionSpec1(:),'c-')
    hold on
    plot(UnionSpec1(:),'c+')
    title('union of networks')
    xlabel(' ps pairs')
    ylabel('% of specificity')


    % Intersection - Specific to corr>0
    subplot(1,2,1)
    plot(InterSpec2(:),'r')
    plot(InterSpec2(:),'ro')

    % Union - Specific to corr>0
    subplot(1,2,2)
    plot(UnionSpec2(:),'m')
    plot(UnionSpec2(:),'mo')

    Position=get(gcf,'position');
    Position(3:4)=[1000,700];
    set(gcf,'position',Position)
    try
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    catch
        mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    end
    saveas(h,sprintf('m%u_fig%uj.png',ChipRank,FigRank),'png')
    close(h)
end

%% FIG35
%Caracteristics of probes

if ~isempty(find(FigRanks==35))
    FigRank=35;

    Colors='rbk';
    Lines={'-','-.'};
    for TypeL=1:2
        h=figure;
        if TypeL==1
            set(h,'name',sprintf('FIG%ua - m%u: ps properties for single',FigRank,ChipRank))
        else
            set(h,'name',sprintf('FIG%ub - m%u: ps properties for multiple',FigRank,ChipRank))
        end
        for SimL=1:2
            subplot(3,4,1)
            hold on
            for GrpL=1:3
                Val=ComGene{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('common genes')
            set(gca,'box','on')
        end
        for SimL=1:2
            subplot(3,4,2)
            hold on
            for GrpL=1:3
                Val=MeanComGene{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('mean common genes')
            set(gca,'box','on')
            set(gca,'ylim',[0,0.2])
        end
        for SimL=1:2
            subplot(3,4,3)
            hold on
            for GrpL=1:3

                Val=UncomGene{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('uncommon genes')
            set(gca,'box','on')
            set(gca,'xlim',[0,5])
        end
        for SimL=1:2
            subplot(3,4,4)
            hold on
            for GrpL=1:3

                Val=ComTscript{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('common transcripts')
            set(gca,'box','on')
            set(gca,'xlim',[0,5])
        end
        for SimL=1:2
            subplot(3,4,5)
            hold on
            for GrpL=1:3

                Val=MeanComTscript{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('mean common transcripts')
            set(gca,'box','on')
            set(gca,'xlim',[0,1])
        end
        for SimL=1:2
            subplot(3,4,6)
            hold on
            for GrpL=1:3

                Val=UncomTscript{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('uncommon transcripts')
            set(gca,'box','on')
        end
        for SimL=1:2
            subplot(3,4,7)
            hold on
            for GrpL=1:3

                Val=CMeanProbe{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('c mean probe')
            set(gca,'box','on')
        end
        for SimL=1:2
            subplot(3,4,8)
            hold on
            for GrpL=1:3

                Val=TMeanProbe{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('t mean probe')
            set(gca,'box','on')
        end
        for SimL=1:2
            subplot(3,4,9)
            hold on
            for GrpL=1:3

                Val=MaxProbeIn1{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('max probe in 1')
            set(gca,'box','on')
            set(gca,'xlim',[0,ProbeNb+1])
        end
        for SimL=1:2
            subplot(3,4,10)
            hold on
            for GrpL=1:3

                Val=MaxProbeIn2{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('max probe in 2')
            set(gca,'box','on')
            set(gca,'xlim',[0,ProbeNb+1])
        end
        for SimL=1:2
            subplot(3,4,11)
            hold on
            for GrpL=1:3

                Val=MinProbeIn1{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('min probe in 1')
            set(gca,'box','on')
            set(gca,'xlim',[0,ProbeNb+1])
        end
        for SimL=1:2
            subplot(3,4,12)
            hold on
            for GrpL=1:3
                Val=MinProbeIn2{TypeL}{SimL}{GrpL};
                plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(GrpL),Lines{SimL}))
            end
            title('min probe in 2')
            set(gca,'box','on')
            set(gca,'xlim',[0,ProbeNb+1])
        end
        legend('InSim Good','InSim Bad','InSim Corr=0','OutSim Good','OutSim Bad','OutSim Corr=0')
        set(gcf,'color',[1,1,1])
        Position=get(gcf,'position');
        Position(3:4)=[1000,800];

        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u_%up',ChipRank),'png'))
        end
        if TypeL==1
            saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
            close(h)
        else
            saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
            close(h)
        end
    end
end


%% FIG36
if ~isempty(find(FigRanks==36))
    FigRank=36;
    %Statistics on all significatives Corr,Anti and Pv values


    h(1)=figure;
    set(h(1),'name',sprintf('FIG%ua - m%u: distribution of min,max,mean,std of CORR on pairs significative in all networks',FigRank,ChipRank))
    h(2)=figure;
    set(h(2),'name',sprintf('FIG%ub - m%u: distribution of min,max,mean,std of ANTI on pairs significative in all networks',FigRank,ChipRank))
    h(3)=figure;
    set(h(3),'name',sprintf('FIG%uc - m%u: distribution of min,max,mean,std of PV on pairs significative in all networks',FigRank,ChipRank))

    for TypeL=1:2
        if TypeL==1
            Title1='single';
        else
            Title1='multiple';
        end
        for SimL=1:2
            if SimL==1
                Title=['InExons - ',Title1];
            else
                Title=['OutOfExons - ',Title1];
            end
            PosIndex=find(RepNb{1}{TypeL}{SimL}==NetNb);
            for ValL=1:3
                figure(h(ValL))
                subplot(2,2,(TypeL-1)*2+SimL)
                title(Title)
                hold on
                switch ValL
                    case 1
                        CurrVal=Corrs{TypeL}{SimL};
                    case 2
                        CurrVal=Antis{TypeL}{SimL};
                    case 3
                        CurrVal=Pvs{TypeL}{SimL};
                end
                CurrMin=zeros(length(PosIndex),1);
                CurrMax=zeros(length(PosIndex),1);
                CurrMean=zeros(length(PosIndex),1);
                CurrStd=zeros(length(PosIndex),1);
                CurrLimit=zeros(length(PosIndex),1);
                for PosL=1:length(PosIndex)
                    CurrMin(PosL)=min(CurrVal{PosIndex(PosL)});
                    CurrMax(PosL)=max(CurrVal{PosIndex(PosL)});
                    CurrMean(PosL)=mean(CurrVal{PosIndex(PosL)});
                    CurrStd(PosL)=std(single(CurrVal{PosIndex(PosL)}));
                    Val=sort(CurrVal{PosIndex(PosL)});
                    if ValL==1
                        CurrLimit(PosL)=Val(round(length(Val)*0.95));
                        %CurrLimit(PosL)=mean(CurrVal{PosIndex(PosL)})-std(single(CurrVal{PosIndex(PosL)}));
                    else
                        CurrLimit(PosL)=Val(round(length(Val)*0.05));
                        %CurrLimit(PosL)=mean(CurrVal{PosIndex(PosL)})+std(single(CurrVal{PosIndex(PosL)}));
                    end
                end
                [CurrLimit,SortIndex]=sort(CurrLimit);
                CurrMin=CurrMin(SortIndex);
                CurrMax=CurrMax(SortIndex);
                CurrMean=CurrMean(SortIndex);
                CurrStd=CurrStd(SortIndex);

                plot(1:length(PosIndex),CurrMin,'b.','markersize',3)
                plot(1:length(PosIndex),CurrMax,'r.','markersize',3)
                plot(1:length(PosIndex),CurrMean,'y.','markersize',3)
                plot(1:length(PosIndex),CurrStd,'k.','markersize',3)
                plot(1:length(PosIndex),CurrLimit,'g.','markersize',3)
                set(gca,'box','on')

                %                     switch ValL
                %                         case 1
                %                             Limits.corr(TypeL,SimL)=CurrLimit(ceil(length(PosIndex)*0.05));
                %                         case 2
                %                             Limits.anti(TypeL,SimL)=CurrLimit(ceil(length(PosIndex)*0.95));
                %                         case 3
                %                             Limits.pv(TypeL,SimL)=CurrLimit(ceil(length(PosIndex)*0.95));
                %                     end
            end
        end
    end
    for FigL=1:3
        set(h(FigL),'color',[1,1,1])
        Position=get(gcf,'position');
        Position(3:4)=[600,450];
        figure(h(FigL))
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
        saveas(h(FigL),sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(FigL)),'png')
        close(h(FigL))
    end
end



%% FIG37
if ~isempty(find(FigRanks==37))
    FigRank=37;
    %Relation between Corr,Anti and Pv values
    h=figure;
    set(h,'name',sprintf('FIG%ub - m%u: In Out Stats',FigRank,ChipRank))
    set(gcf,'color',[1,1,1])
    GoodNb=zeros(2,2);
    BadNb=zeros(2,2);
    SelNb=cell(2,2);
    for i=1:2
        for j=1:2
            SelNb{i,j}=zeros(NetNb+1,1);
        end
    end
    for TypeL=1:2
        if TypeL==1
            Title1='single';
        else
            Title1='multiple';
        end
        for SimL=1:2

            subplot(2,2,(TypeL-1)*2+SimL)
            if SimL==1
                Title=['InExons - ',Title1];
            else
                Title=['OutOfExons - ',Title1];
            end
            title(Title)
            plot(0:NetNb,SelNb{TypeL,SimL},'r')
            hold on
            plot(0:NetNb,SelNb{TypeL,SimL},'r+')
            xlabel('Net nb')
            ylabel('frequency')
            title(Title)
        end
        Position=get(gcf,'position');
        Position(3:4)=[600,450];
        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
    end
    saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
    close(h)

end


%% FIG38
if ~isempty(find(FigRanks==38))
    FigRank=38;
    %Relation between Corr,Anti and Pv values
    h=figure;
    set(h,'name',sprintf('FIG%u - m%u: CORR,ANTI,PV',FigRank,ChipRank))
    set(gcf,'color',[1,1,1])
    GoodNb=zeros(2,2);
    BadNb=zeros(2,2);
    SelNb=cell(2,2);
    for i=1:2
        for j=1:2
            SelNb{i,j}=zeros(NetNb+1,1);
        end
    end
    for TypeL=1:2
        if TypeL==1
            Title1='single';
        else
            Title1='multiple';
        end
        for SimL=1:2
            subplot(2,2,(TypeL-1)*2+SimL)
            hold on
            if SimL==1
                Title=['InExons - ',Title1];
            else
                Title=['OutOfExons - ',Title1];
            end
            title(Title)
            PosIndex=find(RepNb{1}{TypeL}{SimL}==NetNb);

            Corr=[];
            Anti=[];
            Pv=[];
            for PosL=1:length(PosIndex)
                %for the current pair find the networks wherecorr, anti and pv values are good
                GoodPos=find(Corrs{TypeL}{SimL}{PosIndex(PosL)}>=Limit.corr&...
                    Antis{TypeL}{SimL}{PosIndex(PosL)}<=Limit.anti&...
                    Pvs{TypeL}{SimL}{PosIndex(PosL)}<=Limit.pv);
                %sum the total number of networks where values are good
                GoodNb(TypeL,SimL)=GoodNb(TypeL,SimL)+length(GoodPos);

                %distribution of number of networks where values are good
                SelNb{TypeL,SimL}(length(GoodPos)+1)=SelNb{TypeL,SimL}(length(GoodPos)+1)+1;
                %recover allgood values
                if ~isempty(GoodPos)
                    Corr=[Corr,Corrs{TypeL}{SimL}{PosIndex(PosL)}(GoodPos)];
                    Anti=[Anti,Antis{TypeL}{SimL}{PosIndex(PosL)}(GoodPos)];
                    Pv=[Pv,Pvs{TypeL}{SimL}{PosIndex(PosL)}(GoodPos)];
                end
            end
            plot(Corr,Anti,'c+','markersize',3)
            plot(Corr,Pv,'m+','markersize',3)

            ylabel('anti & pv(overlap)')

            %the same but with the complement (bad values)
            Corr=[];
            Anti=[];
            Pv=[];
            for PosL=1:length(PosIndex)
                BadPos=find(Corrs{TypeL}{SimL}{PosIndex(PosL)}<Limit.corr|...
                    Antis{TypeL}{SimL}{PosIndex(PosL)}>Limit.anti|...
                    Pvs{TypeL}{SimL}{PosIndex(PosL)}>Limit.pv);
                BadNb(TypeL,SimL)=BadNb(TypeL,SimL)+length(BadPos);
                if ~isempty(BadPos)
                    Corr=[Corr,Corrs{TypeL}{SimL}{PosIndex(PosL)}(BadPos)];
                    Anti=[Anti,Antis{TypeL}{SimL}{PosIndex(PosL)}(BadPos)];
                    Pv=[Pv,Pvs{TypeL}{SimL}{PosIndex(PosL)}(BadPos)];
                end
            end
            plot(Corr,Anti,'bx','markersize',3)
            plot(Corr,Pv,'rx','markersize',3)
            set(gca,'box','on')
            xlabel(sprintf('corr (%u%% of good)',round(GoodNb(TypeL,SimL)*100/(GoodNb(TypeL,SimL)+BadNb(TypeL,SimL)))))
        end


        Position=get(gcf,'position');
        Position(3:4)=[600,450];

        set(gcf,'position',Position)
        try
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        catch
            mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        end
    end
    saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank),'png')
    close(h)    
end



%% FIG39 - relationships between CORR, ANTI and Pearson's correlation coefficient in first
% network
if ~isempty(find(FigRanks==39))
    FigRank=39;
    h=figure;
    set(h,'name',sprintf('FIG%u - m%u: CORR vs PEARSON',FigRank,ChipRank))
    set(gcf,'color',[1,1,1])
    for TypeL=1:2
        for SimL=1:3
            GoodPos=find(SIM{1}{TypeL,SimL}{1}.corr>=Limit.corr&...
                SIM{1}{TypeL,SimL}{1}.anti<=Limit.anti&...
                SIM{1}{TypeL,SimL}{1}.pv{PvCorrRank}<=Limit.pv);
            subplot(4,3,(TypeL-1)*6+SimL)
            plot(SIM{1}{TypeL,SimL}{1}.corr,Pearson{TypeL}{SimL}.rankCorr,'m+');
            hold on
            plot(SIM{1}{TypeL,SimL}{1}.corr(GoodPos),Pearson{TypeL}{SimL}.rankCorr(GoodPos),'r+');
            xlabel('corr')
            if TypeL==1
                if SimL==1
                    title('In')
                    ylabel('SINGLE')
                elseif SimL==2
                    title('Out')
                else
                    title('Random pairs')
                end
            end
            if TypeL==2
                if SimL==1
                    ylabel('MULTIPLE')
                end
            end
            subplot(4,3,(TypeL-1)*6+3+SimL)
            plot(SIM{1}{TypeL,SimL}{1}.anti,Pearson{TypeL}{SimL}.rankCorr,'c.');
            hold on
            plot(SIM{1}{TypeL,SimL}{1}.anti(GoodPos),Pearson{TypeL}{SimL}.rankCorr(GoodPos),'b+');
            xlabel('anti')
            if TypeL==1
                if SimL==1
                    ylabel('SINGLE')
                end
            end
            if TypeL==2
                if SimL==1
                    ylabel('MULTIPLE')
                end
            end
        end
    end
    set(h,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[800,1200];
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


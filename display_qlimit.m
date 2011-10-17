%============================%
% FUNCTION CALCULATE_LIMITS  %
%============================%

% CALCULATE_LIMITS uses several small networks (1024 comparisons : 32x32 biol cond)
% to find the distributions of p-value of overlaping between neighbourhood and of
% positive and negative correlation of probe sets that target the same groups of transcripts

%INPUT

% 1        ModelRank: chip model rank
% 2          FigList: indicates figures that are to be displayed (set it to [] in order
% 3 FirstNetRankList: first list of networks
% 4           PvCorrRank: pv(overlap) is calculated for corr limit >[0,40,50,60]. PvCorrRank
%                     indicates the corr limit to be used, by giving its index in
%                     the corr list(,[0,40,50,60])
% varargin:
% 5     NoQlimitFlag: indicates if the second group of network is a set of no qlimit network
% 6   SndNetRankList: second list of networks

% FIGURES
% LimitLim is a list of limits calculated for each pair of probe sets which have a significative
% corr within x networks (5th percentile for corr and 95th percentile for anti and
% pv(overlap)).
% LimitVal is a list of all corr, anti and pv(overlap) values for all probe sets which have a
% significative corr within x networks.
% FIG1 - Distibution of corr, anti and pv(overlap) in LimitLim for pairs of probe sets that
%        target a single gene
% FIG2 - Distibution of corr, anti and pv(overlap) in LimitLim for pairs of probe sets that
%        target multiple genes
% FIG3 - Distibution of corr, anti and pv(overlap) in LimitVal for pairs of probe sets that
%        target a single gene
% FIG4 - Distibution of corr, anti and pv(overlap) in LimitVal for pairs of probe sets that
%        target multiple genes



% calculate_limits(8,1,[7:21])
% calculate_limits(8,1,[7:21],0,[24:38])
% without qlimit
% calculate_limits(8,1,[24:38],1,[39:53])

%EXTERNAL FILES

%OUTPUT PARAMETERS
% 1 Limits

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

function [Limit]=calculate_limits(Species,ChipRank,FigRanks,FirstNetRanks,PvCorrRank,NoQlimitFlag,SndNetRanks)
global K

% PS ASSIGNATION
NetNb=[];
NetRankList={};
DataDir=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(DataDir)

NetRankList{1}=FirstNetRanks;
NetNb(1)=length(NetRankList{1});
NetRankListNb=1;

ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName);
PsNb=K.chip.probesetNb(ChipPos);
ProbeNb=K.chip.probeNb(ChipPos);




if NoQlimitFlag
    if length(SndNetRanks)~=length(FirstNetRanks)
        errodlg('NetRankLists do not have the same length')
        error('process canceled')
    end
end
NetRankList{2}=SndNetRanks;
NetNb(2)=length(NetRankList{2});


%load Sim
SIM=cell(NetRankListNb);
for ListL=1:NetRankListNb
    SIM{ListL}=cell(1,2);
    for NetL=1:NetNb(ListL)
        eval(sprintf('load m%u_n%u_nodesim_probenb1_single;',ChipRank,NetRankList{ListL}(NetL)));
        SIM{ListL}{1,1}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb1_single_testout;',ChipRank,NetRankList{ListL}(NetL)));
        SIM{ListL}{1,2}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb1_single_testhigh;',ChipRank,NetRankList{ListL}(NetL)));
        SIM{ListL}{1,3}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb1_multiple;',ChipRank,NetRankList{ListL}(NetL)));
        SIM{ListL}{2,1}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb1_multiple_testout;',ChipRank,NetRankList{ListL}(NetL)));
        SIM{ListL}{2,2}{NetL}=Sim;

        eval(sprintf('load m%u_n%u_nodesim_probenb1_multiple_testhigh;',ChipRank,NetRankList{ListL}(NetL)));
        SIM{ListL}{2,3}{NetL}=Sim;
    end
end

eval(sprintf('load m%u_pearson_probenb1_single;',ChipRank));
p{1}=Pearson;
eval(sprintf('load m%u_pearson_probenb1_multiple;',ChipRank));
p{2}=Pearson;
Pearson={};
Pearson=p;
clear p;
Letters='abcdefghijklmnop';

%load CompInfo
load(sprintf('m%u_compinfo_probenb1_single.mat',ChipRank));
CoupleInfo{1}=CompInfo;
load(sprintf('m%u_compinfo_probenb1_multiple.mat',ChipRank));
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
% recover either the list of intermediate limits for each pair of probe sets,
% of a single list of all the values
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

for TypeL=1:2
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
                        %                         if ValL==3
                        %                             InfPos=isinf(CurrVal{PosIndex(PosL)});
                        %                             NegPos=CurrVal{PosIndex(PosL)}<0;
                        %                             CurrVal{PosIndex(PosL)}(InfPos)=350;
                        %                             CurrVal{PosIndex(PosL)}(InfPos&NegPos)=-350;
                        %                         end
                        [CurrVal{PosIndex(PosL)},SortIndex]=sort(CurrVal{PosIndex(PosL)});
                        if ValL==1
                            %another possibility for calculating CurrLim
                            %CurrLim(PosL)=CurrVal{PosIndex(PosL)}(min(RepL,2));
                            CurrLim(PosL)=mean(CurrVal{PosIndex(PosL)})-std(single(CurrVal{PosIndex(PosL)}));
                            CurrValues((PosL-1)*RepL+1:PosL*RepL)=CurrVal{PosIndex(PosL)};
                        else
                            %another possibility for calculating CurrLim
                            %CurrLim(PosL)=CurrVal{PosIndex(PosL)}(max(RepL-1,1));
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
            end
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



%Use only InSim values
Limit.corr=LimitVal.corr{1,1}(NetNb);
Limit.anti=LimitVal.anti{1,1}(NetNb);
Limit.pv=LimitVal.pv{1,1}(NetNb);


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





    %CALCUL PAR ETAPES SUCCESSIVES DES PARAMETRES STATISTIQUES
    PairNb=zeros(2,3);
    for TypeL=1:2
        for SimL=1:3
            PairNb(TypeL,SimL)=length(MeanPv{1}{TypeL}{SimL});
        end
    end




    %% FIG17
    % caractère très discriminant de pv et de corr
    % on distingue les couples dont la corrélation est nulle dans tous les réseaux (marqués "NoCorr")
    % et les couples dont la corrélation n'est pas nulle dans au moins un réseau (marqués "IsCorr").
    % Pour chacun des couples de ce dernier type on distingue les réseaux dans lesquelles la corrélation est positive
    % et marqués "IsCorr Corr>0" et les réseaux dans lesquels la corrélation est nulle et marqués IsCorr Corr=0.
    % distribution de mean(pv(overlap))
    % Sont aussi considérés les couples de probe set ciblant un même gène et localisés soit dans
    % les exons (sIn), soit en dehors (sOut). Il existe également des
    % couples de probe sets formés aléatoirement, ciblant chacun un gène différent (sRand).
    % Ces trois catégories existent également pour les probe sets
    % ciblant plusieurs gènes et sont notées respectivement mIn, mOut et mRand.
    % Enfin, chaque réseau peut être représenté éventuellement sous
    % deux formes. La forme "qlimit" est celle utilisée habituellement
    % qui a subi l'étape d'elimination des valeurs corr et anti non
    % significatives. La forme no "qlimit" est le réseau ayant conservé
    % toutes les valeurs calculées de corr et de anti ce qui permet de
    % savoir quelles sont les caractéristiques repectives des valeurs
    % qui ont été éliminées ou conservées.
    % Pour ces différentes catégories, plusieurs distributions sont
    % représentées :
    % CADRANT A : mean(pv(overlap)) : moyenne des p-values de la similarité de voisinage de deux probe
    % sets dans plusieurs réseaux (soit tous pour les couples NoCorr soit IsCorr Corr=0 ou IsCorrCorr>0 pour les couples
    % IsCorr)
    % => peu de différences entre catégories s (single, un seul gène
    % cible) et m (multiple, plusieurs gènes ciblés)
    % => grande différence entre IsCorr Corr>0 (pv pratiquement toutes
    % <=Limit.pv) et NoCorr (sigmoide centréé sur pv=0)
    % => situation intermédaire pour IsCorr Corr=0 (pv <0 mais seulement
    % 25<= Limit.pv
    % La distribution des pv(overlap) est identique pour les couples pris au hasard
    % et qui ont une corrélation supérieure à 0.
    % CADRANT B : mean(corr) pour IsCorr Corr>0
    % => ditribution pour In et Out décalée fortement par rapport a
    % Rand (médiane à >60 contre ~40 ; résultats attendus,la corrélation entre probe sets
    % ciblant le même gène doit être plus forte que celle existant
    % entre deux probe sets pris au hasard).
    % Mais
    % CADRANT C : mean(corr) IsCorr Corr=0 et NoCorr dans no qlimit
    % => pour IsCorr même position relative que dans cadrant B, mais
    % courbes décalées vers la gauche (mediane In,Out <=50 et Rand ~35)
    % =>pour NoCorr les médianes sont encore plus basses (~25}
    % CADRANT D : mean(corr) IsCorr Corr>0 dans qlimit vs mean(corr)
    % IsCorr Corr=0 dans no qlimit
    % => autre représentation donnant plus d'informations que dans
    % cadrant C: les moyennes des valeurs éliminées par qlimit sont dans leur
    % grande majorité inférieures aux moyennes des valeurs
    % conservées. Mais la corrélation est toujours forte même et surtout pour les couples Rand ce qui est
    % inattendu. En effet s'il n'est pas surprenant que des couples
    % de probe set ciblant le même gènes aient la même valeur de
    % corrélation quel que soit le réseau (même si certaines sont
    % considérées comme non significatives par l'application de
    % qlimit), il est est en revanche étonnant que quel que soit le
    % réseau considéré, les valeurs brut de corrélation soient
    % très similaire pour des couples de probe sets aléatoires.
    % CONCLUSION
    % Le processus de sélection par qlimit semble très efficace dans le
    % sens où pour les probe sets ciblant le même gène, lorsque les
    % conditions choisies pour construire le réseau font que la
    % corrélation est plus faible (bruit de fond plus important masquant partiellement la corrélation ?)
    % celle ci n'est pas considérée comme significative. Et ceci se
    % répercute sur le voisinage avec des valeurs de pv(overlap) > Limit.pv.
    % La sélection par qlimit est donc sans doute nécessaire pour
    % faire apparaître la structure du réseau (Si l'on se basait sur la
    % forte corrélation des valeurs de corr pour les couples pris au
    % hasard on n'éliminerait aucune valeur et le réseau serait saturé
    % de liens).
    if ~isempty(find(FigRanks==17))
        FigRank=17;
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
                    l=plot(0,0,sprintf('%c%c',Colors(SimL)),Symbols(TypeL));
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
    end


    %% FIG18
    % Vérification des hypothèses de la figure 1 en groupant les couples en fonction du nombre de replicates significatifs
    % (nombre de réseaux pour lesquels un couple a des valeurs corr significatives).
    % On compare MeanCorr calculé dans les réseaux qlimit et OutMeancorr calculé dans les réseaux no qlimit
    % (on marque de manière spécifique les points pour lesquels pv>-5 dans les réseaux qlimit).
    % Même conclusions que pour Cadrant D de la ficure 1 => très bonne corrélation des valeurs Corr entre les
    % réseaux qlimit et no qlimit, mais valeurs systématiquement plus basses dans réseaux no qlimit.
    % On compare également StdCorr et OutStdCorr : pour StdCorr la moyenne de std diminue quand PosNb diminue
    % (de 4.5 pour PosNb=14 à 2.3 pour PosNb=1), et l'écart type augmente très légèrement (1.5 à 2.1).
    % En revanche dans no qlimit, les variations sont bien plus fortes (mean : 5.5 à 17; std : 3.2 à 1.4).
    % ce qui montre la pertinence du filtrage des valeurs significatives.
    if ~isempty(find(FigRanks==18))& NetRankListNb==2
        Colors='rkbcygm';
        CNb=length(Colors);
        SubNb1=floor(sqrt(NetNb(1)));
        SubNb2=ceil(NetNb(1)/SubNb1);
        Colors='rmg';
        Lines={'-',':'};

        SubNb1=floor(sqrt(NetNb(1)));
        SubNb2=ceil(NetNb(1)/SubNb1);
        Colors='rmg';
        Lines={'-',':'};
        LetterRank=0;
        for PlotL=1:2
            for TypeL=1:2
                for SimL=1:3
                    h=figure;
                    LetterRank=LetterRank+1;
                    if TypeL==1
                        if PlotL==1
                            set(h,'name',sprintf('FIG%u%c - m%u: SINGLE (MeanVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        else
                            set(h,'name',sprintf('FIG%u%c - m%u: SINGLE (StdVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        end
                    else
                        if PlotL==1
                            set(h,'name',sprintf('FIG%u%c - m%u: MULTIPLE (MeanVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        else
                            set(h,'name',sprintf('FIG%u%c - m%u: MULTIPLE (StdVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        end
                    end
                    if PlotL==1
                        Val1=MeanCorr{1}{TypeL}{SimL};
                        Val2=OutMeanCorr{2}{TypeL}{SimL};
                    else
                        Val1=StdCorr{1}{TypeL}{SimL};
                        Val2=OutStdCorr{2}{TypeL}{SimL};
                    end
                    PvPos=MeanPv{1}{TypeL}{SimL}>-5;
                    for RepL=1:NetNb(1)-1
                        RepPos=RepNb{1}{TypeL}{SimL}==RepL;
                        CurrPvPos=find(PvPos&RepPos);
                        if ~isempty(find(RepPos))
                            subplot(SubNb1,SubNb2,RepL)
                            hold on
                            %plot(Val1(PosIndex),Val2(PosIndex),sprintf('%c+',Colors(SimL)),'markersize',3)
                            plot(Val1(RepPos),Val2(RepPos),sprintf('%c+',Colors(SimL)),'markersize',3)
                            plot(Val1(CurrPvPos),Val2(CurrPvPos),'bo')
                            line([0,100],[0,100])
                            Pears=corrcoef(Val1(RepPos),Val2(RepPos));
                            CurrMean1=mean(Val1(RepPos));
                            CurrStd1=std(Val1(RepPos));
                            CurrMean2=mean(Val2(RepPos));
                            CurrStd2=std(Val2(RepPos));
                            if length(Pears)==2
                                Pears=Pears(2);
                            end

                            set(gca,'box','on')
                            if PlotL==1
                                switch SimL
                                    case 1
                                        title(sprintf('MeanCorr IN EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 2
                                        title(sprintf('MeanCorr OUT EXONS RepNb=%u (corr=%.2f\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 3
                                        title(sprintf('MeanCorr RAND RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                end
                                set(gca,'xlim',[0,100])
                                set(gca,'ylim',[0,100])
                                xlabel('mean corr in qlimit')
                                ylabel('mean corr in no qlimit')

                            else
                                switch SimL
                                    case 1
                                        title(sprintf('StdCorr IN EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 2
                                        title(sprintf('StdCorr OUT EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 3
                                        title(sprintf('StdCorr RAND RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                end
                                set(gca,'xlim',[0,20])
                                set(gca,'ylim',[0,20])
                                xlabel('std corr in qlimit')
                                ylabel('std corr in no qlimit')
                            end
                        end
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
                    saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(letterRank)),'png')
                end
            end
        end
    end




    %% FIG19
    % Vérification des hypothèses de la figure 1 en groupant les couples en fonction du nombre de replicates significatifs
    % (nombre de réseaux pour lesquels un couple a des valeurs corr significatives).
    % On compare MeanPv et OutMeanPv (TOUS DEUX CALCULÉS SUR LES RÉSEAUX NO QLIMIT CAR L'ON
    % NE PEUT PAS COMPARER DIRECTEMENT PV(OVERLAP) ENTRE LES DEUX TYPES DE RÉSEAUX (dans les
    % réseaux no qlimit le nombre de liens est bien plus grand et de ce fait les pv sont supérieures)
    % (on marque de manière spécifique les points pour lesquels pv>-5 dans les réseaux qlimit).
    % => DANS LES RÉSEAUX NO QLIMIT TOUS LES VALEURS DE PV SEMBLE ÊTRE NÉGATIVES (VOIR AUSSI RAND).
    % => La corrélation est nettement moins bonne que pour Corr entre les réseaux qlimit et no qlimit.
    %
    % => Il existe un effet d'échatillonnage qui complique l'interprétation:
    % On observe une rotation de l'ensemble des points lorque l'on passe de RepNb=1 à RepNb=NetNb(1)
    % Cela s'explique par une dispersion des valeurs qui dépend de repNb (alors que la moyenne des mean(Pv) reste remarquablement
    % constante entre -6.8 et -7.1 pour les réseaux qlimit et -6.3 à -7.1 pour les réseaux no qlimit)
    % Pour RepNb passant de 14 à 1 , std(mean(Pv)) passe de 1.1 à 3.6 dans les réseaux qlimit
    % et de 1.3 à 3.5 dans les réseaux no qlimit.
    % En fait la valeur de std(mean(Pv)) augmente principalement en dessous de RepNb=4 (dépasse alors 2);
    % ET L'ON PEUT SE DEMANDER SI LES CORRÉLATIONS OBSERVÉES SONT SIGNIFICATIVES.
    % Mêmes observations avec std(Pv) concernant l'effet d'échantillonnage et la rotation de
    % l'ensemble des points.
    % pour RepNb passant de 14 à 2 :
    % moyenne des mean(StdPv)
    % qlimit : 3.2 à 1.7; no qlimit: 3.4 à 1.8
    % moyenne des std(StdPv)
    % qlimit : 0.5 à 1.5; no qlimit: 0.6 à 1.6

    if ~isempty(find(FigRanks==19))&NetRankListNb==2

        Colors='rkbcygm';
        CNb=length(Colors);
        SubNb1=floor(sqrt(NetNb(1)));
        SubNb2=ceil(NetNb(1)/SubNb1);
        Colors='rmg';
        Lines={'-',':'};

        SubNb1=floor(sqrt(NetNb(1)));
        SubNb2=ceil(NetNb(1)/SubNb1);
        Colors='rmg';
        Lines={'-',':'};
        LetterRank=0;
        for PlotL=1:2
            if NetRankListNb==2
                for TypeL=1:2
                    for SimL=1:3
                        LetterRank=LetterRank+1;
                        h=figure;
                        if TypeL==1
                            if PlotL==1
                                set(h,'name',sprintf('FIG%u%c - m%u: SINGLE (MeanVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                            else
                                set(h,'name',sprintf('FIG%u%c - m%u: SINGLE (StdVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                            end
                        else
                            if PlotL==1
                                set(h,'name',sprintf('FIG%u%c - m%u: MULTIPLE (MeanVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                            else
                                set(h,'name',sprintf('FIG%u%c - m%u: MULTIPLE (StdVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                            end
                        end
                        if PlotL==1
                            Val1=MeanPv{2}{TypeL}{SimL};
                            Val2=OutMeanPv{2}{TypeL}{SimL};
                        else
                            Val1=StdPv{2}{TypeL}{SimL};
                            Val2=OutStdPv{2}{TypeL}{SimL};
                        end
                        PvPos=MeanPv{1}{TypeL}{SimL}>Limit.pv;
                        for RepL=1:NetNb(1)-1
                            RepPos=RepNb{1}{TypeL}{SimL}==RepL;
                            CurrPvPos=find(PvPos&RepPos);
                            length(RepPos)
                            if ~isempty(RepPos)
                                subplot(SubNb1,SubNb2,RepL)
                                hold on
                                plot(Val1(RepPos),Val2(RepPos),sprintf('%c+',Colors(SimL)),'markersize',3)
                                %plot(Val1(RepPos),Val2(RepPos),sprintf('%c+',Colors(SimL)))
                                plot(Val1(CurrPvPos),Val2(CurrPvPos),'b+','markersize',3)
                                line([-15,15],[-15,15])
                                Pears=corrcoef(Val1(RepPos),Val2(RepPos));
                                CurrMean1=mean(Val1(RepPos));
                                CurrStd1=std(Val1(RepPos));
                                CurrMean2=mean(Val2(RepPos));
                                CurrStd2=std(Val2(RepPos));
                                if length(Pears)==2
                                    Pears=Pears(2);
                                end

                                set(gca,'box','on')
                                if PlotL==1
                                    switch SimL
                                        case 1
                                            title(sprintf('MeanPv IN EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                        case 2
                                            title(sprintf('MeanPv OUT EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                        case 3
                                            title(sprintf('MeanPv RAND RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    end
                                    set(gca,'xlim',[-15,15])
                                    set(gca,'ylim',[-15,15])
                                    xlabel('mean pv IsCorr Corr>0')
                                    ylabel('mean pv IsCorr Corr=0')
                                else
                                    switch SimL
                                        case 1
                                            title(sprintf('StdPv IN EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                        case 2
                                            title(sprintf('StdPv OUT EXONS RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                        case 3
                                            title(sprintf('StdPv RAND RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    end
                                    set(gca,'xlim',[0,10])
                                    set(gca,'ylim',[0,10])
                                    xlabel('std pv IsCorr Corr>0')
                                    ylabel('std pv IsCorr Corr=0')

                                end
                            end
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
                        saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(letterRank)),'png')
                    end
                end
            end
        end
    end
    %% FIG20
    % Pour répondre on prend des réseaux choisi aléatoirement mais en respectant le nombre de réseaux dans lesquels
    % corr est significatifs (PosNb).
    % Cela montre que quand le couple est très reproductible (par ex PosNb=14) si l'on prend un réseau de manière alétoire, on a
    % de grande chance de tomber sur un des PosNb réseaux et la valeur est comprise entre Limit.pv et -5 et non pas entre Limit.pv et 0.
    % De fait mean(mean(Pv))=-7.0 aussi bien pour Kept que pour NotKept et std(mean(Pv)) est stable (passe de 1.4 à 1.6).
    % En revanche quand PosNb est faible cela veut dire que les couples concernés sont peu 'fiables' et effectivement
    % la valeur de pv oscille dans ce cas entre Limit.pv et 0 même pour le ou les réseaux non choisi aléatoirement (figure 3).
    % De fait mean(mean(Pv))=-6.6 aussi bien pour Kept que pour NotKept mais std(mean(Pv)) varie beaucoup et passe de 2.8 pour Kept à 1.4 pour NoKept.
    % Il y a donc superposition d'un effet d'échantillonage (la moyenne de plusieurs valeurs est toujours comprise entre Limit.pv et -5
    % ET UN EFFET DE CORRÉLATION RÉEL et de valeurs différentes pour les couples fiables (avec PosNb pas trop petit). Dans ce dernier
    % cas les valeurs individuelles sont comprises entre Limit.pv et -5 (A VERIFIER).
    % Même conclusions avec std(Pv).
    if ~isempty(find(FigRanks==20)) & NetRankListNb==2
        FigRank=20;
        Colors='rkbcygm';
        CNb=length(Colors);
        SubNb1=floor(sqrt(NetNb(1)));
        SubNb2=ceil(NetNb(1)/SubNb1);
        Colors='rmg';
        Lines={'-',':'};

        SubNb1=floor(sqrt(NetNb(1)));
        SubNb2=ceil(NetNb(1)/SubNb1);
        Colors='rmg';
        Lines={'-',':'};
        for PlotL=1:2
            for TypeL=1:2
                for SimL=1:3
                    h=figure;
                    if TypeL==1
                        if PlotL==1
                            set(h,'name',sprintf('FIG%u%c - m%u: SINGLE (MeanVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        else
                            set(h,'name',sprintf('FIG%u%c - m%u: SINGLE (StdVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        end
                    else
                        if PlotL==1
                            set(h,'name',sprintf('FIG%u%c - m%u: MULTIPLE (MeanVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        else
                            set(h,'name',sprintf('FIG%u%c - m%u: MULTIPLE (StdVal,InSim%u)',FigRank,Letters(LetterRank),ChipRank),SimL)
                        end
                    end
                    if PlotL==1
                        Val1=RandMeanPv{2}{TypeL}{SimL};
                        Val2=RandOutMeanPv{2}{TypeL}{SimL};
                    else
                        Val1=RandStdPv{2}{TypeL}{SimL};
                        Val2=RandOutStdPv{2}{TypeL}{SimL};
                    end
                    for RepL=1:NetNb(1)-1
                        RepPos=RepNb{1}{TypeL}{SimL}==RepL;
                        length(RepPos)
                        if ~isempty(RepPos)
                            subplot(SubNb1,SubNb2,RepL)
                            hold on
                            plot(Val1(RepPos),Val2(RepPos),sprintf('%c+',Colors(SimL)),'markersize',3)

                            %plot(Val1(RepPos),Val2(RepPos),sprintf('%c+',Colors(SimL)))
                            Pears=corrcoef(Val1(RepPos),Val2(RepPos));
                            CurrMean1=mean(Val1(RepPos));
                            CurrStd1=std(Val1(RepPos));
                            CurrMean2=mean(Val2(RepPos));
                            CurrStd2=std(Val2(RepPos));
                            if length(Pears)==2
                                Pears=Pears(2);
                            end
                            set(gca,'box','on')
                            if PlotL==1
                                switch SimL
                                    case 1
                                        title(sprintf('Mean Pv IN RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 2
                                        title(sprintf('Mean Pv OUT RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 3
                                        title(sprintf('Mean Pv RAND RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                end
                                line([-15,5],[-15,5])
                                line([-15,-5],[-5,-5])
                                line([-5,-5],[-15,-5])
                                xlabel('mean pv IsCorr Rand Corr>0')
                                ylabel('mean pv IsCorr Rand Corr=0')
                            else
                                switch SimL
                                    case 1
                                        title(sprintf('Std Pv IN RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 2
                                        title(sprintf('Std Pv OUT RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                    case 3
                                        title(sprintf('Std Pv RAND RepNb=%u (corr=%.2f)\nmeanX=%.1f, stdX==%.1f,meanY=%.1f, stdY==%.1f',RepL,Pears,CurrMean1,CurrStd1,CurrMean2,CurrStd2))
                                end
                                line([0,10],[0,10])
                                xlabel('std pv IsCorr Rand Corr>0')
                                ylabel('std pv IsCorr Rand Corr=0')
                            end
                        end
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
                    saveas(h,sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(letterRank)),'png')
                end
            end
        end
    end

    %% FIG21
    if ~isempty(find(FigRanks==21))
        FigRank=21;
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
    end


    %% FIG22
    if ~isempty(find(FigRanks==22))
        FigRank=22;
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
                CurrNetRank=NetRankList{ListL}(NetL);
                CNetFile=sprintf('c_m%u_n%u.4mat',ChipRank,CurrNetRank);
                ANetFile=sprintf('a_m%u_n%u.4mat',ChipRank,CurrNetRank);
                NetDir=fullfile(K.dir.net,sprintf('m%03u',ChipRank),sprintf('n%05u',CurrNetRank));
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


        %FIG22a - FIG22b
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
            else
                saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
            end
        end


        %FIG22c - FIG22d
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
            else
                saveas(h,sprintf('m%u_fig%ud.png',ChipRank,FigRank),'png')
            end
        end


        %FIG22e - FIG22f
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
            else
                saveas(h,sprintf('m%u_fig%uf.png',ChipRank,FigRank),'png')
            end
        end



        %FIG22g - FIG22h
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
            else
                saveas(h,sprintf('m%u_fig%uh.png',ChipRank,FigRank),'png')
            end
        end




        %FIG22i
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



        %FIG22j
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
    end

    %% FIG23
    %Caracteristics of probes

    if ~isempty(find(FigRanks==23))
        FigRank=23;

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
                title('mean common gene')
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
                title('common transcrips')
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
                set(gca,'xlim',[0,5])
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
                cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
            end
            if FigL==1
                saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
            else
                saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
            end
        end
    end





    %% FIG24
    if ~isempty(find(FigRanks==24))
        FigRank=24;
        %comparaisons de deux séries de réseaux qlimit
        Colors='rb';
        Lines={'-','-.'};
        h=figure;
        if NetRankListNb==2 & NoQlimitFlag==0
            h=figure;
            set(h,'name',sprintf('FIG%u - m%u: comparison of two series of networks',FigRank,ChipRank))
            for TypeL=1:2
                for SimL=1:2
                    subplot(2,3,(TypeL-1)*3+1)
                    hold on
                    for ListL=1:2
                        PosIndex=find(RepNb{ListL}{TypeL}{SimL}==NetNb(ListL));
                        Val=MeanCorr{ListL}{TypeL}{SimL}(PosIndex);
                        plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(SimL),Lines{ListL}))
                    end
                    set(gca,'box','on')
                    title('mean corr')
                    subplot(2,3,(TypeL-1)*3+2)
                    hold on
                    for ListL=1:2
                        PosIndex=find(RepNb{ListL}{TypeL}{SimL}==NetNb(ListL));
                        Val=MeanAnti{ListL}{TypeL}{SimL}(PosIndex);
                        plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(SimL),Lines{ListL}))
                    end
                    set(gca,'box','on')
                    title('mean anti')
                    subplot(2,3,(TypeL-1)*3+3)
                    hold on
                    for ListL=1:2
                        PosIndex=find(RepNb{ListL}{TypeL}{SimL}==NetNb(ListL));
                        Val=MeanPv{ListL}{TypeL}{SimL}(PosIndex);
                        plot(sort(Val),1/length(Val):1/(length(Val)+1):1,sprintf('%c%s',Colors(SimL),Lines{ListL}))
                    end
                    set(gca,'box','on')
                    title('mean pv')
                end
            end
            set(gcf,'color',[1,1,1])
        end
    end

    %% FIG25
    if ~isempty(find(FigRanks==25))
        FigRank=25;
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
        end
    end



    %% FIG26
    if ~isempty(find(FigRanks==26))
        FigRank=26;
        %Relation between Corr,Anti and Pv values


        h(1)=figure;
        set(h(1),'name',sprintf('FIG%ua - m%u: CORR,ANTI,PV',FigRank,ChipRank))
        set(gcf,'color',[1,1,1])
        h(2)=figure;
        set(h(2),'name',sprintf('FIG%ub - m%u: In Out Stats',FigRank,ChipRank))
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
                figure(h1)
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

                figure(h2)
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
            for FigL=1:2
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
            end
        end
    end
end

%% FIG27 - relationships between CORR, ANTI and Pearson's correlation coefficient in first
% network
if ~isempty(find(FigRanks==27))
    FigRank=27;
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
end

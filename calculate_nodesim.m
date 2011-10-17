%============================%
% FUNCTION CALCULATE_NODESIM %
%============================%

%CALCULATE_NODESIM calculates node neighbourhood similarity between probesets that target
% the same genes (duplicates)

%INPUT PARAMETERS
% 1     TestFlag: if =1 indicates that similarity is calculated to find limits on corr, anti and
%                 pv(node neighbourhood similarity); if =0 indicates that similarity is calculated to merge probe
%                 sets
% 2 NotFoundFlag: indicates that duplicates not found in a first round are processed
%                 (TestFlag is equal to 0 in this case)
% 3 ProbeNbLimit: minimal number of probes targeting a gene (used to make
%                 a bipartition of genes)
% 4     ChipRank: rankf of chip model
% 5  NetRankList: a list of nets, designed by their rank, used to calculate node smilarity
% 6 PvCorrLimits: indicates the corr values that must be used to select the
%                 probe sets used to calculate Pv of node neighbourhood similarity
% 7     AceFlag: indicates if AcView data are used
% 8    IdemFlag: indicates if probe set order is identical in file used by PsawnPy et in
%                 networks (CVM)

%OUTPUT FILES

%INTERNAL FUNCTION

%''''''''''''''''''%
% FUNCTION NODESIM %
%''''''''''''''''''%

%INPUT PARAMETERS
% 1    ChipRank: rank of chip model
% 2      NetRank: rank of the current net
% 3 ProbeNbLimit: minimal number of probes targeting a gene (used to make a 
%                 bipartition of genes)
% 4          Dup: samples of pairs of probesets targeting the same gene(s)
% 5 MultipleFlag: used if ProbeNbLImit==1 => allows to consider probeset with
%                 only one target (MultipleFlag=0) or with several targets (MultipleFlag=1)
% 6      DupType: indicates different type of duplicate

%OUTPUT FILES
% For each type of duplicate a file is written. Example: 
% sprintf('m%u_n%u_nodesim_probenb%u_multiple.mat',ChipRank,NetRank,ProbeNbLimit)
% This file contains the structured variable Sim with the following fields.
%
%            Sim.corr: positive correlations values for each pair of probesets
%            Sim.anti: negative correlations values for each pair of probesets
%     Sim.firstPsRank: probeset ranks of the first probeset in each pair
%       Sim.sndPsRank: probeset ranks of the second probeset in each pair
% Sim.commonNodeNb{i}: common number of neighbouhrs
%  Sim.firstNodeNb{i}: number of neighbouhrs for the first probeset
%    Sim.sndNodeNb{i}: number of neighbouhrs for the second probeset
%      Sim.overlap{i}: percentage overlap (common number * 100 / min(firstNodeNb,sndNodeNb))
%           Sim.pv{i}: p-value calculated with hypergeometric distribution;



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

function calculate_nodesim(TestFlag,NotFoundFlag,ProbeNbLimit,ChipRank,NetRankList,PvCorrLimits,AceFlag,IdemFlag)
global K

ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
if isempty(ChipPos)
    h=errordlg(sprintf('chip m%u does not exist',ChipRank));
    waitfor(h)
    error('process canceled')
end
Species=K.chip.species{ChipPos};
ProbeNb=K.chip.probeNb(ChipPos);
PsNb=K.chip.probesetNb(ChipPos);

% load information on probesets by gene (created manually with import_info)
% for a given gene are listed all the probesets that have a given
% number of probes inside its exons
%also are given targeted and not targeted list of transcripts

%load probesets with Ensembl gene information
cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)))
if exist(fullfile(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)),sprintf('m%u_probesets_by_ensembl_gene.mat',ChipRank)),'file')
    eval(sprintf('load m%u_probesets_by_ensembl_gene',ChipRank))
    TargetedGenes{1}=ETargetedGenes;
    TargetingPsRanks{1}=ETargetingPsRanks;
    ExonNames{1}=EExonNames;
    ExonProbeNbs{1}=EExonProbeNbs;
    GroupRanks{1}=EGroupRanks;
    GroupProbeNbs{1}=EGroupProbeNbs;
    TargetedTs{1}=ETargetedTs;
    NotTargetedTs{1}=ENotTargetedTs;
    PsInfo{1}=EPsInfo;
    clear EGeneName EGeneNames ETargetedGenes ETargetingPsRanks EPsNames EPsRanks EExonNames EExonProbeNbs EGroupRanks EGroupProbeNbs ETargetedTs ENotTargetedTs EPsInfo
end

%load probesets with AceView gene information
if exist(fullfile(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)),sprintf('m%u_probesets_by_aceview_gene.mat',ChipRank)),'file')
    eval(sprintf('load m%u_probesets_by_aceview_gene',ChipRank))
    TargetedGenes{2}=ATargetedGenes;
    TargetingPsRanks{2}=ATargetingPsRanks;
    ExonNames{2}=AExonNames;
    ExonProbeNbs{2}=AExonProbeNbs;
    GroupRanks{2}=AGroupRanks;
    GroupProbeNbs{2}=AGroupProbeNbs;
    TargetedTs{2}=ATargetedTs;
    NotTargetedTs{2}=ANotTargetedTs;
    PsInfo{2}=APsInfo;
    clear AGeneName AGeneNames ATargetedGenes ATargetingPsRanks APsNames APsRanks AExonNames AExonProbeNbs AGroupRanks AGroupProbeNbs ATargetedTs ANotTargetedTs APsInfo
    clear GeneName GeneNames PsNames PsRanks TargetedTs NotTargetedTs
end


%construct different samples of pairs of probesets that target the same gene(s) with make_pspairs function
SaveDir=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
if ~exist(SaveDir,'dir')
    mkdir(SaveDir)
end
if TestFlag
    SaveFile=sprintf('m%u_dup_probenb%u_single.mat',ChipRank,ProbeNbLimit);
else
    SaveFile=sprintf('m%u_dup_probenb%u.mat',ChipRank,ProbeNbLimit);
end
cd(SaveDir)
if ~exist(SaveFile,'file')
    if TestFlag==0
        SingleFlag=[];        
        [DupStat,DupRank,GeneRankDuplicate,EnsDuplicateOut,AceDuplicateOut,EnsGeneNameOut,AceGeneNameOut,Duplicate,DuplicateOut,DuplicateLowHigh,DuplicateLow,DuplicateHigh]=make_pspairs(ProbeNbLimit,TargetedGenes,TargetingPsRanks,TestFlag,SingleFlag,PsInfo,ProbeNb,AceFlag);
        if AceFlag
            eval(sprintf('save %s DupStat DupRank GeneRankDuplicate EnsDuplicateOut AceDuplicateOut EnsGeneNameOut AceGeneNameOut Duplicate DuplicateOut DuplicateHigh',SaveFile))
        else
            eval(sprintf('save %s DupStat DupRank GeneRankDuplicate EnsDuplicateOut EnsGeneNameOut Duplicate DuplicateOut DuplicateHigh',SaveFile))
        end
    else
        %probe set targeting a single gene
        SaveFile=sprintf('m%u_dup_probenb%u_single.mat',ChipRank,ProbeNbLimit);
        SingleFlag=1;
        [DupStat,DupRank,GeneRankDuplicate,EnsDuplicateOut,AceDuplicateOut,EnsGeneNameOut,AceGeneNameOut,Duplicate,DuplicateOut,DuplicateLowHigh,DuplicateLow,DuplicateHigh]=make_pspairs(ProbeNbLimit,TargetedGenes,TargetingPsRanks,TestFlag,SingleFlag,PsInfo,ProbeNb,AceFlag);
        eval(sprintf('save %s DupStat DupRank GeneRankDuplicate EnsDuplicateOut EnsGeneNameOut Duplicate DuplicateOut DuplicateLowHigh DuplicateLow DuplicateHigh',SaveFile))        
        %probe set targeting severalgenes
        SaveFile=sprintf('m%u_dup_probenb%u_multiple.mat',ChipRank,ProbeNbLimit);
        SingleFlag=0;
        [DupStat,DupRank,GeneRankDuplicate,EnsDuplicateOut,AceDuplicateOut,EnsGeneNameOut,AceGeneNameOut,Duplicate,DuplicateOut,DuplicateLowHigh,DuplicateLow,DuplicateHigh]=make_pspairs(ProbeNbLimit,TargetedGenes,TargetingPsRanks,TestFlag,SingleFlag,PsInfo,ProbeNb,AceFlag);
        eval(sprintf('save %s DupStat DupRank GeneRankDuplicate EnsDuplicateOut EnsGeneNameOut Duplicate DuplicateOut DuplicateLowHigh DuplicateLow DuplicateHigh',SaveFile))
    end
    clear TargetedGenes TargetingPsRanks PsInfo
end
%calculate Pearson correlation coefficient
if TestFlag
    SaveFile=sprintf('m%u_pearson_probenb%u_single.mat',ChipRank,ProbeNbLimit);
else
    if NotFoundFlag
        SaveFile=sprintf('m%u_pearson_probenb%u_notfound.mat',ChipRank,ProbeNbLimit);
    else
        SaveFile=sprintf('m%u_pearson_probenb%u.mat',ChipRank,ProbeNbLimit);
    end
end
cd(SaveDir)
if ~exist(SaveFile,'file')
    calculate_pearson(Species,ChipRank,ProbeNbLimit,TestFlag,NotFoundFlag,IdemFlag)
end

%calculate p-values on neighbourhood
if TestFlag
    for SingleL=1:2
        %for SingleL=2
        if SingleL==1
            SaveFile=sprintf('m%u_dup_probenb%u_single.mat',ChipRank,ProbeNbLimit);
            MultipleFlag=0;
        else
            SaveFile=sprintf('m%u_dup_probenb%u_multiple.mat',ChipRank,ProbeNbLimit);
            MultipleFlag=1;
        end
        load(SaveFile)

        for NetL=1:length(NetRankList)
            NetRank=NetRankList(NetL);

            if ~isempty(DuplicateLowHigh)
                DupType=1;
                NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,DuplicateLowHigh,TestFlag,MultipleFlag,DupType,PvCorrLimits,IdemFlag)
            end
            if ~isempty(DuplicateOut)
                DupType=2;
                NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,DuplicateOut,TestFlag,MultipleFlag,DupType,PvCorrLimits,IdemFlag)
            end
            if ~isempty(DuplicateLow)
                DupType=3;
                NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,DuplicateLow,TestFlag,MultipleFlag,DupType,PvCorrLimits,IdemFlag)
            end
            if ~isempty(DuplicateHigh)
                DupType=4;
                NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,DuplicateHigh,TestFlag,MultipleFlag,DupType,PvCorrLimits,IdemFlag)
            end
            DupType=0;
            NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,Duplicate,TestFlag,MultipleFlag,DupType,PvCorrLimits,IdemFlag)
        end
        
    end
else
    if NotFoundFlag==0
        SaveFile=sprintf('m%u_dup_probenb%u.mat',ChipRank,ProbeNbLimit);
        load(SaveFile)
        
        for NetL=1:length(NetRankList)
            NetRank=NetRankList(NetL);

            if ~isempty(DuplicateOut{1})
                DupType=1;
                NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,DuplicateOut{1},TestFlag,[],DupType,PvCorrLimits,IdemFlag)
            end
            if ~isempty(DuplicateOut{2})
                DupType=2;
                NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,DuplicateOut{2},TestFlag,[],DupType,PvCorrLimits,IdemFlag)
            end
            DupType=0;
            NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,Duplicate,TestFlag,[],DupType,PvCorrLimits,IdemFlag)
        end
        
    else
        SaveFile=sprintf('m%u_dup_probenb%u_notfound.mat',ChipRank,ProbeNbLimit);
        load(SaveFile)
        
        for NetL=1:length(NetRankList)
            NetRank=NetRankList(NetL);

            DupType=3;
            NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,NotFound,TestFlag,[],DupType,PvCorrLimits,IdemFlag)
        end
        
    end
end


%% NODESIM
function NODESIM(Species,ChipRank,PsNb,NetRank,ProbeNbLimit,Dup,TestFlag,MultipleFlag,DupType,PvCorrLimits,IdemFlag)

%         % load information on probesets by gene
%         % for a given gene are listed all the probesets that have a given
%         % number of probes inside its exons
%         %also are given targeted and not targeted list of transcripts
%


global K

if IdemFlag==0
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)))
    eval(sprintf('load m%u_ps2net',ChipRank))
    NetPs=Ps2Net;
    clear Ps2Net
else
    NetPs=[1:PsNb]';
end
%reconstruct [ori] list for multiple (indicates only one probeset of the first chipset in case it is multiple)
ColIndex=1:PsNb;

%save directory and save file

SaveDir=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(SaveDir)
if TestFlag
    if MultipleFlag==0
        if DupType==0
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_single.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==1
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_single_testlowhigh.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==2
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_single_testout.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==3
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_single_testlow.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==4
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_single_testhigh.mat',ChipRank,NetRank,ProbeNbLimit);
        end
    else
        if DupType==0
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_multiple.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==1
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testlowhigh.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==2
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testout.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==3
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testlow.mat',ChipRank,NetRank,ProbeNbLimit);
        elseif DupType==4
            SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_multiple_testhigh.mat',ChipRank,NetRank,ProbeNbLimit);
        end
    end
else
    if DupType==0
        SaveFile=sprintf('m%u_n%u_nodesim_probenb%u.mat',ChipRank,NetRank,ProbeNbLimit);
    elseif DupType==1
        SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_testout_ensembl.mat',ChipRank,NetRank,ProbeNbLimit);
    elseif DupType==2
        SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_testout_aceview.mat',ChipRank,NetRank,ProbeNbLimit);
    elseif DupType==3
        SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_notfound.mat',ChipRank,NetRank,ProbeNbLimit);
    elseif DupType==4
        SaveFile=sprintf('m%u_n%u_nodesim_probenb%u_testhigh.mat',ChipRank,NetRank,ProbeNbLimit);
    end
end
if ~exist(SaveFile,'file')

    NetDir=fullfile(K.dir.net,sprintf('m%03u',ChipRank),sprintf('n%05u',NetRank));
    %nb of common links, links in first chip and links in second chip (or duplicate)
    %in various conditions are stored in Sim
    Sim=[];
    % CORR and ANTI values are stored
    AntiVal=[];
    CorrVal=[];

    cd(NetDir)
    cfid=fopen(sprintf('c_m%u_n%u.4mat',ChipRank,NetRank),'rb');
    afid=fopen(sprintf('a_m%u_n%u.4mat',ChipRank,NetRank),'rb');
    for DupL=1:length(Dup)
        A={};
        C={};
        %PsRanks are unique and ordered in Dup{DupL}
        for PsL=1:length(Dup{DupL})
            fseek(afid,(NetPs(Dup{DupL}(PsL))-1)*PsNb,'bof');
            fseek(cfid,(NetPs(Dup{DupL}(PsL))-1)*PsNb,'bof');
            A{PsL}=uint8(fread(afid,[1,PsNb],'uint8'));
            C{PsL}=uint8(fread(cfid,[1,PsNb],'uint8'));
            A{PsL}=A{PsL}(ColIndex);
            C{PsL}=C{PsL}(ColIndex);
            A{PsL}(isnan(A{PsL}))=0;
            C{PsL}(isnan(C{PsL}))=0;
        end

        for PsL1=1:length(C)-1
            for PsL2=PsL1+1:length(C)
                %do not record doublons
                if isempty(Sim)|isempty(find(Sim(:,1)==Dup{DupL}(PsL1)&Sim(:,2)==Dup{DupL}(PsL2)))
                    Pos2=find(ColIndex==Dup{DupL}(PsL2));
                    AntiVal=[AntiVal;A{PsL1}(Pos2)];
                    CorrVal=[CorrVal;C{PsL1}(Pos2)];

                    Res=[Dup{DupL}(PsL1),Dup{DupL}(PsL2)];
                    %
                    for Limit=PvCorrLimits
                        Pos1=find(C{PsL1}>Limit&C{PsL2}>Limit);
                        Pos2=find(C{PsL1}>Limit);
                        Pos3=find(C{PsL2}>Limit);
                        if C{PsL1}(Pos2)>Limit
                            %without weight
                            SimVal=[length(Pos1),length(Pos2),length(Pos3),(length(Pos1)-1)/min(length(Pos2),length(Pos3))];
                        else
                            SimVal=[length(Pos1),length(Pos2),length(Pos3),(length(Pos1))/(min(length(Pos2),length(Pos3))+1)];
                            %without weight
                        end
                        Res=[Res,SimVal];
                        SimVal=SimVal(1:3);
                        if SimVal(2)==0|SimVal(3)==0
                            Res=[Res,0];
                        else
                            Res=[Res,hypergeometric(SimVal(1),SimVal(2),SimVal(3),PsNb)];
                        end
                    end
                    Sim=[Sim;Res];
                end
            end
        end
    end
    fclose(afid);
    fclose(cfid);

    %reorder results on origin probeset (first column)
    [Temp SortOrder]=sort(Sim(:,1));

    Sim=Sim(SortOrder,:);
    CorrVal=CorrVal(SortOrder);
    AntiVal=AntiVal(SortOrder);
    %create structured variable

    NewSim.corr=CorrVal;
    NewSim.anti=AntiVal;
    NewSim.firstPsRank=Sim(:,1);
    NewSim.sndPsRank=Sim(:,2);
    NewSim.pvCorrLimit=PvCorrLimits;
    for i=1:length(PvCorrLimits)
        NewSim.commonNodeNb{i}=Sim(:,3+(i-1)*5);
        NewSim.firstNodeNb{i}=Sim(:,4+(i-1)*5);
        NewSim.sndNodeNb{i}=Sim(:,5+(i-1)*5);
        NewSim.overlap{i}=Sim(:,6+(i-1)*5);
        NewSim.pv{i}=Sim(:,7+(i-1)*5);
    end
    Sim=NewSim;
    cd(SaveDir)
    eval(sprintf('save %s Sim;',SaveFile))
end
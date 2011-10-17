%============================%
% FUNCTION CALCULATE_PEARSON %
%============================%

% CALCULAT_PEARSON calculates the Pearson's correlation coefficient on all pairs of probe
% sets referenced in dup files
 
%INPUT PARAMETERS
% 1     ModelRank:
% 2 ProbeNbLimit:


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

function calculate_pearson(Species,ChipRank,ProbeNbLimit,TestFlag,NotFoundFlag,IdemFlag)
global K

PsDir=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(PsDir)
ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
PsNb=K.chip.probesetNb(ChipPos);

if IdemFlag==0
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)))
    eval(sprintf('load m%u_ps2net',ChipRank))
    NetPs=Ps2Net;
    clear Ps2Net
else
    NetPs=[1:PsNb]';
end

if TestFlag
    RoundNb=2;
else
    RoundNb=1;
end

for RoundL=1:RoundNb
    if RoundL==1
        if TestFlag
            LoadFile=sprintf('m%u_dup_probenb%u_single.mat',ChipRank,ProbeNbLimit);
            SaveFile=sprintf('m%u_pearson_probenb%u_single.mat',ChipRank,ProbeNbLimit);
        else
            if NotFoundFlag
                LoadFile=sprintf('m%u_dup_probenb%u_notfound.mat',ChipRank,ProbeNbLimit);
                SaveFile=sprintf('m%u_pearson_probenb%u_notfound.mat',ChipRank,ProbeNbLimit);
            else
                LoadFile=sprintf('m%u_dup_probenb%u.mat',ChipRank,ProbeNbLimit);
                SaveFile=sprintf('m%u_pearson_probenb%u.mat',ChipRank,ProbeNbLimit);
            end
        end
    else
        LoadFile=sprintf('m%u_dup_probenb%u_multiple.mat',ChipRank,ProbeNbLimit);
        SaveFile=sprintf('m%u_pearson_probenb%u_multiple.mat',ChipRank,ProbeNbLimit);
    end        
    load(LoadFile)
    
    Pearson=cell(1,5);
    for DupL=1:5
        switch DupL
            case 1
                if NotFoundFlag
                    Dup='NotFound';
                else
                    Dup='Duplicate';
                end
            case 2
                Dup='DuplicateOut';
            case 3
                Dup='DuplicateHigh';
            case 4
                Dup='DuplicateLowHigh';
            case 5
                Dup='DuplicateLow';
        end
        if exist(Dup,'var')
            eval(sprintf('CurrDup=%s;',Dup));
            if ~isempty(CurrDup)
                if iscell(CurrDup{1})
                    [Pearson{DupL}{1}.firstPsRank,Pearson{DupL}{1}.sndPsRank,Pearson{DupL}{1}.rankCorr,Pearson{DupL}{1}.signalCorr]=PEARSON(ChipRank,CurrDup{1},NetPs);
                    [Pearson{DupL}{2}.firstPsRank,Pearson{DupL}{2}.sndPsRank,Pearson{DupL}{2}.rankCorr,Pearson{DupL}{2}.signalCorr]=PEARSON(ChipRank,CurrDup{2},NetPs);
                else
                    [Pearson{DupL}.firstPsRank,Pearson{DupL}.sndPsRank,Pearson{DupL}.rankCorr,Pearson{DupL}.signalCorr]=PEARSON(ChipRank,CurrDup,NetPs);
                end
            else
                Pearson{DupL}.firstPsRank=[];
                Pearson{DupL}.sndPsRank=[];
                Pearson{DupL}.rankCorr=[];
                Pearson{DupL}.signalCorr=[];
            end
        end
    end
    cd(PsDir)
    eval(sprintf('save %s Pearson',SaveFile))
end

function [FirstPsRank,SndPsRank,RankCorr,SignalCorr]=PEARSON(ChipRank,Dup,NetPs)
global K
DataDir=fullfile(K.dir.table,sprintf('m%u',ChipRank),sprintf('m%u_data',ChipRank));
cd(DataDir)
load DataRanks
FirstPsRank=[];
SndPsRank=[];
RankCorr=[];
SignalCorr=[];
for DupL=1:length(Dup)
    %PsRanks are unique and ordered in Dup{DupL}
    for PsL1=1:length(Dup{DupL})-1
        for PsL2=PsL1+1:length(Dup{DupL})
            if isempty(FirstPsRank)|isempty(find(FirstPsRank==Dup{DupL}(PsL1)&SndPsRank==Dup{DupL}(PsL2)))
                CurrFirstPsRank=Dup{DupL}(PsL1);
                CurrSndPsRank=Dup{DupL}(PsL2);           
                FirstPsRank(end+1,1)=CurrFirstPsRank;
                SndPsRank(end+1,1)=CurrSndPsRank;              
                FirstRanks=DataRanks(NetPs(CurrFirstPsRank),:);
                SndRanks=DataRanks(NetPs(CurrSndPsRank),:);
                NegPos=find(FirstRanks<=0|SndRanks<=0);
                FirstRanks(NegPos)=[];
                SndRanks(NegPos)=[];
                CorrCoef=corrcoef(FirstRanks,SndRanks);
                RankCorr(end+1,1)=CorrCoef(2);
                CorrCoef=corrcoef(log2(interp1(K.chip.ref.rank,K.chip.ref.signal,FirstRanks)),log2(interp1(K.chip.ref.rank,K.chip.ref.signal,SndRanks)));
                SignalCorr(end+1,1)=CorrCoef(2);                
            end
        end
    end
end

%reorder results on first probeset
[FirstPsRank SortOrder]=sort(FirstPsRank);
SndPsRank=SndPsRank(SortOrder,:);
RankCorr=RankCorr(SortOrder,:);
SignalCorr=SignalCorr(SortOrder,:);



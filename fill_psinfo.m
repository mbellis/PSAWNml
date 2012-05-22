%======================%
% FUNCTION FILL_PSINFO %
%======================%

%FILL_PSINFO: for a given nb of probes (ProbeNbLimit) recover for each probeset information
% on genes targeted by more or equal ProbeNbLimit probes and
% on genes targeted by less than ProbeNbLimit probes
% Then for each couple of probe sets referenced in Sim files, calculates
% common and uncommon quantities (e.g. the nb of probes in common exons and
% the number of probes in uncommon exons)

%INPUT PARAMETERS
% 1       PsInfo: probe set information filled by import_targetinfo
% 2 ProbeNbLimit: minimum number of tageting probes
%   Statistics on different types of paired probe sets
% 3          Sim: same gene(s) targeted inside exons
% 4       OutSim: same gene(s) targeted outside exons
% 5        LHSim: one probe set targeting a gene or multiple genes with a low number of
%                 probes and the other
%                 probe set targeting another gene(s) with a high number (random pairs)
% 6         LSim: probe set targeting different genes with a low number of probes
%                 (random pairs)
% 7         HSim: probe set targeting diferent genes with a high number of probes
%                 (random pairs
% 8  AceviewFlag: if = 1 process Ensembl and AceView genes; if = 0 process only Ensembl genes. 
% 9      GopFlag: if = 1 pairs of probe sets in Sim that target only one genes are tested to
%                 see if the gene is a group of probes (GOP)
% 

%OUTPUT PARAMETERS
% 1 DupInfo is a structure wich keep information about pairs of probe sets tested
%           and has the following fields for each Sim type (DupInfo{Sim}:
%                 psRank1: first probe set rank
%                 psRank2: second probe set rank
%   Information on genes targeted inside exons
%               comGeneIn: list of commonly targeted genes
%           meanComGeneIn: geometric mean of common genes relative to the number of genes
%                          targeted by each probe set
%             uncomGeneIn: list of genes targeted by only one probe set
%            comTscriptIn: list of commonly targeted transcripts
%           meanTscriptIn: geometric mean of common transcripts relative to the number of
%                          transcripts targeted by each probe set
%          uncomTscriptIn: list of transcripts targeted by only one probe set
%             maxProbe1In: greatest number of targeting probe for the first probe set
%             maxProbe2In: greatest number of targeting probe for the second probe set
%             minProbe1In: least number of targeting probe for the first probe set
%             minProbe2In: least number of targeting probe for the second probe set
%       cMeanGroupProbeIn: geometric mean of common probe nb relative to the number of probe
%                          targeting common exons in each probe set
%       tMeanGroupProbeIn: geometric mean of common probe nb relative to the number of probe
%                          targeting all exons in each probe set
%   Information on genes targeted outside exons (same fields with Out in place of In)
%               comGeneOut .... tMeanGroupProbeOut
%                   isGop: if SingleFlag==1, indicates if the single targeted gene is
%                          a group of probes (GOP: colocalized probes, but no gene described)

% 2 Ps is a structure with the following fields (Ps{Type}{PsRank} with Type=1 for Ensembl
%                                               Type=2 for AceView)
%   Information on genes targeted by more than or equal to ProbeNbLimit probes:
%         geneNamesSup: name of targeted genes
%        groupRanksSup: rank of targeted group of transcripts
%     groupProbeNbsSup: number of targeting probes in each group
%       transcriptsSup: list of targeted transcripts
%    notTranscriptsSup: list of not targeted transcripts
%           probeNbSup: number of targeting probes in each gene
%      geneNamesOutSup: name of genes targeted outside exons
%        probeNbOutSup: number of targeting probes in each outside gene
%   Information on genes targeted by less than ProbeNbLimit probes:
%         geneNamesInf: name of targeted genes
%        groupRanksInf: rank of targeted group of transcripts
%     groupProbeNbsInf: number of targeting probes in each group
%       transcriptsInf: list of targeted transcripts
%    notTranscriptsInf: list of not targeted transcripts
%           probeNbInf: number of targeting probes in each gene
%      geneNamesOutInf: name of genes targeted outside exons

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



function [DupInfo,Ps]=fill_psinfo(PsInfo,ProbeNbLimit,Sim,OutSim,HSim,LHSim,LSim,AceviewFlag,GopFlag)
%% FILL PS


if AceviewFlag==0
    %process only Ensembl genes
    TypeEnd=1;
else
    %process Ensembl and AceView genes
    TypeEnd=2;
end
PsNb=length(PsInfo{1}.geneNames);
Ps=cell(2,1);
for TypeL=1:TypeEnd
    for PsL=1:PsNb
        Ps{TypeL}{PsL}.geneNamesSup={};
        Ps{TypeL}{PsL}.groupRanksSup={};
        Ps{TypeL}{PsL}.groupProbeNbsSup={};
        Ps{TypeL}{PsL}.transcriptsSup={};
        Ps{TypeL}{PsL}.notTranscriptsSup={};
        Ps{TypeL}{PsL}.probeNbSup=[];
        Ps{TypeL}{PsL}.geneNamesOutSup={};
        Ps{TypeL}{PsL}.probeNbOutSup=[];
        Ps{TypeL}{PsL}.geneNamesInf={};
        Ps{TypeL}{PsL}.groupRanksInf={};
        Ps{TypeL}{PsL}.groupProbeNbsInf={};
        Ps{TypeL}{PsL}.transcriptsInf={};
        Ps{TypeL}{PsL}.notTranscriptsInf={};
        Ps{TypeL}{PsL}.probeNbInf=[];
        Ps{TypeL}{PsL}.geneNamesOutInf={};
    end
end

%process Ensembl and eventually AceView Genes
for TypeL=1:TypeEnd
    %maximal number of probes in a probe set
    MaxProbeNb=size(PsInfo{TypeL}.geneNames{1},1)-1;
    % fill info for each probesets
    for PsL=1:PsNb
        %information on genes targeted by more or equal ProbeNbLimit probes
        GeneNames={};
        GroupRanks={};
        GroupProbeNbs={};
        Transcripts={};
        NotTranscripts={};
        ProbeNb=[];
        %for PsInfo name ProbeNb is shifted (+1)
        for ProbeNbL=ProbeNbLimit+1:MaxProbeNb+1            
                GeneNameNb=length(PsInfo{TypeL}.geneNames{PsL}{ProbeNbL});
                if GeneNameNb>0
                    GeneNames=[GeneNames,PsInfo{TypeL}.geneNames{PsL}{ProbeNbL}];
                    GroupRanks=[GroupRanks,PsInfo{TypeL}.groupRanks{PsL}{ProbeNbL}];
                    GroupProbeNbs=[GroupProbeNbs,PsInfo{TypeL}.groupProbeNbs{PsL}{ProbeNbL}];
                    Transcripts=[Transcripts,PsInfo{TypeL}.transcripts{PsL}{ProbeNbL}];
                    NotTranscripts=[NotTranscripts,PsInfo{TypeL}.notTranscripts{PsL}{ProbeNbL}];
                    ProbeNb=[ProbeNb,repmat(ProbeNbL-1,1,GeneNameNb)];
                end
        end
        Ps{TypeL}{PsL}.geneNamesSup=GeneNames;
        Ps{TypeL}{PsL}.groupRanksSup=GroupRanks;
        Ps{TypeL}{PsL}.groupProbeNbsSup=GroupProbeNbs;
        Ps{TypeL}{PsL}.transcriptsSup=Transcripts;
        Ps{TypeL}{PsL}.notTranscriptsSup=NotTranscripts;
        Ps{TypeL}{PsL}.probeNbSup=ProbeNb;

        %probesets targetting only genes in up, down or intron with equal or more than ProbeNbLimit probes
        %may target genes in exon or splice but with less than ProbeNbLimit probes
        Pos=find(PsInfo{TypeL}.notInExonProbeNbs{PsL}{1}>=ProbeNbLimit);
        if ~isempty(Pos)
            SelPos=[];
            for PosL=1:length(Pos)
                %don'tlist genes that have already been selected in the
                %previous step
                if isempty(strmatch(PsInfo{TypeL}.geneNames{PsL}{1}{Pos(PosL)},GeneNames,'exact'))
                    SelPos=[SelPos;Pos(PosL)];
                end
            end
            Ps{TypeL}{PsL}.geneNamesOutSup=PsInfo{TypeL}.geneNames{PsL}{1}(SelPos);
            Ps{TypeL}{PsL}.probeNbOutSup=PsInfo{TypeL}.notInExonProbeNbs{PsL}{1}(SelPos);
        end

        if ProbeNbLimit>1
            %information on genes targeted by less than ProbeNbLimit probes
            GeneNames={};
            GroupRanks={};
            GroupProbeNbs={};
            Transcripts={};
            NotTranscripts={};
            ProbeNb=[];

            for ProbeNbL=2:ProbeNbLimit
                GeneNameNb=length(PsInfo{TypeL}.geneNames{PsL}{ProbeNbL});
                if GeneNameNb>0
                    GeneNames=[GeneNames,PsInfo{TypeL}.geneNames{PsL}{ProbeNbL}];
                    GroupRanks=[GroupRanks,PsInfo{TypeL}.groupRanks{PsL}{ProbeNbL}];
                    GroupProbeNbs=[GroupProbeNbs,PsInfo{TypeL}.groupProbeNbs{PsL}{ProbeNbL}];
                    Transcripts=[Transcripts,PsInfo{TypeL}.transcripts{PsL}{ProbeNbL}];
                    NotTranscripts=[NotTranscripts,PsInfo{TypeL}.notTranscripts{PsL}{ProbeNbL}];
                    ProbeNb=[ProbeNb,repmat(ProbeNbL-1,1,GeneNameNb)];
                end
            end
            Ps{TypeL}{PsL}.geneNamesInf=GeneNames;
            Ps{TypeL}{PsL}.groupRanksInf=GroupRanks;
            Ps{TypeL}{PsL}.groupProbeNbsInf=GroupProbeNbs;
            Ps{TypeL}{PsL}.transcriptsInf=Transcripts;
            Ps{TypeL}{PsL}.notTranscriptsInf=NotTranscripts;
            Ps{TypeL}{PsL}.probeNbInf=ProbeNb;

            Pos=find(PsInfo{TypeL}.notInExonProbeNbs{PsL}{1}<ProbeNbLimit&PsInfo{TypeL}.notInExonProbeNbs{PsL}{1}>0);
            if ~isempty(Pos)
                SelPos=[];
                for PosL=1:length(Pos)                   
                    if isempty(strmatch(PsInfo{TypeL}.geneNames{PsL}{1}{Pos(PosL)},GeneNames,'exact'))
                        SelPos=[SelPos;Pos(PosL)];
                    end
                end
                Ps{TypeL}{PsL}.geneNamesOutInf=PsInfo{TypeL}.geneNames{PsL}{1}(SelPos);
            end
        end
    end
end

%% FILL COMPINFO
%fill info for each couple of probesets present in Sim
DupInfo=cell(1,3);
for SimL=1:5
    DupInfo{SimL}=cell(1,2);
    switch SimL
        case 1
            CurrSim=Sim;
        case 2
            CurrSim=OutSim;
        case 3
            CurrSim=HSim;
        case 4
            CurrSim=LHSim;
        case 5
            CurrSim=LSim;
    end
    if ~isempty(CurrSim)
        CompNb=length(CurrSim.firstPsRank);
        for TypeL=1:TypeEnd
            CompInfo.psRank1=zeros(CompNb,1);
            CompInfo.psRank2=zeros(CompNb,1);
            CompInfo.comGeneIn=zeros(CompNb,1);
            CompInfo.meanComGeneIn=zeros(CompNb,1);
            CompInfo.uncomGeneIn=zeros(CompNb,1);
            CompInfo.comTscriptIn=zeros(CompNb,1);
            CompInfo.meanTscriptIn=zeros(CompNb,1);
            CompInfo.uncomTscriptIn=zeros(CompNb,1);
            CompInfo.maxProbe1In=zeros(CompNb,1);
            CompInfo.maxProbe2In=zeros(CompNb,1);
            CompInfo.minProbe1In=zeros(CompNb,1);
            CompInfo.minProbe2In=zeros(CompNb,1);
            CompInfo.cMeanGroupProbeIn=zeros(CompNb,1);
            CompInfo.tMeanGroupProbeIn=zeros(CompNb,1);
            CompInfo.comGeneOut=zeros(CompNb,1);
            CompInfo.meanComGeneOut=zeros(CompNb,1);
            CompInfo.uncomGeneOut=zeros(CompNb,1);
            CompInfo.comTscriptOut=zeros(CompNb,1);
            CompInfo.meanTscriptOut=zeros(CompNb,1);
            CompInfo.uncomTscriptOut=zeros(CompNb,1);
            CompInfo.maxProbe1Out=zeros(CompNb,1);
            CompInfo.maxProbe2Out=zeros(CompNb,1);
            CompInfo.minProbe1Out=zeros(CompNb,1);
            CompInfo.minProbe2Out=zeros(CompNb,1);
            CompInfo.tMeanGroupProbeOut=zeros(CompNb,1);
            CompInfo.cMeanGroupProbeOut=zeros(CompNb,1);
            CompInfo.isGop=zeros(CompNb,1);

            for PsL=1:CompNb
                %info for genes targeted with more than ProbNbLimit - 1 probes
                PsRank1=CurrSim.firstPsRank(PsL);
                PsRank2=CurrSim.sndPsRank(PsL);
                CompInfo.psRank1(PsL)=PsRank1;
                CompInfo.psRank2(PsL)=PsRank2;

                if SimL==2
                    %list of common genes
                    [CommonGenes,Index1,Index2]=intersect(Ps{TypeL}{PsRank1}.geneNamesOutSup,Ps{TypeL}{PsRank2}.geneNamesOutSup);
                    %geometric mean of common genes relative to the number of genes targeted by each probe set
                    CompInfo.meanComGeneIn(PsL)=length(CommonGenes)/sqrt(length(Ps{TypeL}{PsRank1}.geneNamesOutSup)*length(Ps{TypeL}{PsRank2}.geneNamesOutSup));
                    %nb of uncommon genes
                    CompInfo.uncomGeneIn(PsL)=length(setxor(Ps{TypeL}{PsRank1}.geneNamesOutSup,Ps{TypeL}{PsRank2}.geneNamesOutSup));
                else
                    [CommonGenes,Index1,Index2]=intersect(Ps{TypeL}{PsRank1}.geneNamesSup,Ps{TypeL}{PsRank2}.geneNamesSup);
                    CompInfo.meanComGeneIn(PsL)=length(CommonGenes)/sqrt(length(Ps{TypeL}{PsRank1}.geneNamesSup)*length(Ps{TypeL}{PsRank2}.geneNamesSup));                    
                    CompInfo.uncomGeneIn(PsL)=length(setxor(Ps{TypeL}{PsRank1}.geneNamesSup,Ps{TypeL}{PsRank2}.geneNamesSup));
                end
                %nb of common genes
                CompInfo.comGeneIn(PsL)=length(CommonGenes);

                if GopFlag &(SimL==1|SimL==2)
                    if length(CommonGenes)==1
                        if findstr('GOP',CommonGenes{1})>0
                            CompInfo.isGop(PsL)=1;
                        end
                    end
                end

                %nb of targeted transcripts & nb of common transcripts nb of  common probes in common exons
                if SimL==1;
                    %total number of probes in all exons for each probe set
                    MaxProbeNb=zeros(1,2);
                    %total number of probes in common exons for each probe set
                    TCProbeNb=zeros(1,2);
                    %number of probes in common in common exons (minimum of numbers of probes)
                    CProbeNb=0;
                    %total number of transcripts for each probe set
                    TTscriptNb=zeros(1,2);
                    %total number of transcripts for both probe sets
                    TscriptNb=0;
                    %number of transcripts in common to both probe sets
                    CTscriptNb=0;
                    if CompInfo.comGeneIn(PsL)>0
                        for i=1:CompInfo.comGeneIn(PsL)
                            %total number of probes in all exons for probe set 1
                            MaxProbeNb(1)=MaxProbeNb(1)+sum(Ps{TypeL}{PsRank1}.groupProbeNbsSup{Index1(i)});
                            %total number of probes in all exons for probe set 2
                            MaxProbeNb(2)=MaxProbeNb(2)+sum(Ps{TypeL}{PsRank2}.groupProbeNbsSup{Index2(i)});
                            %list of common exons and respective index for each probe set
                            [CommonGroups,GIndex1,GIndex2]=intersect(Ps{TypeL}{PsRank1}.groupRanksSup{Index1(i)},Ps{TypeL}{PsRank2}.groupRanksSup{Index2(i)});
                            %total number of probes in common exons for probe set 1
                            TCProbeNb(1)=TCProbeNb(1)+sum(Ps{TypeL}{PsRank1}.groupProbeNbsSup{Index1(i)}(GIndex1));
                            %total number of probes in common exons for probe set 2
                            TCProbeNb(2)=TCProbeNb(2)+sum(Ps{TypeL}{PsRank2}.groupProbeNbsSup{Index2(i)}(GIndex2));
                            %number of probes in common in common exons                             
                            if ~isempty(Ps{TypeL}{PsRank1}.groupProbeNbsSup{Index1(i)}(GIndex1))
                                CProbeNb=CProbeNb+min(sum(Ps{TypeL}{PsRank1}.groupProbeNbsSup{Index1(i)}(GIndex1)),sum(Ps{TypeL}{PsRank2}.groupProbeNbsSup{Index2(i)}(GIndex2)));
                            end
                            %total number of transcripts for probe set 1
                            TTscriptNb(1)=TTscriptNb(1)+length(Ps{TypeL}{PsRank1}.transcriptsSup{Index1(i)});
                            %total number of transcripts for probe set 2
                            TTscriptNb(2)=TTscriptNb(2)+length(Ps{TypeL}{PsRank2}.transcriptsSup{Index2(i)});
                            %total number of transcripts for both probe sets
                            TscriptNb=TscriptNb+length(union(Ps{TypeL}{PsRank1}.transcriptsSup{Index1(i)},Ps{TypeL}{PsRank2}.transcriptsSup{Index2(i)}));
                            %number of transcripts in common to both probe sets
                            CTscriptNb=CTscriptNb+length(intersect(Ps{TypeL}{PsRank1}.transcriptsSup{Index1(i)},Ps{TypeL}{PsRank2}.transcriptsSup{Index2(i)}));
                        end
                        %geometric mean of common probe nb relative to the number of probe targeting all exons in each probe set
                        CompInfo.tMeanGroupProbeIn(PsL)=CProbeNb/sqrt(MaxProbeNb(1)*MaxProbeNb(2));
                        %geometric mean of common probe nb relative to the number of probe targeting common exons in each probe set
                        CompInfo.cMeanGroupProbeIn(PsL)=CProbeNb/sqrt(TCProbeNb(1)*TCProbeNb(2));
                        CompInfo.comTscriptIn(PsL)=CTscriptNb;
                        %geometric mean of common transcript nb relative to the number of all transcripts targeted by each probe set
                        CompInfo.meanTscriptIn(PsL)=CTscriptNb/sqrt(TTscriptNb(1)*TTscriptNb(2));
                        CompInfo.uncomTscriptIn(PsL)=TscriptNb-CTscriptNb;
                    end

                    if ~isempty(Ps{TypeL}{PsRank1}.probeNbSup)
                        CompInfo.maxProbe1In(PsL)=max(Ps{TypeL}{PsRank1}.probeNbSup);
                        CompInfo.minProbe1In(PsL)=min(Ps{TypeL}{PsRank1}.probeNbSup);
                    end
                    if ~isempty(Ps{TypeL}{PsRank2}.probeNbSup)
                        CompInfo.maxProbe2In(PsL)=max(Ps{TypeL}{PsRank2}.probeNbSup);
                        CompInfo.minProbe2In(PsL)=min(Ps{TypeL}{PsRank2}.probeNbSup);
                    end
                elseif SimL==2
                    if ~isempty(Ps{TypeL}{PsRank1}.probeNbOutSup)
                        CompInfo.maxProbe1In(PsL)=max(Ps{TypeL}{PsRank1}.probeNbOutSup);
                        CompInfo.minProbe1In(PsL)=min(Ps{TypeL}{PsRank1}.probeNbOutSup);
                    end
                    if ~isempty(Ps{TypeL}{PsRank2}.probeNbOutSup)
                        CompInfo.maxProbe2In(PsL)=max(Ps{TypeL}{PsRank2}.probeNbOutSup);
                        CompInfo.minProbe2In(PsL)=min(Ps{TypeL}{PsRank2}.probeNbOutSup);
                    end
                end

                if ProbeNbLimit>1
                    %info for genes targeted with less than ProbNbLimit probes
                    %nb of common genes
                    [CommonGenes,Index1,Index2]=intersect(Ps{TypeL}{PsRank1}.geneNamesInf,Ps{TypeL}{PsRank2}.geneNamesInf);
                    CompInfo.meanComGeneOut(PsL)=length(CommonGenes)/sqrt(length(Ps{TypeL}{PsRank1}.geneNamesInf)*length(Ps{TypeL}{PsRank2}.geneNamesInf));
                    CompInfo.comGeneOut(PsL)=length(CommonGenes);
                    %nb of uncommon genes
                    CompInfo.uncomGeneOut(PsL)=length(setxor(Ps{TypeL}{PsRank1}.geneNamesInf,Ps{TypeL}{PsRank2}.geneNamesInf));
                    %nb of targeted transcripts & nb of common transcripts
                    MaxProbeNb=zeros(1,2);
                    TCProbeNb=zeros(1,2);
                    CProbeNb=0;
                    TTscriptNb=zeros(1,2);
                    TscriptNb=0;
                    CTscriptNb=0;
                    if CompInfo.comGeneOut(PsL)>0
                        for i=1:CompInfo.comGeneOut(PsL)
                            MaxProbeNb(1)=MaxProbeNb(1)+sum(Ps{TypeL}{PsRank1}.groupProbeNbsInf{Index1(i)});
                            MaxProbeNb(2)=MaxProbeNb(2)+sum(Ps{TypeL}{PsRank2}.groupProbeNbsInf{Index2(i)});
                            [CommonGroups,GIndex1,GIndex2]=intersect(Ps{TypeL}{PsRank1}.groupRanksInf{Index1(i)},Ps{TypeL}{PsRank2}.groupRanksInf{Index2(i)});
                            if ~isempty(Ps{TypeL}{PsRank1}.groupProbeNbsInf{Index1(i)}(GIndex1))
       
                                CProbeNb=CProbeNb+min(sum(Ps{TypeL}{PsRank1}.groupProbeNbsInf{Index1(i)}(GIndex1)),sum(Ps{TypeL}{PsRank2}.groupProbeNbsInf{Index2(i)}(GIndex2)));
                            end
                            TCProbeNb(1)=TCProbeNb(1)+sum(Ps{TypeL}{PsRank1}.groupProbeNbsInf{Index1(i)}(GIndex1));
                            TCProbeNb(2)=TCProbeNb(2)+sum(Ps{TypeL}{PsRank2}.groupProbeNbsInf{Index2(i)}(GIndex2));
                            TTscriptNb(1)=TTscriptNb(1)+length(Ps{TypeL}{PsRank1}.transcriptsInf{Index1(i)});
                            TTscriptNb(2)=TTscriptNb(2)+length(Ps{TypeL}{PsRank2}.transcriptsInf{Index2(i)});
                            TscriptNb=TscriptNb+length(union(Ps{TypeL}{PsRank1}.transcriptsInf{Index1(i)},Ps{TypeL}{PsRank2}.transcriptsInf{Index2(i)}));
                            CTscriptNb=CTscriptNb+length(intersect(Ps{TypeL}{PsRank1}.transcriptsInf{Index1(i)},Ps{TypeL}{PsRank2}.transcriptsInf{Index2(i)}));
                        end
                        CompInfo.tMeanGroupProbeOut(PsL)=CProbeNb/sqrt(MaxProbeNb(1)*MaxProbeNb(2));
                        CompInfo.cMeanGroupProbeOut(PsL)=CProbeNb/sqrt(TCProbeNb(1)*TCProbeNb(2));
                        CompInfo.comTscriptOut(PsL)=CTscriptNb;
                        CompInfo.meanTscriptOut(PsL)=CTscriptNb/sqrt(TTscriptNb(1)*TTscriptNb(2));
                        CompInfo.uncomTscriptOut(PsL)=TscriptNb-CTscriptNb;
                    end

                    if ~isempty(Ps{TypeL}{PsRank1}.probeNbInf)
                        CompInfo.maxProbe1Out(PsL)=max(Ps{TypeL}{PsRank1}.probeNbInf);
                        CompInfo.minProbe1Out(PsL)=min(Ps{TypeL}{PsRank1}.probeNbInf);
                    end
                    if ~isempty(Ps{TypeL}{PsRank2}.probeNbInf)
                        CompInfo.maxProbe2Out(PsL)=max(Ps{TypeL}{PsRank2}.probeNbInf);
                        CompInfo.minProbe2Out(PsL)=min(Ps{TypeL}{PsRank2}.probeNbInf);
                    end
                end
            end
            DupInfo{SimL}{TypeL}=CompInfo;
        end
    end
end
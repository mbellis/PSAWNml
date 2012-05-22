%==============================%
% FUNCTION IMPORT_TARGETINFO   %
%==============================%

%IMPORT_TARGETINFO read a series of text files that indicate for each gene the list of
% probesets that target it, with detailed information about the exons, group of exons
% and transcripts that are targeted.
%
%INPUT PARAMETERS
% ChipRank: chip model rank
%
%EXTERNAL FILES
% Read files 'ensembl_probesets_by_gene_%ProbeNb_m%ChipRank.txt'
% (and eventually files 'aceview_probesets_by_gene_%ProbeNb_m%ChipRank.txt')
% with ProbeNb in range (0,n). File with probe nb 0 lists the genes that are targeted out of
% their exons (i.e. in their introns, or in the 2kb upwards and downwards sequence). n is
% the maximum number of probes found targeting a single gene.
%
% File format: [Gene ID,...
%              list of probesets targeting this gene, 
%              list of corresponding rank,..
%              lists of targeted exons (one for each targeting probeset,...
%              lists of number of probes in each exon,...
%              last group rank,... 
%              list of grouped targeted exons (some exons overlap each other),...
%              number of probes in each group,...
%              list of targeted transcripts,...
%              list of not targeted transcripts,...
%              list of number of probe outside exons
%              list of number of probe inside gene
%
%          ex: ENSMUSG00000000263
%              {'100385_at'}
%              [9079]
%              {{'ENSMUSE00000662385' 'ENSMUSE00000653296' 'ENSMUSE00000307025'}
%              {[6 6 2 6 6]}
%              7 
%              {[8 6 7]}
%              {[6 6 2]}
%              {[1,2,3]}
%              {[]}
%              {[0]}
%              {[15]}
%
%OUTPUT PARAMETERS
% Write 'm%ChipRank_probeset_by_ensembl_gene .mat' (and eventually
% 'm%ChipRank_probeset_by_aceview_gene .mat')
%
% This file contains:
%
% EPsInfo for Ensembl genes (APsInfo for Aceview genes)
% a PsNb x n cell with the following structure:
%
% EPsInfo{PsRank}{ProbeNb}.exonNames
% EPsInfo{PsRank}{ProbeNb}.exonRanks
% EPsInfo{PsRank}{ProbeNb}.exonProbeNbs
% EPsInfo{PsRank}{ProbeNb}.lastExon
% EPsInfo{PsRank}{ProbeNb}.lastGroup
% EPsInfo{PsRank}{ProbeNb}.groupRanks
% EPsInfo{PsRank}{ProbeNb}.groupProbeNbs
% EPsInfo{PsRank}{ProbeNb}.transcripts
% EPsInfo{PsRank}{ProbeNb}.notTranscripts
% EPsInfo{PsRank}{ProbeNb}.inGeneProbeNbs
% EPsInfo{PsRank}{ProbeNb}.notInExonProbeNbs
%
% and the following variables for Ensembl genes (the same prefixed with A for Aceview genes):
%
%         EGeneName: cell(n+1,1) containing the list of Gene ID of each input file
%                    (genes targeted by n- 1 probes)
%        EGeneNames: a list of unique Gene ID contained in EGeneName
%    ETargetedGenes: a PsNbxn matrix containing at position i,j the position of GeneName
%                    in EGeneName{j}
% ETargetingPsRanks: a PsNbxn matrix containing at position i,j the ranks of probesets
%                    that target EGeneName{j} with n-1 probes
%          EPsNames: cell(n+1,1) containing the list of probeset names of each input file
%  EPsRanks=PsRanks: cell(n+1,1) containing the list of probeset ranks of each input file
%        EExonNames: cell(n+1,1) containing the list of exon names of each input file
%        EExonRank: cell(n+1,1) containing the list of exon ranks of each input file
%     EExonProbeNbs: cell(n+1,1) containing the list of probeset names of each input file
%        ELastExon: a PsNbx1 matrix containing the rank of the last exon
%        ELastGroup: a PsNbx1 matrix containing the rank of the last group of exons
%       EGroupRanks: cell(n+1,1) containing the list of grouped exons
%    EGroupProbeNbs: cell(n+1,1) containing the list of number of probes targeting grouped
%                    exons (overlaping exons)
%       ETargetedTs: cell(n+1,1) containing the list of targeted transcripts
%    ENotTargetedTs: cell(n+1,1) containing the list of not targeted transcripts


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



function import_targetinfo(ChipRank)
global K
ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
if isempty(ChipPos)
    h=errordlg(sprintf('chip m%u does not exist',ChipRank));
    waitfor(h)
    error('process canceled')
end
Species=K.chip.species{ChipPos};
ProbeNb=K.chip.probeNb(ChipPos);
NetPsNb=K.chip.probesetNb(ChipPos);
PROBE_NB=25;
for TypeL=1:2
    if TypeL==1
        Type='ensembl';
    else
        Type='aceview';
    end
    ExistFile=0;
    for ProbeNbL=1:PROBE_NB
        if exist(fullfile(K.dir.pydata,Species,'txt',sprintf('%s_m%u_probesets_by_gene_%02u.txt',Type,ChipRank,ProbeNbL-1)),'file')
            ExistFile=1;
            break
        end
    end
    if ExistFile
        cd(fullfile(K.dir.pydata,Species,'txt'))
        %load informations written from python organized by gene and by nb or targeting probesets (ProbeNb)    %
        %real value of ProbeNb equal ProbeNb-1
        %for ProbeNb=1 the data concerns the genes that have probes not in their exons (up, down or intron)
        GeneName=cell(PROBE_NB,1);   
        PsNames=cell(PROBE_NB,1);
        PsRanks=cell(PROBE_NB,1);
        TargetedT=cell(PROBE_NB,1);
        NotTargetedT=cell(PROBE_NB,1);
        InGeneProbeNb=cell(PROBE_NB,1);
        NotInExonProbeNb=cell(PROBE_NB,1);
        for ProbeNbL=1:PROBE_NB
            %cell array of gene names targeted by ProbeNb-1 probes
            GeneName{ProbeNbL}={};
            %cell array of probeset names targeting each gene by ProbeNb probes
            PsNames{ProbeNbL}={};
            %matrix of probeset ranks targeting each gene by ProbeNb probes
            PsRanks{ProbeNbL}={};
            %matrix of transcript ranks targeted by each probeset by ProbeNb probes
            TargetedT{ProbeNbL}={};
            %matrix of transcript ranks not targeted by each probeset by ProbeNb probes
            NotTargetedT{ProbeNbL}={};
            %matrix of nb of probes targeting the gene (up + down + inexon + inintron + insplice)
            InGeneProbeNb{ProbeNbL}={};
            %matrix of nb of probes not targeting the exons (up + down + inintron )
            NotInExonProbeNb{ProbeNbL}={};
            if exist(fullfile(K.dir.pydata,Species,'txt',sprintf('%s_m%u_probesets_by_gene_%02u.txt',Type,ChipRank,ProbeNbL-1)),'file')
                if TypeL==1
                    [GeneName{ProbeNbL},PsNames{ProbeNbL},PsRanks{ProbeNbL},...
                        ExonName{ProbeNbL},ExonRanks{ProbeNbL},ExonProbeNbs{ProbeNbL},LastExon{ProbeNbL},LastGroup{ProbeNbL},GroupRanks{ProbeNbL},GroupProbeNbs{ProbeNbL},...
                        TargetedT{ProbeNbL},NotTargetedT{ProbeNbL},NotInExonProbeNb{ProbeNbL},InGeneProbeNb{ProbeNbL}]=textread(sprintf('%s_m%u_probesets_by_gene_%02u.txt',Type,ChipRank,ProbeNbL-1),'%s%s%s%s%s%s%u%u%s%s%s%s%s%s','delimiter','\t','bufsize',40000);
                    ProbeNb=min(PROBE_NB,ProbeNbL);
                else
                    [GeneName{ProbeNbL},PsNames{ProbeNbL},PsRanks{ProbeNbL},...
                        ExonName{ProbeNbL},ExonProbeNbs{ProbeNbL},LastExon{ProbeNbL},LastGroup{ProbeNbL},GroupRanks{ProbeNbL},GroupProbeNbs{ProbeNbL},...
                        TargetedT{ProbeNbL},NotTargetedT{ProbeNbL},NotInExonProbeNb{ProbeNbL},InGeneProbeNb{ProbeNbL}]=textread(sprintf('%s_m%u_probesets_by_gene_%02u.txt',Type,ChipRank,ProbeNbL-1),'%s%s%s%s%s%u%u%s%s%s%s%s%s','delimiter','\t','bufsize',40000);
                    ProbeNb=min(PROBE_NB,ProbeNbL);
                end
            end
        end
        %transform probeset indexes into ranks
        for ProbeNbL=1:ProbeNb        
            if ~isempty(PsRanks{ProbeNbL})
                for GeneL=1:length(PsRanks{ProbeNbL})
                    PsRanks{ProbeNbL}{GeneL}=eval(PsRanks{ProbeNbL}{GeneL})+1;
                end
            end
        end

        % recover information in matrix form with eval
        % which transform string into cell
        for ProbeNbL=1:ProbeNb
            if ~isempty(ExonProbeNbs{ProbeNbL})
                for GeneL=1:length(ExonProbeNbs{ProbeNbL})
                    ExonProbeNbs{ProbeNbL}{GeneL}=eval(ExonProbeNbs{ProbeNbL}{GeneL});
                end
            end
        end
        
        ExonNames=cell(size(ExonName));
        for ProbeNbL=1:ProbeNb
            if ~isempty(ExonName{ProbeNbL})
                ExonNames{ProbeNbL}=cell(size(ExonName{ProbeNbL}));
                for GeneL=1:length(ExonName{ProbeNbL})
                    ExonNames{ProbeNbL}{GeneL}=eval(ExonName{ProbeNbL}{GeneL});
                end
            end
        end
        
        if TypeL==2
            ExonRanks=cell(size(ExonName));
        end
        for ProbeNbL=1:ProbeNb
            if TypeL==1
                if ~isempty(ExonRanks{ProbeNbL})
                    for GeneL=1:length(ExonRanks{ProbeNbL})
                        ExonRanks{ProbeNbL}{GeneL}=eval(ExonRanks{ProbeNbL}{GeneL});
                    end
                end
            else
                if ~isempty(ExonNames{ProbeNbL})
                    for GeneL=1:length(ExonNames{ProbeNbL})
                        CurrExonNames=ExonNames{ProbeNbL}{GeneL};
                        for PsL=1:length(CurrExonNames)
                            CurrExonRanks=zeros(1,length(CurrExonNames{PsL}));
                            for ExonL=1:length(CurrExonNames{PsL})
                                CurrExonRank=regexp(CurrExonNames{PsL}{ExonL},'(?<=.exon)\d+','match');
                                try
                                    CurrExonRanks(ExonL)=str2num(CurrExonRank{1});
                                catch
                                    CurrExonRanks(ExonL)=-1;
                                end
                            end
                            ExonRanks{ProbeNbL}{GeneL}{PsL}=CurrExonRanks;
                        end
                    end
                else
                    ExonRanks{ProbeNbL}{GeneL}=[];
                end
            end
        end

        for ProbeNbL=1:ProbeNb
            if ~isempty(GroupRanks{ProbeNbL})
                for GeneL=1:length(GroupRanks{ProbeNbL})
                    GroupRanks{ProbeNbL}{GeneL}=eval(GroupRanks{ProbeNbL}{GeneL});
                end
            end
        end

        for ProbeNbL=1:ProbeNb
            if ~isempty(GroupProbeNbs{ProbeNbL})
                for GeneL=1:length(GroupProbeNbs{ProbeNbL})
                    GroupProbeNbs{ProbeNbL}{GeneL}=eval(GroupProbeNbs{ProbeNbL}{GeneL});
                end
            end
        end        

        TargetedTs=cell(size(TargetedT));
        for ProbeNbL=1:ProbeNb
            if ~isempty(TargetedT{ProbeNbL})
                TargetedTs{ProbeNbL}=cell(size(TargetedT{ProbeNbL}));
                for GeneL=1:length(TargetedT{ProbeNbL})
                    TargetedTs{ProbeNbL}{GeneL}=eval(TargetedT{ProbeNbL}{GeneL});
                end
            end
        end

        NotTargetedTs=cell(size(NotTargetedT));
        for ProbeNbL=1:ProbeNb
            if ~isempty(NotTargetedT{ProbeNbL})
                NotTargetedTs{ProbeNbL}=cell(size(NotTargetedT{ProbeNbL}));
                for GeneL=1:length(NotTargetedT{ProbeNbL})
                    NotTargetedTs{ProbeNbL}{GeneL}=eval(NotTargetedT{ProbeNbL}{GeneL});
                end
            end
        end

        InGeneProbeNbs=cell(size(InGeneProbeNb));
        for ProbeNbL=1:ProbeNb
            if ~isempty(InGeneProbeNb{ProbeNbL})
                InGeneProbeNbs{ProbeNbL}=cell(size(InGeneProbeNb{ProbeNbL}));
                for GeneL=1:length(InGeneProbeNb{ProbeNbL})
                    InGeneProbeNbs{ProbeNbL}{GeneL}=eval(InGeneProbeNb{ProbeNbL}{GeneL});
                    InGeneProbeNbs{ProbeNbL}{GeneL}=InGeneProbeNbs{ProbeNbL}{GeneL}{1};
                end
            end
        end

        NotInExonProbeNbs=cell(size(NotInExonProbeNb));
        for ProbeNbL=1:ProbeNb
            if ~isempty(NotInExonProbeNb{ProbeNbL})
                NotInExonProbeNbs{ProbeNbL}=cell(size(NotInExonProbeNb{ProbeNbL}));
                for GeneL=1:length(NotInExonProbeNb{ProbeNbL})
                    NotInExonProbeNbs{ProbeNbL}{GeneL}=eval(NotInExonProbeNb{ProbeNbL}{GeneL});
                    NotInExonProbeNbs{ProbeNbL}{GeneL}=NotInExonProbeNbs{ProbeNbL}{GeneL}{1};
                end
            end
        end

        %synthetize results by indexing on genes
        %construct the list of all represented genes
        %and their rank in the different lists indexed on probe nb
        GeneNames={};
        %find PsNb (can be different from probestNb field which indicates the number of
        %probesets in networks
        PsNb=1;
        for ProbeNbL=1:ProbeNb
            if ~isempty(PsRanks{ProbeNbL})
                for PsL=1:length(PsRanks{ProbeNbL})
                    PsNb=max(PsNb,max(PsRanks{ProbeNbL}{PsL}));
                end
            end
        end
        if PsNb~=NetPsNb
            sprintf('PsNb in networks: %u , in PSAWNpy: %u',NetPsNb,PsNb)
        end
        %synthetize results by indexing on probesets

        PsInfo=[];
        for PsL=1:PsNb
            PsInfo.geneNames{PsL}=cell(ProbeNb,1);
            PsInfo.exonNames{PsL}=cell(ProbeNb,1);
            PsInfo.exonRanks{PsL}=cell(ProbeNb,1);
            PsInfo.exonProbeNbs{PsL}=cell(ProbeNb,1);
            PsInfo.lastExon{PsL}=cell(ProbeNb,1);
            PsInfo.lastGroup{PsL}=cell(ProbeNb,1);
            PsInfo.groupRanks{PsL}=cell(ProbeNb,1);
            PsInfo.groupProbeNbs{PsL}=cell(ProbeNb,1);
            PsInfo.transcripts{PsL}=cell(ProbeNb,1);
            PsInfo.notTranscripts{PsL}=cell(ProbeNb,1);
            PsInfo.inGeneProbeNbs{PsL}=cell(ProbeNb,1);
            PsInfo.notInExonProbeNbs{PsL}=cell(ProbeNb,1);
        end
            
        for PsL=1:PsNb
            for ProbeNbL=1:ProbeNb
                PsInfo.geneNames{PsL}{ProbeNbL}={};
                PsInfo.exonNames{PsL}{ProbeNbL}={};
                PsInfo.exonRanks{PsL}{ProbeNbL}={};
                PsInfo.exonProbeNbs{PsL}{ProbeNbL}={};
                PsInfo.lastExon{PsL}{ProbeNbL}=[];
                PsInfo.lastGroup{PsL}{ProbeNbL}=[];
                PsInfo.groupRanks{PsL}{ProbeNbL}={};
                PsInfo.groupProbeNbs{PsL}{ProbeNbL}={};
                PsInfo.transcripts{PsL}{ProbeNbL}={};
                PsInfo.notTranscripts{PsL}{ProbeNbL}={};
                PsInfo.inGeneProbeNbs{PsL}{ProbeNbL}=[];
                PsInfo.notInExonProbeNbs{PsL}{ProbeNbL}=[];
            end
        end

        %recover GeneNames
        for ProbeNbL=1:ProbeNb
            if ~isempty(GeneName{ProbeNbL})
                GeneNames=[GeneNames;GeneName{ProbeNbL}];
            end
        end
        GeneNames=unique(GeneNames);
        GeneNb=length(GeneNames);

        %gene nb x probe nb
        TargetedGenes=zeros(GeneNb,ProbeNb);
        TargetingPsRanks=cell(GeneNb,ProbeNb);
        for ProbeNbL=1:ProbeNb
            if ~isempty(GeneName{ProbeNbL})
                for GeneL=1:length(GeneName{ProbeNbL})
                    %find the position of the current gene in the unique list
                    GenePos=strmatch(GeneName{ProbeNbL}{GeneL},GeneNames,'exact');
                    TargetedGenes(GenePos,ProbeNbL)=GeneL;
                    TargetingPsRanks{GenePos,ProbeNbL}=PsRanks{ProbeNbL}{GeneL};
                    %process all the probeset targeting the current gene
                    for PsL=1:length(PsRanks{ProbeNbL}{GeneL})
                        CurrPsRank=PsRanks{ProbeNbL}{GeneL}(PsL);                        
                        PsInfo.geneNames{CurrPsRank}{ProbeNbL}{end+1}=GeneName{ProbeNbL}{GeneL};
                        try
                            PsInfo.exonNames{CurrPsRank}{ProbeNbL}{end+1}=ExonNames{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.exonNames{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end
                        try
                            PsInfo.exonRanks{CurrPsRank}{ProbeNbL}{end+1}=ExonRanks{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.exonRanks{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end

                        try
                            PsInfo.exonProbeNbs{CurrPsRank}{ProbeNbL}{end+1}=ExonProbeNbs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.exonProbeNbs{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end
                        try
                            PsInfo.lastExon{CurrPsRank}{ProbeNbL}(end+1)=LastExon{ProbeNbL}(GeneL);
                        catch
                            PsInfo.lastExon{CurrPsRank}{ProbeNbL}(end+1)=-1;
                        end
                        try
                            PsInfo.lastGroup{CurrPsRank}{ProbeNbL}(end+1)=LastGroup{ProbeNbL}(GeneL);
                        catch
                            PsInfo.lastGroup{CurrPsRank}{ProbeNbL}(end+1)=-1;
                        end
                        try
                            PsInfo.groupRanks{CurrPsRank}{ProbeNbL}{end+1}=GroupRanks{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.groupRanks{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end
                        try
                            PsInfo.groupProbeNbs{CurrPsRank}{ProbeNbL}{end+1}=GroupProbeNbs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.groupProbeNbs{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end
                        try
                            PsInfo.transcripts{CurrPsRank}{ProbeNbL}{end+1}=TargetedTs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.transcripts{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end
                        try
                            PsInfo.notTranscripts{CurrPsRank}{ProbeNbL}{end+1}=NotTargetedTs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo.notTranscripts{CurrPsRank}{ProbeNbL}{end+1}=[];
                        end
                        PsInfo.inGeneProbeNbs{CurrPsRank}{ProbeNbL}(end+1)=InGeneProbeNbs{ProbeNbL}{GeneL}(PsL);
                        PsInfo.notInExonProbeNbs{CurrPsRank}{ProbeNbL}(end+1)=NotInExonProbeNbs{ProbeNbL}{GeneL}(PsL);
                    end
                end
            end
        end

        %keep notInExonProbeNb only if the gene is not targeted
        for PsL=1:PsNb
            %recover GeneNames that are targeted at least in one exon
            TargetedGeneNames={};
            for ProbeNbL=2:ProbeNb
                TargetedGeneNames=[TargetedGeneNames;PsInfo.geneNames{PsL}{ProbeNbL}'];
            end
            TargetedGeneNames=unique(TargetedGeneNames);
            for ProbeNbL=2:ProbeNb
                if ~isempty(PsInfo.geneNames{PsL}{ProbeNbL})
                    for GeneL=1:length(PsInfo.geneNames{PsL}{ProbeNbL})
                        CurrGeneName=PsInfo.geneNames{PsL}{ProbeNbL}{GeneL};
                        Pos=strmatch(CurrGeneName,TargetedGeneNames,'exact');
                        if ~isempty(Pos)
                            Pos=strmatch(CurrGeneName,PsInfo.geneNames{PsL}{ProbeNbL},'exact');
                            PsInfo.notInExonProbeNbs{PsL}{ProbeNbL}(Pos)=0;
                        end
                    end
                end
            end
        end
       
        %save data
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)))
        if TypeL==1
            EGeneName=GeneName;
            EGeneNames=GeneNames;
            ETargetedGenes=TargetedGenes;
            ETargetingPsRanks=TargetingPsRanks;
            EPsNames=PsNames;
            EExonRanks=ExonRanks;
            EPsRanks=PsRanks;
            EExonNames=ExonNames;
            EExonProbeNbs=ExonProbeNbs;
            ELastExon=LastExon;
            ELastGroup=LastGroup;
            EGroupRanks=GroupRanks;
            EGroupProbeNbs=GroupProbeNbs;
            ETargetedTs=TargetedTs;
            ENotTargetedTs=NotTargetedTs;
            EPsInfo=PsInfo;
            eval(sprintf('save m%u_probesets_by_%s_gene EGeneName EGeneNames ETargetedGenes ETargetingPsRanks EPsNames EPsRanks EExonNames EExonRanks EExonProbeNbs ELastExon ELastGroup EGroupRanks EGroupProbeNbs ETargetedTs ENotTargetedTs',ChipRank,Type))
            eval(sprintf('save m%u_%s_psinfo EPsInfo',ChipRank,Type))
            clear GeneName GeneNames TargetedGenes TargetingPsRanks PsNames PsRanks ExonNames ExonRanks ExonProbeNbs LastExon LastGroup GroupRanks GroupProbeNbs TargetedTs NotTargetedTs PsInfo TargetedT NotTargetedT           
            clear EGeneName EGeneNames ETargetedGenes ETargetingPsRanks EPsNames EPsRanks EExonNames EExonRanks EExonProbeNbs ELastExon ELastGroup EGroupRanks EGroupProbeNbs ETargetedTs ENotTargetedTs EPsInfo

        else
            AGeneName=GeneName;
            AGeneNames=GeneNames;
            ATargetedGenes=TargetedGenes;
            ATargetingPsRanks=TargetingPsRanks;
            APsNames=PsNames;
            APsRanks=PsRanks;
            AExonNames=ExonNames;
            AExonRanks=ExonRanks;
            AExonProbeNbs=ExonProbeNbs;
            ALastExon=LastExon;            
            ALastGroup=LastGroup;            
            AGroupRanks=GroupRanks;
            AGroupProbeNbs=GroupProbeNbs;
            ATargetedTs=TargetedTs;
            ANotTargetedTs=NotTargetedTs;
            APsInfo=PsInfo;
            eval(sprintf('save m%u_probesets_by_%s_gene AGeneName AGeneNames ATargetedGenes ATargetingPsRanks APsNames APsRanks AExonNames AExonRanks AExonProbeNbs ALastExon ALastGroup AGroupRanks AGroupProbeNbs ATargetedTs ANotTargetedTs',ChipRank,Type))
            eval(sprintf('save m%u_%s_psinfo APsInfo',ChipRank,Type))
            clear GeneName GeneNames TargetedGenes TargetingPsRanks PsNames PsRanks ExonNames ExonRanks ExonProbeNbs LastExon LastGroup GroupRanks GroupProbeNbs TargetedTs NotTargetedTs PsInfo TargetedT NotTargetedT
            clear AGeneName AGeneNames ATargetedGenes ATargetingPsRanks APsNames APsRanks AExonNames AExonRanks AExonProbeNbs ALastExon ALastGroup AGroupRanks AGroupProbeNbs ATargetedTs ANotTargetedTs APsInfo
        end
    end
end



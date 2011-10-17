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
% EPsInfo{PsRank}{ProbeNb}.exonProbeNbs
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
%     EExonProbeNbs: cell(n+1,1) containing the list of probeset names of each input file
%       EGroupRanks: cell(n+1,1) containing the list of grouped exons
%    EGroupProbeNbs: cell(n+1,1) containing the list of number of probes targeting grouped
%                    exons
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
                [GeneName{ProbeNbL},PsNames{ProbeNbL},PsRanks{ProbeNbL},...
                    ExonName{ProbeNbL},ExonProbeNbs{ProbeNbL},GroupRanks{ProbeNbL},GroupProbeNbs{ProbeNbL},...
                    TargetedT{ProbeNbL},NotTargetedT{ProbeNbL},NotInExonProbeNb{ProbeNbL},InGeneProbeNb{ProbeNbL}]=textread(sprintf('%s_m%u_probesets_by_gene_%02u.txt',Type,ChipRank,ProbeNbL-1),'%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t','bufsize',40000);
                ProbeNb=min(PROBE_NB,ProbeNbL);
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
        % which transcform string into cell
        for ProbeNbL=1:ProbeNb
            if ~isempty(ExonProbeNbs{ProbeNbL})
                for GeneL=1:length(ExonProbeNbs{ProbeNbL})
                    ExonProbeNbs{ProbeNbL}{GeneL}=eval(ExonProbeNbs{ProbeNbL}{GeneL});
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

        ExonNames=cell(size(ExonName));
        for ProbeNbL=1:ProbeNb
            if ~isempty(ExonName{ProbeNbL})
                ExonNames{ProbeNbL}=cell(size(ExonName{ProbeNbL}));
                for GeneL=1:length(ExonName{ProbeNbL})
                    ExonNames{ProbeNbL}{GeneL}=eval(ExonName{ProbeNbL}{GeneL});
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
        PsInfo=cell(PsNb,1);
        for PsL=1:PsNb
            for ProbeNbL=1:ProbeNb
                PsInfo{PsL}{ProbeNbL}.geneNames={};
                PsInfo{PsL}{ProbeNbL}.exonNames={};
                PsInfo{PsL}{ProbeNbL}.exonProbeNbs={};
                PsInfo{PsL}{ProbeNbL}.groupRanks={};
                PsInfo{PsL}{ProbeNbL}.groupProbeNbs={};
                PsInfo{PsL}{ProbeNbL}.transcripts={};
                PsInfo{PsL}{ProbeNbL}.notTranscripts={};
                PsInfo{PsL}{ProbeNbL}.inGeneProbeNbs=[];
                PsInfo{PsL}{ProbeNbL}.notInExonProbeNbs=[];
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
                        PsInfo{CurrPsRank}{ProbeNbL}.geneNames{end+1}=GeneName{ProbeNbL}{GeneL};
                        try
                            PsInfo{CurrPsRank}{ProbeNbL}.exonNames{end+1}=ExonNames{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo{CurrPsRank}{ProbeNbL}.exonNames{end+1}=[];
                        end
                        try
                            PsInfo{CurrPsRank}{ProbeNbL}.exonProbeNbs{end+1}=ExonProbeNbs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo{CurrPsRank}{ProbeNbL}.exonProbeNbs{end+1}=[];
                        end
                        try
                            PsInfo{CurrPsRank}{ProbeNbL}.groupRanks{end+1}=GroupRanks{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo{CurrPsRank}{ProbeNbL}.groupRanks{end+1}=[];
                        end
                        try
                            PsInfo{CurrPsRank}{ProbeNbL}.groupProbeNbs{end+1}=GroupProbeNbs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo{CurrPsRank}{ProbeNbL}.groupProbeNbs{end+1}=[];
                        end
                        try
                            PsInfo{CurrPsRank}{ProbeNbL}.transcripts{end+1}=TargetedTs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo{CurrPsRank}{ProbeNbL}.transcripts{end+1}=[];
                        end
                        try
                            PsInfo{CurrPsRank}{ProbeNbL}.notTranscripts{end+1}=NotTargetedTs{ProbeNbL}{GeneL}{PsL};
                        catch
                            PsInfo{CurrPsRank}{ProbeNbL}.notTranscripts{end+1}=[];
                        end
                        PsInfo{CurrPsRank}{ProbeNbL}.inGeneProbeNbs(end+1)=InGeneProbeNbs{ProbeNbL}{GeneL}(PsL);
                        PsInfo{CurrPsRank}{ProbeNbL}.notInExonProbeNbs(end+1)=NotInExonProbeNbs{ProbeNbL}{GeneL}(PsL);
                    end
                end
            end
        end

        %keep notInExonProbeNb only if the gene is not targeted
        for PsL=1:PsNb
            %recover GeneNames that are targeted at least in one exon
            TargetedGeneNames={};
            for ProbeNbL=2:ProbeNb
                TargetedGeneNames=[TargetedGeneNames;PsInfo{PsL}{ProbeNbL}.geneNames'];
            end
            TargetedGeneNames=unique(TargetedGeneNames);
            for ProbeNbL=2:ProbeNb
                if ~isempty(PsInfo{PsL}{ProbeNbL}.geneNames)
                    for GeneL=1:length(PsInfo{PsL}{ProbeNbL}.geneNames)
                        CurrGeneName=PsInfo{PsL}{ProbeNbL}.geneNames{GeneL};
                        Pos=strmatch(CurrGeneName,TargetedGeneNames,'exact');
                        if ~isempty(Pos)
                            Pos=strmatch(CurrGeneName,PsInfo{PsL}{ProbeNbL}.geneNames,'exact');
                            PsInfo{PsL}{ProbeNbL}.notInExonProbeNbs(Pos)=0;
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
            EPsRanks=PsRanks;
            EExonNames=ExonNames;
            EExonProbeNbs=ExonProbeNbs;
            EGroupRanks=GroupRanks;
            EGroupProbeNbs=GroupProbeNbs;
            ETargetedTs=TargetedTs;
            ENotTargetedTs=NotTargetedTs;
            EPsInfo=PsInfo;
            clear GeneName GeneNames TargetedGenes TargetingPsRanks PsNames PsRanks ExonNames ExonProbeNbs GroupRanks GroupProbeNbs TargetedTs NotTargetedTs PsInfo TargetedT NotTargetedT
            eval(sprintf('save m%u_probesets_by_%s_gene EGeneName EGeneNames ETargetedGenes ETargetingPsRanks EPsNames EPsRanks EExonNames EExonProbeNbs EGroupRanks EGroupProbeNbs ETargetedTs ENotTargetedTs EPsInfo',ChipRank,Type))
            clear EGeneName EGeneNames ETargetedGenes ETargetingPsRanks EPsNames EPsRanks EExonNames EExonProbeNbs EGroupRanks EGroupProbeNbs ETargetedTs ENotTargetedTs EPsInfo

        else
            AGeneName=GeneName;
            AGeneNames=GeneNames;
            ATargetedGenes=TargetedGenes;
            ATargetingPsRanks=TargetingPsRanks;
            APsNames=PsNames;
            APsRanks=PsRanks;
            AExonNames=ExonNames;
            AExonProbeNbs=ExonProbeNbs;
            AGroupRanks=GroupRanks;
            AGroupProbeNbs=GroupProbeNbs;
            ATargetedTs=TargetedTs;
            ANotTargetedTs=NotTargetedTs;
            APsInfo=PsInfo;
            clear GeneName GeneNames TargetedGenes TargetingPsRanks PsNames PsRanks ExonNames ExonProbeNbs GroupRanks GroupProbeNbs TargetedTs NotTargetedTs PsInfo TargetedT NotTargetedT
            eval(sprintf('save m%u_probesets_by_%s_gene AGeneName AGeneNames ATargetedGenes ATargetingPsRanks APsNames APsRanks AExonNames AExonProbeNbs AGroupRanks AGroupProbeNbs ATargetedTs ANotTargetedTs APsInfo',ChipRank,Type))
            clear AGeneName AGeneNames ATargetedGenes ATargetingPsRanks APsNames APsRanks AExonNames AExonProbeNbs AGroupRanks AGroupProbeNbs ATargetedTs ANotTargetedTs APsInfo
        end
    end
end



%=======================%
% FUNCTION MAKE_PSPAIRS %
%=======================%

%MAKE_PSPAIRS read info on probe sets and construct different kinds of pairs of probe sets
% targeting eventually the same gene(s) (duplicates)

%INPUT PARAMETERS

%INPUT
%
% 1     ProbeNbLimit: is the minimal number of probes that a probeset must have in a gene
% 2    TargetedGenes: contains matrix (nb of targeted genes x nb  of probes in the targeted
%                     gene) in cell (Ensembl & AceView)
%                     The matrix line rank corresponds to a gene of same rank in a list of
%                     genes (GeneName, not loaded here)
%                     The number contained in the matrix, if not null, is the rank of
%                     the gene in partial lists of genes (geneNames, not loaded here)
%                     example : if 56 is found at position (45,4) of the matrix,
%                     it means that the gene GeneName{45} is the 56th in GeneNames{4}
%                     If ProbeNbLimit>1, TargetedGenes is not used
% 3 TargetingPsRanks: cel with same dimensions than TargetedGenes
%                     Indicates the rank(s) of the probe set(s) that target the gene at
%                     a given probe nb. If ProbeNbLimit>1, TargetingPsRanks is not used
% 4         TestFlag: if =1 indicates that similarity is calculated to find limits on corr,
%                     anti and pv(node neighbourhodd similarity)); if =0 indicates that
%                     similarity is calculated to merge probe sets
% 5       SingleFlag: used if TestFlag==1 and indicates if one makes a difference between
%                     probe sets which target a single gene and those which target several
%                     genes (SingleFlag=1 for single genes or 0 for multiple genes)or if one
%                     does not make this difference (TestFlag=0 <=> SingleFlag=[]);
% 6           PsInfo: PsInfo{Type}{PsRank}{ProbeNb+1} gives the genes
%                     Type=1 => ENSEMBL, Type=2 => ACEVIEW)targetted by ProbeNb probes
%                     of the probeset of rank PsRank
%                     PsInfo{Type}{PsRank}{1} GIVES THE GENES THAT ARE TARGETED OUTSIDE EXONS
% 7        PsProbeNb: the number of probes in a normal probe set
%                     (few probe sets may have more probes than PsProbeNb).

%OUTPUT
% If TestFlag==0 no difference is made between single of multiple targets (SingleFalg is not
% used                
% If TestFlag==1 and SingleFlag=1 only probe sets targeting a single gene are considered
% If TestFlag==1 and SingleFlag=0 only probe sets targeting several genes are considered
%
% 1            DupStat : distribution of the number of targeted genes
% 2            DupRank : probe set ranks belonging to each partition on the number of
%                        targeted genes
% 3  GeneRankDuplicate : a list of gene ranks (referending GeneName list) that are repeated
%                        as much as they are targeted by different probe sets. If GeneName 45
%                        is targeted by 3 probe set, 45,45,45, is present in this list
%                        (Ensembl and AceView genes are respectively at the begining
%                        and at the end ot the list)
% 4    EnsDuplicateOut : list of couple of probe sets targeting the same Ensembl gene(s)
%                        outside of their exons
% 5    AceDuplicateOut : list of couple of probe sets targeting the same AceView gene(s)
%                        outside of their exons
% 6     EnsGeneNameOut : list of the Ensembl gene names targeted outside of their exons
% 7     AceGeneNameOut : list of the AceView gene names targeted outside of their exons
% 8          Duplicate : list of couple of probe sets targeting the same gene(s)
%                        in their exons
% 9       DuplicateOut : list of couple of probe sets targeting the same gene(s)
%                        out of their exons
% If ProbeNbLimit==1:
% 11      DuplicateLow : 10000 couples of randomly matched probe sets present in Duplicate
%                        and targeting with less than 3 probes
% 10  DuplicateLowHigh : 10000 couples of randomly matched probe sets,
%                        one present in DuplicateHigh and the other in DuplicateLow
% 12     DuplicateHigh : 10000 couples of randomly matched probe sets present in Duplicate
%                        and targeting with more or = than max(1,ProbeNbLimit-2 probes)


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

function [DupStat,DupRank,GeneRankDuplicate,EnsDuplicateOut,AceDuplicateOut,EnsGeneNameOut,AceGeneNameOut,Duplicate,DuplicateOut,DuplicateLowHigh,...
    DuplicateLow,DuplicateHigh]=make_pspairs(ProbeNbLimit, TargetedGenes, TargetingPsRanks,TestFlag,SingleFlag, PsInfo, PsProbeNb,AceFlag)


if AceFlag
    TypeNb=2;
else
    TypeNb=1;
end
EnsDuplicateOut=[];
AceDuplicateOut=[];
DuplicateLowHigh={};
DuplicateLow={};
DuplicateHigh={};
GeneRankDuplicate=[];
EnsGeneNameOut=[];
AceGeneNameOut=[];


%construct Duplicate
%Read input
PsNb=length(PsInfo{1});
%DupByBlocOut=cell(2,1);
PsRankHigh=[];


if TestFlag==0
    %recover all probeset targetting a gene with at least ProbeNbLimit probes
    %eliminate non used information (probeset with less than ProbeNbLimit
    %probes targetting genes
    %PROBE NB IS SHIFTED BY +1 IN PsInfo, TargettingPsRanks & TargetedGenes
    %=> (:,1:ProbeNbLimit)=[] clear Probe sets with <ProbeNbLimit
    %probes targeting a gene
    %the first column (i.e. probe nb=0)is always deleted
    for TypeL=1:TypeNb        
            TargetingPsRanks{TypeL}(:,1:ProbeNbLimit)=[];
        TargetedGenes{TypeL}(:,1:ProbeNbLimit)=[];
    end
    %cell array of list of probe sets targetting the same gene (either ensembl of aceview)
    Duplicate={};
    DuplicateOut=cell(2,1);
    GeneRankDuplicate=[];
    for TypeL=1:TypeNb
        %process each gene
        for GeneL=1:length(TargetedGenes{TypeL})
            %position of cell containing the rank(s) of probe set(s)
            %that target the current gene with at least ProbeNbLimit
            %probes
            PsPos=find(TargetedGenes{TypeL}(GeneL,:)~=0);
            if ~isempty(PsPos)
                %recover all the indexes of the probe sets targeting
                %the current gene with at least ProbeNbLimit probes
                CurrDuplicate=[];
                %construct a list of N identical gene rank (GeneL) that correspond to the N probe sets registered in CurrDuplicate
                CurrGeneRank=[];
                for PosL=1:length(PsPos)
                    CurrDuplicate=[CurrDuplicate,TargetingPsRanks{TypeL}{GeneL,PsPos(PosL)}];
                    CurrGeneRank=[CurrGeneRank;ones(length(TargetingPsRanks{TypeL}{GeneL,PsPos(PosL)}),1)*GeneL];
                    if TypeL==1
                        %assignation to high for random pairing
                        PsRankHigh=[PsRankHigh;TargetingPsRanks{TypeL}{GeneL,PsPos(PosL)}'];
                    end
                end
                if length(CurrDuplicate)>1
                    Duplicate{end+1,1}=sort(CurrDuplicate);
                    GeneRankDuplicate=[GeneRankDuplicate;CurrGeneRank];
                end

            end
        end

        


        %find duplicates for Ps targetting up, down or intron with at lest probeNbLimit probes
        % and not exons with more than probeNbLimit probes
        % (COULD TARGET EXONS WITH LESS)

        %find Ps without gene targeted in their exons with at least ProbNbLimit probes
        NotTarget=ones(PsNb,1);
        for DupL=1:length(Duplicate)
            NotTarget(Duplicate{DupL})=0;
        end
        %search among these genes those which have more than ProbeNbLimit probes
        %inside gene but not in exons these genes are in position one of PsInfo that is probenb=0 in
        %exons
        NotTarget=find(NotTarget);
        if ~isempty(NotTarget)
            %construct matched list of gene names  and  ps ranks
            GeneName={};
            PsRank=[];
            for PsL=1:length(NotTarget)
                CurrPsRank=NotTarget(PsL);
                Pos=find(PsInfo{TypeL}{CurrPsRank}{1}.notInExonProbeNbs>=ProbeNbLimit);
                if ~isempty(Pos)
                    for PosL=1:length(Pos)
                        PsRank=[PsRank;CurrPsRank];
                        GeneName=[GeneName;PsInfo{TypeL}{CurrPsRank}{1}.geneNames{Pos(PosL)}];
                    end
                end
            end
            %search genes targeted by at least two probe sets and group
            %these probe sets
            if ~isempty(GeneName)
                [GeneName,SortIndex]=sort(GeneName);
                PsRank=PsRank(SortIndex);
                MemGeneName=GeneName{1};
                CurrDuplicate=PsRank(1);

                GeneNames={};
                for GeneL=2:length(GeneName)
                    if isequal(GeneName{GeneL},MemGeneName)
                        CurrDuplicate=[CurrDuplicate,PsRank(GeneL)];
                    else
                        if length(CurrDuplicate)>1
                            DuplicateOut{TypeL}{end+1,1}=sort(CurrDuplicate);
                            GeneNames{end+1}=GeneName{GeneL};
                        end
                        MemGeneName=GeneName{GeneL};
                        CurrDuplicate=PsRank(GeneL);
                    end
                end
                %process the last gene
                if length(CurrDuplicate)>1
                    DuplicateOut{TypeL}{end+1,1}=sort(CurrDuplicate);
                    GeneNames{end+1}=GeneName{GeneL};
                end

                if TypeL==1
                    EnsDuplicateOut=DuplicateOut{1};
                    EnsGeneNameOut=GeneNames;
                else
                    AceDuplicateOut=DuplicateOut{2};
                    AceGeneNameOut=GeneNames;
                end
            else
                if TypeL==1
                    EnsDuplicateOut={};
                    EnsGeneNameOut={};
                else
                    AceDuplicateOut={};
                    AceGeneNameOut={};
                end
            end
        else
            if TypeL==1
                EnsDuplicateOut={};
                EnsGeneNameOut={};
            else
                AceDuplicateOut={};
                AceGeneNameOut={};
            end
        end
    end

    %construct random couples of ps targeting each an a priori different
    %gene in its exon with more or equal ProbeNbLimit probes
    PsRank=unique(PsRankHigh);
    PsRank=PsRank(randperm(length(PsRank)));
    if length(PsRank)>20000
        Grp1=PsRank(1:10000);
        Grp2=PsRank(10001:20000);
    else
        Limit=floor(length(PsRank)/2);
        GrpNb=ceil(10000/Limit);
        Grp1=PsRank(1:Limit);
        Grp2=PsRank(Limit:end);
        MGrp1=Grp1;
        MGrp2=Grp2;
        for GrpL=1:GrpNb
            Grp1=[Grp1;MGrp1(randperm(Limit))];
            Grp2=[Grp2;MGrp2(randperm(Limit))];
        end
    end
    for RankL=1:10000
        RefPsRank=Grp1(RankL);
        TestPsRank=Grp2(RankL);
        DuplicateHigh{end+1,1}=[RefPsRank,TestPsRank];
    end
else
    %recover Ps targetting exactly one gene with at least ProbeNbLimit probes in exons (SingleFlag=1)
    %or Ps targeting more than one gene  with at least ProbeNbLimit probes in exons (SingleFlag=0)
    %construct three control set (only on Ensembl genes and POG):
    %CONTROL SET1:
    %one Ps targetting exactly one gene but with 1,2 or 3 probes
    %CONTROL SET2:
    %one Ps targetting exactly one gene with max to max-2 probe nb (e.g. : 16,15 or 14 probes)
    %from control set1 & 2 are constructed random pairs of probe sets,
    %two belonging to SET1
    %two belonging to SET2
    %one belonging to SET1, the other to SET2.
    %CONTROL SET3:
    %Ps targetting only one gene with at least ProbeNbLimit probes inside the up, down or intron (nothing in exons or splice)


    %CONSTRUCT DUPLICATE AND CONTROL SET1, SET2 & SET3
    %ps targetting exactly one gene  with at least ProbeNbLimit probes
    Duplicate={};
    GeneNameHigh={};
    GeneNameLow={};
    PsRankHigh=[];
    PsRankLow=[];
    GeneNameOut={};
    PsRankOut=[];
    if SingleFlag
        PsRank=[];
        GeneName={};
        ProbeBinNb=length(PsInfo{1}{1});
        %process each probeset
        for PsL=1:length(PsInfo{1})
            GeneNb=0;
            %process each ProbeNb
            %ProbeNb in PsInfo is shifted +1
            %detect if the current probe set target only one gene in
            %exons => start at position 2
            for ProbeNbL=2:ProbeBinNb
                if ~isempty(PsInfo{1}{PsL}{ProbeNbL}.geneNames)
                    if length(PsInfo{1}{PsL}{ProbeNbL}.geneNames)==1
                        GeneNb=GeneNb+1;
                        CurrGeneName=PsInfo{1}{PsL}{ProbeNbL}.geneNames{1};
                        CurrProbeNb=ProbeNbL-1;
                        if GeneNb>1
                            break
                        end
                    else
                        GeneNb=2;
                        break
                    end
                end
            end
            %keep only probeset targetting exactly one gene
            if GeneNb==1
                if CurrProbeNb>=ProbeNbLimit
                    PsRank=[PsRank;PsL];
                    GeneName=[GeneName;CurrGeneName];
                end
                %assignation to low or high according to the probe nb (made only for Ensembl genes)
                if CurrProbeNb<=3
                    %CONTROL SET1
                    GeneNameLow=[GeneNameLow;CurrGeneName];
                    PsRankLow=[PsRankLow;PsL];
                elseif CurrProbeNb>=PsProbeNb-2
                    %CONTROL SET2
                    GeneNameHigh=[GeneNameHigh;CurrGeneName];
                    PsRankHigh=[PsRankHigh;PsL];
                end
            end

            %CONTROL SET3
            %probesets targetting only one gene in up, down or intron (only Ensembl genes) with equal or more than ProbeNbLimit
            if GeneNb==0
                if length(PsInfo{1}{PsL}{1}.geneNames)==1 & PsInfo{1}{PsL}{1}.notInExonProbeNbs(1)>=ProbeNbLimit
                    GeneNameOut=[GeneNameOut;PsInfo{1}{PsL}{1}.geneNames{1}];
                    PsRankOut=[PsRankOut;PsL];
                end
            end
        end % of PsL
    else
        %search probeset with multiple targets

        PsRank=[];
        GeneName={};
        ProbeBinNb=length(PsInfo{1}{1});

        %process each probeset
        for PsL=1:length(PsInfo{1})
            GeneNb=0;
            CurrGeneNames={};
            CurrProbeNbs=[];
            %process each ProbeNb
            %ProbeNb in PsInfo is shifted +1
            %detect if the current probe set target several genes in
            %exons => start at position ProbeNb
            for ProbeNbL=ProbeNbLimit+1:ProbeBinNb
                if ~isempty(PsInfo{1}{PsL}{ProbeNbL}.geneNames)
                    for GeneL=1:length(PsInfo{1}{PsL}{ProbeNbL}.geneNames)
                        GeneNb=GeneNb+1;
                        CurrGeneNames{end+1}=PsInfo{1}{PsL}{ProbeNbL}.geneNames{GeneL};
                        CurrProbeNbs=[CurrProbeNbs,ProbeNbL-1];
                    end
                end
            end
            %keep only probeset targetting several genes
            if GeneNb>1
                for GeneL=1:length(CurrGeneNames)
                    PsRank=[PsRank;PsL];
                    GeneName=[GeneName;CurrGeneNames{GeneL}];
                    %assignation to low or high according to the probe nb (made only for Ensembl genes)
                    if CurrProbeNbs(GeneL)<=3
                        %CONTROL SET1
                        GeneNameLow=[GeneNameLow;CurrGeneNames{GeneL}];
                        PsRankLow=[PsRankLow;PsL];
                    elseif CurrProbeNbs(GeneL)>=PsProbeNb-2
                        %CONTROL SET2
                        GeneNameHigh=[GeneNameHigh;CurrGeneNames{GeneL}];
                        PsRankHigh=[PsRankHigh;PsL];
                    end
                end
            end

            %CONTROL SET3
            %probesets targetting more tan one gene in up, down or intron (only Ensembl genes) with equal or more than ProbeNbLimit
            if GeneNb==0
                if length(PsInfo{1}{PsL}{1}.geneNames)>1
                    KeepIt=0;
                    KeepNb=0;
                    for GeneL=1:length(PsInfo{1}{PsL}{1}.geneNames)
                        if PsInfo{1}{PsL}{1}.notInExonProbeNbs(GeneL)>=ProbeNbLimit
                            KeepNb=KeepNb+1;
                            if KeepNb>1
                                KeepIt=1;
                                break;
                            end
                        end
                    end
                    %recover Gene iff there is at least two that are targeted by more than ProbeNbLimit
                    if KeepIt
                        for GeneL=1:length(PsInfo{1}{PsL}{1}.geneNames)
                            if PsInfo{1}{PsL}{1}.notInExonProbeNbs(GeneL)>=ProbeNbLimit
                                GeneNameOut=[GeneNameOut;PsInfo{1}{PsL}{1}.geneNames{GeneL}];
                                PsRankOut=[PsRankOut;PsL];
                            end
                        end
                    end
                end
            end
        end % of PsL
    end


    %PROCESS DUPLICATE
    %find duplicates for Ps targetting one (SingleFlag=1) or several (SingleFlag=0) gene(s) with at least ProbeNbLimit probes
    [GeneName,SortIndex]=sort(GeneName);
    PsRank=PsRank(SortIndex);   
    MemGeneName=GeneName{1};    
    CurrDuplicate=PsRank(1);
    %cell array of list of probesets targetting the same gene (either ensembl of aceview)
    for GeneL=2:length(GeneName)
        if isequal(GeneName{GeneL},MemGeneName)
            CurrDuplicate=[CurrDuplicate,PsRank(GeneL)];
        else
            CurrDuplicate=unique(CurrDuplicate);
            if length(CurrDuplicate)>1
                Duplicate{end+1,1}=sort(CurrDuplicate);
            end
            MemGeneName=GeneName{GeneL};
            CurrDuplicate=PsRank(GeneL);
        end
    end
    CurrDuplicate=unique(CurrDuplicate);
    if length(CurrDuplicate)>1
        Duplicate{end+1,1}=sort(CurrDuplicate);
    end



    %PROCESS CONTROL SET1, SET2 & SET3
    %find duplicates for Ps targetting one (SingleFlag=1) or several (SingleFlag=0) genes with at least ProbeNbLimit probes in up, down or intron

    %CONTROL SET1 & SET2
    DuplicateLowHigh={};
    if length(GeneNameLow)>100
        HighRank=randperm(length(GeneNameHigh));
        LowRank=randperm(length(GeneNameLow));
        RefNb=1;
        PsPairs=[];
        %try to construct ~ 10000 pairs of probe sets if possible
        while RefNb<min(101,length(GeneNameHigh))
            RefGeneName=GeneNameHigh{HighRank(RefNb)};
            RefPsRank=PsRankHigh(HighRank(RefNb));
            TestNb=1;
            while TestNb<min(101,length(GeneNameLow))
                if ~isequal(GeneNameLow{LowRank(TestNb)},RefGeneName)
                    TestPsRank=PsRankLow(LowRank(TestNb));
                    PsPairs=[PsPairs;sort([RefPsRank,TestPsRank])];
                end
                    TestNb=TestNb+1;
            end
            RefNb=RefNb+1;
        end
        PsPairs=unique(PsPairs,'rows');
        [Temp,SortIndex]=sort(PsPairs(:,1));
        PsPairs=PsPairs(SortIndex,:);
        for PairL=1:size(PsPairs,1)
            DuplicateLowHigh{end+1,1}=PsPairs(PairL,:);
        end
    end

    %CONTROL SET1
    DuplicateLow={};
    if length(GeneNameLow)>100
        PsRank=unique(PsRankLow);
        PsRank=PsRank(randperm(length(PsRank)));
        PsPairs=[];
        if length(PsRank)>20000
            Grp1=PsRank(1:1000);
            Grp2=PsRank(10001:20000);
        else
            Limit=floor(length(PsRank)/2);
            GrpNb=ceil(10000/Limit);
            Grp1=PsRank(1:Limit);
            Grp2=PsRank(Limit:end);
            MGrp1=Grp1;
            MGrp2=Grp2;
            for GrpL=1:GrpNb
                Grp1=[Grp1;MGrp1(randperm(Limit))];
                Grp2=[Grp2;MGrp2(randperm(Limit))];
            end
        end
        for RankL=1:10000
            RefPsRank=Grp1(RankL);
            TestPsRank=Grp2(RankL);
            PsPairs=[PsPairs;sort([RefPsRank,TestPsRank])];
        end
        PsPairs=unique(PsPairs,'rows');
        [Temp,SortIndex]=sort(PsPairs(:,1));
        PsPairs=PsPairs(SortIndex,:);
        for PairL=1:size(PsPairs,1)
            DuplicateLow{end+1,1}=PsPairs(PairL,:);
        end
    end


    %CONTROL SET2
    DuplicateHigh={};
    PsRank=unique(PsRankHigh);
    PsRank=PsRank(randperm(length(PsRank)));
    PsPairs=[];
    if length(PsRank)>20000
        Grp1=PsRank(1:10000);
        Grp2=PsRank(10001:20000);
    else
        Limit=floor(length(PsRank)/2);
        GrpNb=ceil(10000/Limit);
        Grp1=PsRank(1:Limit);
        Grp2=PsRank(Limit:end);
        MGrp1=Grp1;
        MGrp2=Grp2;
        for GrpL=1:GrpNb
            Grp1=[Grp1;MGrp1(randperm(Limit))];
            Grp2=[Grp2;MGrp2(randperm(Limit))];
        end
    end
    for RankL=1:10000
        RefPsRank=Grp1(RankL);
        TestPsRank=Grp2(RankL);
        PsPairs=[PsPairs;sort([RefPsRank,TestPsRank])];
    end
    PsPairs=unique(PsPairs,'rows');
    [Temp,SortIndex]=sort(PsPairs(:,1));
    PsPairs=PsPairs(SortIndex,:);
    for PairL=1:size(PsPairs,1)
        DuplicateHigh{end+1,1}=PsPairs(PairL,:);
    end


    %CONTROL SET3
    DuplicateOut={};
    if ~isempty(GeneNameOut)
        [GeneNameOut,SortIndex]=sort(GeneNameOut);
        PsRankOut=PsRankOut(SortIndex);
        MemGeneNameOut=GeneNameOut{1};
        CurrDuplicate=PsRankOut(1);        
        GeneNames={};
        for GeneL=2:length(GeneNameOut)
            if isequal(GeneNameOut{GeneL},MemGeneNameOut)
                CurrDuplicate=[CurrDuplicate,PsRankOut(GeneL)];
            else
                CurrDuplicate=unique(CurrDuplicate);
                if length(CurrDuplicate)>1
                    DuplicateOut{end+1,1}=sort(CurrDuplicate);
                    GeneNames{end+1}=GeneNameOut{GeneL};
                end
                MemGeneNameOut=GeneNameOut{GeneL};
                CurrDuplicate=PsRankOut(GeneL);
            end
        end
        CurrDuplicate=unique(CurrDuplicate);
        if length(CurrDuplicate)>1
            DuplicateOut{end+1,1}=sort(CurrDuplicate);
            GeneNames{end+1}=GeneNameOut{GeneL};
        end
        EnsDuplicateOut=DuplicateOut;
        EnsGeneNameOut=GeneNames;
    end
end

% OPTIM_CHIPDUP
%find the maximal number of duplicates
GeneNb=length(Duplicate);
MaxDup=0;
for GeneL=1:GeneNb
    MaxDup=max(MaxDup,length(Duplicate{GeneL}));
end

%recover statistics on the number of genes targeted by x duplicated probesets
DupStat=zeros(MaxDup-1,2);
DupStat(:,1)=[2:MaxDup]';
%recover the probeset ranks wich belongs to the group of x replicated
%probesets
DupRank=cell(MaxDup-1,1);
for GeneL=1:GeneNb
    DupNb=length(Duplicate{GeneL});
    DupStat(DupNb-1,2)=DupStat(DupNb-1,2)+1;
    DupRank{DupNb-1}=[DupRank{DupNb-1},Duplicate{GeneL}];
end
ZIndex=find(DupStat(:,2)==0);
DupStat(ZIndex,:)=[];
DupRank(ZIndex)=[];

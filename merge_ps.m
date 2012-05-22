%===================%
% FUNCTION MERGE_PS %
%===================%
%
% MERGE_PS finds the group of probe sets which target common transcripts on the basis
% of their similar behavior in several networks : high positive correlation, low
% negative correlation and a highly similar neighbourhood, as measured by pv(overlap) which
% is the p-value of observing a given number of common neighbors under a hypergeometric
% distribution
%
%
%
%INPUT PARAMETERS
%  1        Species: species
%  2       ChipRank: chip rank
%  3       NetRanks: ranks of networks used
%  4   ProbeNbLimit: minimum number of probes targeting a gene
%  5     PvCorrRank: pv(overlap) is calculated for corr limit >[0,40,50,60]. PvCorrRank
%                    indicates the corr limit to be used, by giving its index in
%                    the corr list([0,40,50,60])
%  6      StepRanks: list of merge_ps steps to be processed
%  7 NetFrequencies: list of net frequencies that mus be considered ([1,25,50,75,100],
%                    recommended)
%  8    DisplayFlag: indicates if figures must be displayed
%  9        SumFlag: indicates if positive networks are those where probe set pair 
%                    have positive correlation and are similar (=1) or only those where 
%                    probe sets are similar i.e. satisfy all the tested conditions 
%                    (positive and negative correlation and p-value of overlaping of
%                    their neighbourhood) (=0, recommended)
% 10        ValFlag: indicates if all values of corr, anti and pv of each pair are used 
%                    (=1, recommended), or only a single derived value (mean - std) (=0) to
%                    calculate limits used for testing probe set pairs
% 11       MeanFlag: indicates if limits are calculated from mean of four values (single and
%                    multiple targeted genes, and InSim and OutSim) (=1) or only from single
%                    targeted genes and InSim (=0, recommended)
% 12        AceFlag: indicates if AceView data are available
% 13       IdemFlag: indicates if probe set order is identical in file used by PsawnPy et in
%                    networks (=1 in general case, if =0 a m%u_net.txt file must exist 
%                    in K.dir.rawdata)
%
% MERGE_PS STEPS (StepRanks parameter)
% 1 CONSTRUCT NEWPS
% 2 CALCULATE LIMITS
% 3 COMPLETE NEWPS
% 4 CONSTRUCT PSBY
% 5 STAT ON THE THIRD EDGE IN PROBE SET TRIANGLES
% 6 SUMMARIZE (CONSTRUCT PSMATRIX)
% 7 FIGURES 51 TO 53
% 8 WRITE TXT FILES:
%
%   Write PsMatrix files:
%
%	 1: first probeset ID
% 	 2: gene ID assigned to the probesets (either Ensembl ID or eventually AceView or GOP id)
%	 3: gene name
%	 4: rank of the assigned gene to the current probeset
%	 5: position of the assigned gene in NewPs.geneNames
%	 6: target type of assigned gene (1 for probeset in exon, 0 otherwise)
%	 7: source type of assigned gene (1 for Ensembl, 2 for AceView, 3 for GOP)
%	 8: nb of probe targeting the assigned gene
%	 9: nb of not assigned genes targeted with the same nb of probes
%	10: nb of not assigned genes targeted with less nb of probes
%	11: nb of groups of transcripts corresponding to the assigned gene
%	12: rank of the parent probeset
%	13: Rank of the group of transcripts targeted by the current
%	    probeset in the assigned gene
%	14: Rank of the group(s) of transcripts targeted by the current
%	    probeset, if it is a pivot
%	15: [0,1] indicates if the current probeset is a pivot
%	16: [0,1] indicates if the current probeset is paired with a pivot
%	17: nb of probesets that do not target the assigned gene but target
%	    a common gene with the current probeset        
%	18: nb of genes that are targeted by probesets that do not target the
%	    assigned gene
%	19: nb of genes that are targeted by other probeset with a nb
%	    of probes higher than the number of probes of the current probeset that
%	    target the assigned gene
%	20: ClassRank
%	21: and beyond: nb of  other genes targeted with all possible nb of probes (from 1 to ProbeNb)
%
%   Write a file describing all the existing probe set pairs
%    1: gene ID assigned to the two probesets forming a pair (either Ensembl ID or eventually AceView or GOP id)
%	 2: gene name
%	 3: first probeset ID
%	 4: second probeset ID
%	 5: first probeset rank in PsMatrix
%	 6: second probeset rank in PsMatrix
%	 7: indicates if the probesets are similar in 1% PsMatrix
%	 8: indicates if the probesets are similar in 25% PsMatrix
%	 9: indicates if the probesets are similar in 50% PsMatrix
%	10: indicates if the probesets are similar in 75% PsMatrix
%	11: indicates if the probesets are similar in 100% PsMatrix
%	     pivot information for fields 7 to 10, if the two probesets are similar:
%	     if first and second probesets are not a pivot => 1
%	     if only one of them is a pivot => 2
%	     if both are pivots => 3
%	12: probeset class (3=MS, 4=MM, 5=CX, 6=HX)
%	13: number of probes of the first probeset targeting the assigned gene
%	14: number of genes targeted by the first probeset with the same number of probes
%	15: number of genes targeted by the first probeset with an inferior number of probes
%	16: number of probes of the second probeset targeting the assigned gene
%	17: number of genes targeted by the second probeset with the same number of probes
%	18: number of genes targeted by the second probeset with an inferior number of probes
%	19: class rank of first probe set
%	20: class rank of snd probe set
%	21: v: the probeset pair is tested in all the networks,
%	    *: the probeset pair is absent of at leat one networks (~1% of all pairs)
%	22: total number of transcripts targeted by the two probesets
%	23: number of transcripts targeted in common by the two probesets
%	24: number of transcripts specifically targeted by the first probeset
%	25: number of transcripts specifically targeted by the second probeset
%	26: total number of exons targeted by the two probesets
%	27: number of probes located in exons targeted in common by the two probesets
%	28: number of probes located in exons specifically targeted by the first probeset 
%	29: number of probes located in exons specifically targeted by the second probeset 
%	30: total number of groups of exons targeted by the two probesets
%	31: number of probes located in groups of exons targeted in common by the two probesets
%	32: number of probes located in groups of exons specifically targeted by the first probeset 
%	33: number of probes located in groups of exons specifically targeted by the second probeset 
%	34: overlaping score for transcripts (column 23*100/column 22)
%	35: overlaping score for exons 
%		for each targeted exon a local weighted overlap score, using the number of probes of the
%		first (PNb1) and the second (PNb2) probeset targeting this exon, and the total nb
%		of probes targeting an exon (PNb)
%		(PNb1+PNb2)/PNb)*(min(PNb1,PNb2)/max(PNb1,PNb2)
%		The overlaping score is the mean of all local scores
%	36: overlaping score for groups of exons
%		same method used for the overlaping score for exons
%	37: percentage of first probeset probes located in the last exon
%	38: percentage of second probeset probes located in the last exon
%	39: indicates (1/0) if all probes of the first probeset are located in a single exon
%	40: indicates (1/0) if all probes of the second probeset are located in a single exon
%	41: percentage of first probeset probes located in the last group
%	42: percentage of second probeset probes located in the last group
%	43: indicates (1/0) if all probes of the first probeset are located in a single group
%	44: indicates (1/0) if all probes of the second probeset are located in a single group
%	45 to 67: iAceView information, idem to 22 to 44
%
%  Plots FIGURE 56
%
%SUB FUNCTIONS
%
%   LOADSIM
%     Load sim 
%     INPUT PARAMETERS
%      1   ChipRank : the rank of chip set model
%      2   NetRanks : a list of network rank
%      3 PvCorrRanks: five series of p-values are calculated (on edges with
%                     corr>=0,10,20,30,40,50)
%      4   FileName : the common part of or file name to be loaded
%
%     OUTPUT
%      AllVal : all CORR,ANTI and PV (p-values) of the networks
%         Sim : The last loaded Sim
%
%   CONSTRUCT_NEWPS
%      construct a new structure (NewPs) containing information about relationships
%      between probe sets
%      one record per probe set within this structure:
%         NewPs{PsL,1}.geneNames: EnsGeneID or Ace gene names targeted by
%                                   the probeset
%           NewPs{PsL,1}.probeNb: nb of probes in each targeted gene
%           NewPs{PsL,1}.psRanks: rank ot other probe set which target the
%                                 same genes
%            NewPs{PsL,1}.source: 1=Ensembl genes; 2=AceView genes (not found
%                                 in Ensembl)
%            NewPs{PsL,1}.target: 0= in GOP (group of probe sets that do not overlap a
%                                 1=in exon or splice; 2=in up, intron or down
%      constructs also
%           Genes.name: all the gene names targetted by at least one probe set;
%         Genes.source: source of the gene (1=Ensembl; 2=AceView)
%      and
%         PsBy.gene: for each position of Genes, indcates the ranks of the
%                    probe sets that target that gene
%
%   FILL_PAIRED
%      add information in NewPs
%         NewPs{PsL,1}.psRanks: other probeset targetting the sames genes
%            NewPs{PsL,1}.corr: mean corr between the current probeset and
%                               the other probe sets
%            NewPs{PsL,1}.anti: mean anti ...
%              NewPs{PsL,1}.pv: mean pv ...
%           NewPs{PsL,1}.repnb: number of network in which there are
%                               significative corr or anti values
%         NewPs{PsL,1}.stdcorr: std corr ...
%         NewPs{PsL,1}.stdanti: std anti ...
%           NewPs{PsL,1}.stdpv: std pv ...
%       new pairs of probe sets are recovered at this step which needs
%       to calculate their corr, anti, pv, repnb, stdcorr, stdanti & stdpv
%       in the networks used
%
%   CONSTRUCT_PSBY
%       for each targeted gene recover all the probe set that target this
%       gene and construct the matrix of their interactions (corr, anti,
%       pv, repnb, stdcorr, stdanti, stdpv)
%
%   STAT
%      OUTPUT

%        Stat, a structure containing information about genes targeted by probe sets
%           Stat.singleTargNb: The number of genes that are targetted only 
%                              by the current probe set
%           Stat.doubleTargNb: The number of genes that are targetted by the current probe
%                              set plus a single other one
%         Stat.multipleTargNb: The number of genes that are are targetted by the current 
%                              probe set plus two or more another probe sets
%          Stat.maxPsNb=zeros: The maximum nb of probe sets targetting a gene targeted 
%                              by the current probe set
%               Stat.grpSizes: the distribution of ps group sizes
%              Stat.linkTypes: the type of link (0: no corr, 1:don't pass the test,
%                              2 pass the test)
%               Stat.badLinks: Properties of pairs that don't pass the test but are in 
%                              a ps group:
%                              [PsL,GeneL,Node1,Node2,Nb of bad links,
%                              Nb of links,Nb of not significative corr,
%                              Nb of significative corr];
%                   Stat.hubs: During merging process of triangle, some transitory groups
%                              of probe sets are split into two new  groups. Single probe
%                              sets that are common to these two groups are taken away and
%                              considered as forming a hub which is in relation with 
%                              the two groups.

%        LinkedPs gives information on pairs of probe set that belong to the same group
%           LinkedPs{1}: they target a gene with a number of probe smaller than 
%                        the number of probes of the current ps targeting the assigned gene
%           LinkedPs{2}: they target a gene with a number of probe equal to the number 
%                        of probes of the current ps targeting the assigned gene
%                        [Ps1,Ps2,Target type of targeted gene,
%                        Nb of probes of Ps1 targeting the gene
%                        Nb of probes of Ps2 targeting the gene,
%                        Source type of the targeted gene,
%                        maximum nb of probes of Ps1 that target a gene,
%                        maximum nb of probes of Ps1 that target a gene];
%           LinkedPs{1}(PairPos,:): [Ps1,Ps2,NewPs{Ps1}.target(GenePos1),NewPs{Ps1}.
%                                    probeNb(GenePos1), NewPs{Ps2}.probeNb(GenePos2),
%                                    NewPs{Ps1}.source(GenePos1),max(NewPs{Ps1}.probeNb),
%                                    max(NewPs{Ps2}.probeNb)];
%
%   TRIANGLE_STAT
%      calculate statistics on probe set triangles
%
%   SUMMARIZE
%       fill PsMatrix:
%          1: rank of the assigned gene to the current probe set
%          2: position of the assigned gene in NewPs.geneNames
%          3: target type of assigned gene
%          4: source type of assigned gene
%          5: nb of probe targetting the assigned gene
%          6: nb of not assigned genes targetted with the same nb of probes
%          7: nb of not assigned genes targetted with less nb of probes
%          8: nb of groups of transcripts corresponding to the assigned gene
%          9: rank of the parent probe set
%         10: Rank of the group of transcripts targetted by the current
%             probe set in the assigned gene
%         11: Rank of the group(s) of transcripts targetted by the current
%             probe set, if it is a pivot
%         12: [0,1] indicates if the current probe set is a pivot
%         13: [0,1] indicates if the current probe set is paired with a pivot
%         14: nb of probe sets that do not target the assigned gene but target
%             a common gene with the current probe set        
%         15: nb of genes that are targetted by probe sets that do not target the
%             assigned gene
%         16: nb of genes that are targetted by other probe set with a nb
%             of probes higher than the number of probes of the current probe set that
%             target the assigned gene
%         17: ClassRank
%         18: and beyond: nb of  other genes targetted with all possible nb of probes
%
%   DISPLAY_PS
%      display FIG25a,FIG25b
%
%   DIST_PROBENB
%      calculate the percentage of different source and gene types
%


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

function merge_ps(Species,ChipRank,NetRanks,ProbeNbLimit,PvCorrRanks,StepRanks,NetFrequencies,DisplayFlag,SumFlag,ValFlag,MeanFlag,AceFlag,IdemFlag)
global K


%% MAIN
ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
ProbeNb=K.chip.probeNb(ChipPos);
NetNb=length(NetRanks);
DirName=fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank));
cd(DirName)
%LOAD NewPs IF EXISTS OR CONSTRUCT IT
FileName=sprintf('m%u_n%u_netnb%u_probenb%u_newps.mat',ChipRank,NetRanks(1),NetNb,ProbeNbLimit);

%% STEP1 CONSTRUCT NEWPS
if isempty(find(StepRanks==1))
    if exist(FileName,'file')
        eval(sprintf('load %s',FileName))
    else
        h=errordlg(sprintf('File %s does not exists => step 1 must be in StepRanks',FileName));
        waitfor(h)
        error('process canceled')
    end
else
    %nb of genes targetted by x probes (PsNb x ProbeNb matrix)
    eval(sprintf('load m%u_probesets_ensembl',ChipRank))
    if AceFlag
        eval(sprintf('load m%u_probesets_aceview',ChipRank))
        'DISPLAY_PS'
        DISPLAY_PS(Species,ChipRank,EnsExonGeneNbs,AceExonGeneNbs,ProbeNbLimit)
    else
        DISPLAY_PS(Species,ChipRank,EnsExonGeneNbs,[],ProbeNbLimit)
    end
    cd(DirName)
    if AceFlag
        %RECOVER ENSEMBL GENE NAMES CORRESPONDING TO ACEVIEW GENES
        cd(fullfile(K.dir.pydata,Species,'txt'))
        [AceGenesDict,EnsGDict]=textread(sprintf('%s_ens_by_ace_gene.txt',Species),'%s%s','delimiter','\t');
        EnsGenesDict=cell(length(AceGenesDict),1);
        AceGeneNb=0;
        %Ensembl Genes without AceView corresponding genes are at the end of the files
        for AceL=1:length(AceGenesDict)
            if ~isequal(AceGenesDict{AceL},'0')
                AceGeneNb=AceGeneNb+1;
                EnsGenesDict{AceL}=eval(EnsGDict{AceL});
            else
                break
            end
        end
        AceGenesDict=AceGenesDict(1:AceGeneNb);   
    end

    %recover Ps
    cd(DirName)
    CurrFileName=sprintf('m%u_compinfo_probenb%u.mat',ChipRank,ProbeNbLimit);
    load(CurrFileName)
    clear CompInfo

    'CONSTRUCT NEWPS'
    if AceFlag
        [NewPs,PsBy,Genes]=CONSTRUCT_NEWPS(Ps,EnsGenesDict,AceGenesDict);
    else
        [NewPs,PsBy,Genes]=CONSTRUCT_NEWPS(Ps,[],[]);
    end
    cd(DirName)
    eval(sprintf('save %s NewPs PsBy Genes',FileName))
end

%% STEP2 CALCULATE LIMITS
if ~isempty(find(StepRanks==2))    
    %RECOVER THE CURRENT CORR,ANTI &PV VALUES FOR CURRENT SIM
    cd(DirName)
    CurrFileName=sprintf('nodesim_probenb%u.mat',ProbeNbLimit);
    'LOADSIM'
    [AllVal,Sim]=LOADSIM(ChipRank,NetRanks,PvCorrRanks,CurrFileName);
    clear Temp    

    CurrFileName=sprintf('limits_m%u_n%u_netnb%u_rank%u_v%u_m%u.mat',ChipRank,NetRanks(1),length(NetRanks),PvCorrRanks,ValFlag,MeanFlag);
    if exist(CurrFileName,'file')
        'LOAD LIMITS'
        eval(sprintf('load %s',CurrFileName))
    else
        'CALCULATE LIMITS'
        Limit=calculate_limits(Species,ChipRank,ProbeNbLimit,[],NetRanks,PvCorrRanks,ValFlag,MeanFlag);        
    end        
    'FILL PAIRED (FIRST ROUND)'
    [NewPs,NotFound]=FILL_PAIRED(Limit,NewPs,Sim,AllVal,PsBy,Genes,1);
    cd(DirName)
    eval(sprintf('save %s NewPs Genes  PsBy',FileName))
    eval(sprintf('save %s NotFound',sprintf('m%u_dup_probenb%u_notfound.mat',ChipRank,ProbeNbLimit)));
end

%% STEP3 COMPLETE NEWPS
if ~isempty(find(StepRanks==3))     
    ProcessNet=zeros(NetNb,1);
    for NetL=1:NetNb
        CurrFileName=sprintf('m%u_n%u_nodesim_probenb%u_notfound.mat',ChipRank,NetRanks(NetL),ProbeNbLimit);
        if ~exist(CurrFileName,'file')
            ProcessNet(NetL)=1;
        end        
    end
    CurrNetRankList=NetRanks(find(ProcessNet));
    if ~isempty(CurrNetRankList)
        TestFlag=0;
        NotFoundFlag=1;
        'CALCULATE NODE SIM'
        calculate_nodesim(TestFlag,NotFoundFlag,ProbeNbLimit,ChipRank,CurrNetRankList,[0,40,50,60],AceFlag,IdemFlag)
    end

    %RECOVER THE CURRENT CORR,ANTI &PV VALUES FOR NOT FOUND
    %keep the last SIM loaded in LOADSIM to have PsRank information
    cd(DirName)
    CurrFileName=sprintf('nodesim_probenb%u_notfound.mat',ProbeNbLimit);
    'LOADSIM'
    [AllVal,NotFoundSim]=LOADSIM(ChipRank,NetRanks,PvCorrRanks,CurrFileName);

    
    CurrFileName=sprintf('limits_m%u_n%u_netnb%u_rank%u_v%u_m%u.mat',ChipRank,NetRanks(1),length(NetRanks),PvCorrRanks,ValFlag,MeanFlag);
    if exist(CurrFileName,'file')
        'LOAD LIMITS'
        eval(sprintf('load %s',CurrFileName))
    else        
        'CALCULATE LIMITS'    
        Limit=calculate_limits(Species,ChipRank,ProbeNbLimit,[],NetRanks,PvCorrRanks,ValFlag,MeanFlag);        
    end
   
    'FILL PAIRED (SECOND ROUND)'
    %FILL PV AND CORR VALUES (SECOND ROUND)
    [NewPs,Tmp]=FILL_PAIRED(Limit,NewPs,NotFoundSim,AllVal,PsBy,Genes,2);   
    cd(DirName)
    eval(sprintf('save %s NewPs Genes PsBy',FileName))
end

%% STEP4 CONSTRUCT PSBY
if ~isempty(find(StepRanks==4))    
    %process ps with targets only in out
    'CONSTRUCT_PSBY'
    %SYNTHESIS OF ALL RESULTS IN PS BY GENES
    PsBy=CONSTRUCT_PSBY(NewPs,PsBy,Genes);
    cd(DirName)
    eval(sprintf('save %s NewPs Genes PsBy',FileName))
end

%% STEP5 STAT ON THE THIRD EDGE IN PROBE SET TRIANGLE
if ~isempty(find(StepRanks==5))
     'STAT ON THIRD EDGE IN TRIANGLE'
     if DisplayFlag
         TRIANGLE_STAT(Species,ChipRank,NewPs,PsBy,NetNb)     
     end
     MemNewPs=NewPs;
     for FreqL=1:length(NetFrequencies)         
         TestLimit=max(1,round(NetNb*NetFrequencies(FreqL)/100));         
         StatFileName=sprintf('m%u_n%u_netnb%u_probenb%u_newps_stat_netprc%03u_pvcorr%u',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,NetFrequencies(FreqL),PvCorrRanks);         
         %COMPLETE THE INFORMATION
         %CATEGORIZE THE PS ACCORDING TO THE NUMBER OF GENES TARGETED
         'STAT'
         [NewPs,Stat,LinkedPs]=STAT(Species,ChipRank,MemNewPs,PsBy,NetNb,ProbeNb,TestLimit,SumFlag,DisplayFlag);
         cd(DirName)
         eval(sprintf('save %s NewPs Genes  PsBy Stat LinkedPs',StatFileName))
     end
end

%% STEP6 CONSTRUCT PSMATRIX
if ~isempty(find(StepRanks==6))
    for FreqL=1:length(NetFrequencies)        
        StatFileName=sprintf('m%u_n%u_netnb%u_probenb%u_newps_stat_netprc%03u_pvcorr%u',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,NetFrequencies(FreqL),PvCorrRanks);        
        cd(DirName)
        eval(sprintf('load %s',StatFileName))
        %COMPLETE THE INFORMATION
        'SUMMARIZE'
        [NewPs,PsMatrix,PsType]=SUMMARIZE(ProbeNb,NewPs,Stat);
        cd(DirName)
        eval(sprintf('save %s NewPs Genes PsBy Stat LinkedPs PsType PsMatrix',StatFileName))
    end
end

%% STEP7 FIGURES 51 TO 53
if ~isempty(find(StepRanks==7))

    FigRank=51;
    Colors=colors(colormap,5);
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u - m%u -STAT CORR',FigRank,ChipRank))

    for FreqL=1:length(NetFrequencies)        
        StatFileName=sprintf('m%u_n%u_netnb%u_probenb%u_newps_stat_netprc%03u_pvcorr%u',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,NetFrequencies(FreqL),PvCorrRanks);
        cd(DirName)
        eval(sprintf('load %s',StatFileName))
        Res={};
        GeneNb=zeros(1,3);
        for ClassL=1:6
            Pos=find(PsMatrix(:,17)==ClassL);
            GeneNb(ClassL)=length(unique(PsMatrix(Pos,1)));
        end
        %process genes targeted by several probe sets
        for ClassL=3:6
            Pos=find(PsMatrix(:,17)==ClassL);
            PsRanks=unique(PsMatrix(Pos,9));
            Res{ClassL}={};
            for PsL=1:length(PsRanks)
                Pos=find(PsMatrix(:,9)==PsRanks(PsL));
                %nb of got
                GotNb=PsMatrix(Pos(1),8);
                %got used
                Got=unique(PsMatrix(Pos,10));
                %distribution of ps in got
                Dist=histc(PsMatrix(Pos,10),Got);
                %distribution of ps grp size
                Dist=histc(Dist,[1:max(Dist)]);
                if size(Dist,1)<size(Dist,2)
                    Dist=Dist';
                end
                if length(Res{ClassL})<GotNb
                    Res{ClassL}{GotNb}=Dist;
                else
                    if length(Res{ClassL}{GotNb})>=length(Dist)                        
                        Res{ClassL}{GotNb}(1:length(Dist))=Res{ClassL}{GotNb}(1:length(Dist))+Dist;
                    else
                        Res{ClassL}{GotNb}=[Res{ClassL}{GotNb};repmat(0,length(Dist)-length(Res{ClassL}{GotNb}),1)];                        
                        Res{ClassL}{GotNb}=Res{ClassL}{GotNb}+Dist;                                                
                    end
                end
            end
        end
        Sum{FreqL}={};
        for ClassL=3:6
            Sum{FreqL}{ClassL}=Res{ClassL}{1};
            if length(Res{ClassL})>1
                for NbL=2:length(Res{ClassL})
                    if length(Sum{FreqL}{ClassL})==length(Res{ClassL}{NbL})
                        Sum{FreqL}{ClassL}=Sum{FreqL}{ClassL}+Res{ClassL}{NbL};
                    elseif length(Sum{FreqL}{ClassL})<length(Res{ClassL}{NbL})
                        Sum{FreqL}{ClassL}=[Sum{FreqL}{ClassL};zeros(length(Res{ClassL}{NbL})-length(Sum{FreqL}{ClassL}),1)];
                        Sum{FreqL}{ClassL}=Sum{FreqL}{ClassL}+Res{ClassL}{NbL};
                    else
                        Sum{FreqL}{ClassL}=Sum{FreqL}{ClassL}+[Res{ClassL}{NbL};zeros(length(Sum{FreqL}{ClassL})-length(Res{ClassL}{NbL}),1)];
                    end
                end
            end
        end

        load(sprintf('m%u_pearson_probenb%u',ChipRank,ProbeNbLimit))
        Pearson1=Pearson;
        load(sprintf('m%u_pearson_probenb%u_notfound',ChipRank,ProbeNbLimit))
        %recover information on corr
        %PsMatrix:
        %  1: rank of the assigned gene to the current probe set
        %  2: position of the assigned gene in NewPs.geneNames
        %  3: target type of assigned gene
        %  4: source type of assigned gene
        %  5: nb of probe targetting the assigned gene
        %  6: nb of not assigned genes targetted with the same nb of probes
        %  7: nb of not assigned genes targetted with less nb of probes
        %  8: nb of groups of transcripts corresponding to the assigned gene
        %  9: rank of the parent probe set
        % 10: Rank of the group of transcripts targetted by the current
        %     probe set in the assigned gene
        % 11: Rank of the group(s) of transcripts targetted by the current
        %     probe set, if it is a pivot
        % 12: [0,1] indicates if the current probe set is a pivot
        % 13: [0,1] indicates if the current probe set is paired with a pivot
        % 14: nb of probe sets that do not target the assigned gene but target
        %     a common gene with the current probe set
        % 15: nb of genes that are targetted by probe sets that do not target the
        %     assigned gene
        % 16: nb of genes that are targetted by other probe set with a nb
        %     of probes higher than the number of probes of the current probe set that
        %     target the assigned gene       
        % 17: ClassRank
        % 18: and beyond: nb of  other genes targetted with all possible nb of probes

        GeneNb=zeros(7,1);
        for ClassL=1:7
            GeneNb(ClassL)=length(unique(PsMatrix(find(PsMatrix(:,17)==ClassL),1)));
        end
        if GeneNb(7)==1
            GeneNb(7)=0;
        else
            h=errordlg(sprintf('class 7 has %u genes !',geneNb(7)));
            waitfor(h)
            error('process canceled')
        end

        InCorr{1}=cell(1,6);
        OutCorr{1}=cell(1,6);
        NotFoundCorr{1}=cell(1,6);
        NotFound{1}=cell(1,6);
        InCorr{2}=cell(1,6);
        OutCorr{2}=cell(1,6);
        NotFoundCorr{2}=cell(1,6);
        NotFound{2}=cell(1,6);
        for ClassL=3:6
            ClassPos=find(PsMatrix(:,17)==ClassL);
            % multi targeting probe sets
            PsPos=ClassPos(find(PsMatrix(ClassPos,9)));
            if ~isempty(PsPos)
                ParentRank=unique(PsMatrix(PsPos,9));
                for PsL=1:length(ParentRank)
                    %current Ps family
                    CurrPsPos=find(PsMatrix(:,9)==ParentRank(PsL));
                    %current group of transcripts
                    GotList=unique(PsMatrix(CurrPsPos,10));
                    %recover corr between probe sets that do not target the same group of
                    %transcripts
                    if length(GotList)>1
                        for GotL1=1:length(GotList)-1
                            %find the probe set that target the first group of transcripts
                            PsRanks1=find(PsMatrix(:,9)==ParentRank(PsL)&PsMatrix(:,10)==GotList(GotL1));
                            for GotL2=GotL1+1:length(GotList)
                                %find the probe set that target the snd group of transcripts
                                PsRanks2=find(PsMatrix(:,9)==ParentRank(PsL)&PsMatrix(:,10)==GotList(GotL2));
                                %recover corr for pair of probe sets
                                for PsL1=1:length(PsRanks1)
                                    PsRank1=PsRanks1(PsL1);
                                    for PsL2=1:length(PsRanks2)
                                        PsRank2=PsRanks2(PsL2);
                                        CorrPos=find(Pearson1{1}.firstPsRank==min(PsRank1,PsRank2)&Pearson1{1}.sndPsRank==max(PsRank1,PsRank2));
                                        if ~isempty(CorrPos)
                                            InCorr{2}{ClassL}=[InCorr{2}{ClassL};round(100*Pearson1{1}.rankCorr(CorrPos))];
                                        else
                                            CorrPos=find(Pearson1{2}{1}.firstPsRank==min(PsRank1,PsRank2)&Pearson1{2}{1}.sndPsRank==max(PsRank1,PsRank2));
                                            if ~isempty(CorrPos)
                                                OutCorr{2}{ClassL}=[OutCorr{2}{ClassL};round(100*Pearson1{2}{1}.rankCorr(CorrPos))];
                                            else
                                                CorrPos=find(Pearson1{2}{2}.firstPsRank==min(PsRank1,PsRank2)&Pearson1{2}{2}.sndPsRank==max(PsRank1,PsRank2));
                                                if ~isempty(CorrPos)
                                                    OutCorr{2}{ClassL}=[OutCorr{2}{ClassL};round(100*Pearson1{2}{2}.rankCorr(CorrPos))];
                                                else
                                                    CorrPos=find(Pearson{1}.firstPsRank==min(PsRank1,PsRank2)&Pearson{1}.sndPsRank==max(PsRank1,PsRank2));
                                                    if ~isempty(CorrPos)
                                                        NotFoundCorr{2}{ClassL}=[NotFoundCorr{2}{ClassL};round(100*Pearson{1}.rankCorr(CorrPos))];
                                                    else
                                                        NotFound{2}{ClassL}=[NotFound{2}{ClassL};[PsRank1,PsRank2]];
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    %recover corr between probe sets that target the same group of transcripts
                    for GotL1=1:length(GotList)-1
                        PsRanks1=find(PsMatrix(:,9)==ParentRank(PsL)&PsMatrix(:,10)==GotList(GotL1));
                        if length(PsRanks1)>1
                            for PsL1=1:length(PsRanks1)-1
                                PsRank1=PsRanks1(PsL1);
                                for PsL2=PsL1+1:length(PsRanks1)
                                    PsRank2=PsRanks1(PsL2);
                                    CorrPos=find(Pearson1{1}.firstPsRank==min(PsRank1,PsRank2)&Pearson1{1}.sndPsRank==max(PsRank1,PsRank2));
                                    if ~isempty(CorrPos)
                                        InCorr{1}{ClassL}=[InCorr{1}{ClassL};round(100*Pearson1{1}.rankCorr(CorrPos))];
                                    else
                                        CorrPos=find(Pearson1{2}{1}.firstPsRank==min(PsRank1,PsRank2)&Pearson1{2}{1}.sndPsRank==max(PsRank1,PsRank2));
                                        if ~isempty(CorrPos)
                                            OutCorr{1}{ClassL}=[OutCorr{1}{ClassL};round(100*Pearson1{2}{1}.rankCorr(CorrPos))];
                                        else
                                            CorrPos=find(Pearson1{2}{2}.firstPsRank==min(PsRank1,PsRank2)&Pearson1{2}{2}.sndPsRank==max(PsRank1,PsRank2));
                                            if ~isempty(CorrPos)
                                                OutCorr{1}{ClassL}=[OutCorr{1}{ClassL};round(100*Pearson1{2}{2}.rankCorr(CorrPos))];
                                            else
                                                CorrPos=find(Pearson{1}.firstPsRank==min(PsRank1,PsRank2)&Pearson{1}.sndPsRank==max(PsRank1,PsRank2));
                                                if ~isempty(CorrPos)
                                                    NotFoundCorr{1}{ClassL}=[NotFoundCorr{1}{ClassL};round(100*Pearson{1}.rankCorr(CorrPos))];
                                                else
                                                    NotFound{1}{ClassL}=[NotFound{1}{ClassL};[PsRank1,PsRank2]];
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        %% FIG51

        Good=[];
        Bad=[];
        for i=1:4
            Good=[Good;InCorr{1}{2+i}];
            Good=[Good;OutCorr{1}{2+i}];
            Good=[Good;NotFoundCorr{1}{2+i}];
            Bad=[Bad;InCorr{2}{2+i}];
            Bad=[Bad;OutCorr{2}{2+i}];
            Bad=[Bad;NotFoundCorr{2}{2+i}];
        end

        subplot(1,2,1)
        Val=histc(Good,[-100:100]);
        hold on
        if ~isempty(Good)
            plot([-100:100],Val*100/length(Good),'color',Colors(FreqL,:))
        end
        subplot(1,2,2)
        Val=histc(Bad,[-100:100]);
        hold on
        if ~isempty(Bad)
            plot([-100:100],Val*100/length(Bad),'color',Colors(FreqL,:))
        end
    end



    subplot(1,2,1)
    set(gca,'box','on')
    ylabel('density')
    xlabel('correlation')
    title ('ps targeting the same group')
    Legend={};
    for FreqL=1:5
        Legend{end+1}=sprintf('net frequency %03u%%',NetFrequencies(FreqL));
    end
    legend(Legend,'location','northwest')

    subplot(1,2,2)
    set(gca,'box','on')
    ylabel('density')
    xlabel('correlation')
    title ('ps targeting different groups')

    Position=get(gcf,'position');
    Position(3:4)=[700,300];
    set(gcf,'position',Position)
    try
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    catch
        mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    end
    saveas(h,sprintf('m%u_fig%u_n%03u.png',ChipRank,FigRank,NetFrequencies(FreqL)),'png')
    close(h)

    CatName={'mu','mm','cx','hx'};
    %% FIG52
    FigRank=52;
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u - m%u -SIZE OF GROUPED PROBE SETS BY CLASS',FigRank,ChipRank))
    for ClassL=3:6
        %find maximum size of probe set groups
        MaxSize=0;
        for FreqL=1:length(NetFrequencies)
            MaxSize=max(MaxSize,length(Sum{FreqL}{ClassL}));
        end
        %Construct matrix
        GrpNb=zeros(MaxSize,length(NetFrequencies));
        for FreqL=1:length(NetFrequencies)
            GrpNb(1:length(Sum{FreqL}{ClassL}),FreqL)=Sum{FreqL}{ClassL};
        end
        Colors=colors(colormap,MaxSize);
        subplot(2,2,ClassL-2)
        plot(NetFrequencies,GrpNb')
        hold on
        plot(NetFrequencies,GrpNb','+')
        xlabel('net frequency')
        ylabel('probe set group frequency')
        title(sprintf('class %s',CatName{ClassL-2}))
    end
    Position=get(gcf,'position');
    Position(3:4)=[800,600];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
    close(h)

    %% FIG53
    FigRank=53;
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u - m%u -SIZE OF GROUPED PROBE SETS',FigRank,ChipRank))
    MaxSize=0;
    for ClassL=3:6
        %find maximum size of probe set groups

        for FreqL=1:length(NetFrequencies)
            MaxSize=max(MaxSize,length(Sum{FreqL}{ClassL}));
        end
    end
    %Construct matrix
    GrpNb=zeros(MaxSize,length(NetFrequencies));
    for ClassL=3:6
        for FreqL=1:length(NetFrequencies)
            GrpNb(1:length(Sum{FreqL}{ClassL}),FreqL)=GrpNb(1:length(Sum{FreqL}{ClassL}),FreqL)+Sum{FreqL}{ClassL};
        end
    end
    Colors=colors(colormap,MaxSize);
    %         XPlot=ceil(sqrt(MaxSize));
    %         YPlot=round(MaxSize/XPlot);
    %         if XPlot*YPlot<MaxSize
    %             XPlot=XPlot+1;
    %         end
    UsedSize=min(6,size(GrpNb,1));
    XPlot=ceil(sqrt(UsedSize));
    YPlot=round(UsedSize/XPlot);
    if XPlot*YPlot<UsedSize
        XPlot=XPlot+1;
    end
    for SizeL=1:UsedSize
        subplot(YPlot,XPlot,SizeL)       
        plot(NetFrequencies,GrpNb(SizeL,:),'color',Colors(SizeL,:))
        hold on
        plot(NetFrequencies,GrpNb(SizeL,:),'+','color',Colors(SizeL,:))
        title(sprintf('size: %u',SizeL))
        xlabel('net frequency')
        ylabel('probe set group frequency')
    end
     Position=get(gcf,'position');
    Position(3:4)=[700,400];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
    close(h)
end

%% STEP8 WRITE TXT FILES
if ~isempty(find(StepRanks==8))

    if AceFlag
        %load file which make correspondance between Ensembl and AceView genes
        cd(fullfile(K.dir.pydata,Species,'txt'))
        [AceId,EnsId]=textread(sprintf('./%s_ens_by_ace_gene.txt',Species),'%s%s');
        EnsGeneName={};
        AceGeneName={};
        for AceL=1:length(AceId)
            %don't use Ensembl gene name without corresponding AceView gene
            if ~isequal(AceId{AceL},'0')
                CurrEnsId=eval(EnsId{AceL});
                if ~isempty(CurrEnsId)
                    %one AceView gene can correspond to several Ensembl genes
                    for EnsL=1:length(CurrEnsId)
                        AceGeneName{end+1,1}=AceId{AceL};
                        EnsGeneName{end+1,1}=CurrEnsId{EnsL};
                    end
                end
            end
        end
        AceFlag=1;
    end


    %load gene name and their definition corresponding to Ensembl ID
    cd(K.dir.rawdata)
    [EnsId,Def,GeneName]=textread(fullfile('.',sprintf('%s_ensembl.txt',Species)),'%s%s%s','delimiter','\t');

    %load information files listing targeted transcripts and exons for Ensembl genes
    cd(DirName)
    PsRanks={};
    eval(sprintf('load m%u_probesets_by_ensembl_gene',ChipRank))
    GeneNames{1}=EGeneName;
    PsRanks{1}=EPsRanks;
    %TargetingPsRanks{1}=ETargetingPsRanks;
    ExonNames{1}=EExonNames;
    ExonRanks{1}=EExonRanks;
    ExonProbeNbs{1}=EExonProbeNbs;
    LastExons{1}=ELastExon;
    LastGroups{1}=ELastGroup;
    GroupRanks{1}=EGroupRanks;
    GroupProbeNbs{1}=EGroupProbeNbs;
    TargetedTs{1}=ETargetedTs;
    NotTargetedTs{1}=ENotTargetedTs;
    clear EGeneName EGeneNames ETargetedGenes ETargetingPsRanks EPsNames EPsRanks EExonNames EExonProbeNbs EGroupRanks EGroupProbeNbs ETargetedTs ENotTargetedTs

    if AceFlag
        %load information files listing targeted transcripts and exons for AceView genes
        eval(sprintf('load m%u_probesets_by_aceview_gene',ChipRank))
        GeneNames{2}=AGeneName;
        PsRanks{2}=APsRanks;
        %TargetingPsRanks{2}=ATargetingPsRanks;
        ExonNames{2}=AExonNames;
        ExonRanks{2}=AExonRanks;
        ExonProbeNbs{2}=AExonProbeNbs;
        LastExons{2}=ALastExon;
        LastGroups{2}=ALastGroup;
        GroupRanks{2}=AGroupRanks;
        GroupProbeNbs{2}=AGroupProbeNbs;
        TargetedTs{2}=ATargetedTs;
        NotTargetedTs{2}=ANotTargetedTs;
        clear AGeneName AGeneNames ATargetedGenes ATargetingPsRanks APsNames APsRanks AExonNames AExonProbeNbs AGroupRanks AGroupProbeNbs ATargetedTs ANotTargetedTs
    end


    %load probe set names
    cd(K.dir.rawdata)
    PsName=textread(sprintf('m%u_probeset.txt',ChipRank),'%s');
    AllPsPair=[];
    PsPair=[];
    GeneId={};
    NetFrequencies=[1,25,50,75,100];
    cd(DirName)
    %find GeneName
    PsGeneName=cell(size(Genes.name,1),1);
    for GeneL=1:size(Genes.name,1)
        CurrGeneId=Genes.name{GeneL};
        EnsPos=strmatch(CurrGeneId,EnsId,'exact');
        if ~isempty(EnsPos)
            PsGeneName{GeneL}=GeneName{EnsPos};
        else
            PsGeneName{GeneL}=CurrGeneId;
        end
    end

    PsNb={};
    OrphPs={};
    for FreqL=1:length(NetFrequencies)
        StatFileName=sprintf('m%u_n%u_netnb%u_probenb%u_newps_stat_netprc%03u_pvcorr%u',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,NetFrequencies(FreqL),PvCorrRanks);
        eval(sprintf('load %s',StatFileName))

        %write PsMatrix   
        fid=fopen(sprintf('m%u_n%u_netnb%u_probenb%u_netprc%03u_pvcorr%u_psmatrix.txt',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,NetFrequencies(FreqL),PvCorrRanks),'w');
        Out=[repmat('%s\t',1,3),repmat('%u\t',1,size(PsMatrix,2)-1),'%u\n'];
        for PsL=1:size(PsMatrix,1)
            if PsMatrix(PsL,1)>0
                fprintf(fid,Out,PsName{PsL},Genes.name{PsMatrix(PsL,1)},PsGeneName{PsMatrix(PsL,1)},PsMatrix(PsL,:));
            else
                fprintf(fid,Out,PsName{PsL},'-','-',PsMatrix(PsL,:));
            end
        end
        fclose(fid)

        PsMade=zeros(size(PsMatrix,1),1);
        %process only classes where exist probe set pairs
        Pos=find(PsMatrix(:,17)>=3& PsMatrix(:,17)<7);  
        %probe set without redundancy (possible for class CX and HX since assignation does
        %not necessarily conserve redundancy) => select class MM probe sets
        OrphPs{FreqL}=[];
        for PsL=1:length(Pos)
            if PsMade(Pos(PsL))==0
                %recover all the probe sets related to the parent of the currently proccessed
                %probe set => allow to determine if exist similar probe set pairs
                CurrPsRanks=find(PsMatrix(:,9)==PsMatrix(Pos(PsL),9));        
                %scan all the possible pairs of these probe sets
                if length(CurrPsRanks)>1
                    for PsL1=1:length(CurrPsRanks)-1
                        PsRank1=CurrPsRanks(PsL1);
                        IsPivot1=0;
                        if PsMatrix(PsRank1,12)==1
                            %groups of transcripts targetted by the first probe set if it is a
                            %pivot probe set
                            Grps1=find(dec2bin(PsMatrix(PsRank1,11))=='1');
                            IsPivot1=1;
                        end
                        for PsL2=PsL1+1:length(CurrPsRanks)
                            PsRank2=CurrPsRanks(PsL2);
                            %test if the current pair exist
                            try
                                PairPos=find(PsPair(:,1)==min(PsRank1,PsRank2)&PsPair(:,2)==max(PsRank1,PsRank2));
                            catch
                                PairPos=[];
                            end
                            if isempty(PairPos)
                                PairPos=size(PsPair,1)+1;
                                GeneId{PairPos}=Genes.name{PsMatrix(Pos(PsL),1)};
                                PsPair(PairPos,:)=[min(PsRank1,PsRank2),max(PsRank1,PsRank2),zeros(1,19)];
                                PsPair(PairPos,13)=PsMatrix(PsRank1,17);
                                PsPair(PairPos,14)=PsMatrix(PsRank1,5);
                                PsPair(PairPos,15)=PsMatrix(PsRank1,6);
                                PsPair(PairPos,16)=PsMatrix(PsRank1,7);
                                PsPair(PairPos,17)=PsMatrix(PsRank2,5);
                                PsPair(PairPos,18)=PsMatrix(PsRank2,6);
                                PsPair(PairPos,19)=PsMatrix(PsRank2,7);
                                PsPair(PairPos,20)=PsMatrix(PsRank1,3);
                                PsPair(PairPos,21)=PsMatrix(PsRank2,3);
                            end
                            %Presence of ps pair in the current frequency
                            PsPair(PairPos,2+(FreqL-1)*2+1)=1;
                            IsPivot2=0;
                            if PsMatrix(PsRank2,12)==1
                                %groups of transcripts targetted by the second probe set if it
                                %is a pivot probe set
                                Grps2=find(dec2bin(PsMatrix(PsRank2,11))=='1');
                                IsPivot2=1;
                            end
                            %find if the two current probe set are a pair which target the same
                            %group of transcripts
                            if IsPivot1==0&IsPivot2==0
                                if PsMatrix(PsRank1,10)==PsMatrix(PsRank2,10)
                                    PsPair(PairPos,2+(FreqL-1)*2+2)=1;
                                end
                            elseif IsPivot1==1 & IsPivot2==0
                                if ~isempty(find(Grps1==PsMatrix(CurrPsRanks(PsL2),10)))
                                    PsPair(PairPos,2+(FreqL-1)*2+2)=2;
                                end
                            elseif IsPivot2==1 & IsPivot1==0
                                if ~isempty(find(Grps2==PsMatrix(CurrPsRanks(PsL1),10)))
                                    PsPair(PairPos,2+(FreqL-1)*2+2)=2;
                                end
                            else
                                if ~isempty(intersect(Grps1,Grps2))
                                    PsPair(PairPos,2+(FreqL-1)*2+2)=3;
                                end
                            end
                        end
                    end
                else
                    if PsMatrix(Pos(PsL),17)<5
                        OrphPs{FreqL}(end+1,1)=Pos(PsL);
                    end
                end                
                PsMade(CurrPsRanks)=1;
            end
        end
    end

    %write results (PsPair)
    %Order on Gene ID
    [GeneId,SortOrder]=sort(GeneId);
    PsPair=PsPair(SortOrder,:);
    %mark pair not present in all net frequencies (~1%)
    MarkPos=find(PsPair(:,3)==0|PsPair(:,5)==0|PsPair(:,7)==0|PsPair(:,9)==0|PsPair(:,11)==0);
    Mark=repmat('v',size(PsPair,1),1);
    Mark(MarkPos)='*';
    %don't keep columns 3:2:11 information on similarity is in columns 4:2:12
    PsPair=PsPair(:,[1,2,4:2:12,13:21]);


    %complete PsPair
    %recover information on transcripts and exons
    if AceFlag
        TypeNb=2;
    else
        TypeNb=1;
    end

    for TypeL=1:TypeNb
        TargetInfo=zeros(size(PsPair,1),26);
        for PairL=1:size(PsPair,1)
            if isempty(strmatch('GOP',GeneId{PairL}))
                %find the number of targeting probes to interrogate the right cell files in
                %information files
                ProbeNb1=PsPair(PairL,9)+1;
                ProbeNb2=PsPair(PairL,12)+1;
                %find the position of the targeted gene in information files
                GenePos1=[];
                if TypeL==1
                    if ~isempty(strmatch('ENS',GeneId{PairL}))|~isempty(strmatch('GOP',GeneId{PairL}))
                        GenePos1=strmatch(GeneId{PairL},GeneNames{TypeL}{ProbeNb1},'exact');
                        if isempty(GenePos1)
                            GenePos1=strmatch(GeneId{PairL},GeneNames{TypeL}{1},'exact');
                            ProbeNb1=1;
                        end
                    end
                else
                    if ~isempty(strmatch('ENS',GeneId{PairL}))
                        EnsPos=strmatch(GeneId{PairL},EnsGeneName,'exact');
                        if ~isempty(EnsPos)
                            %use the first AceId if exists several
                            AceId=AceGeneName{EnsPos(1)};
                            GenePos1=strmatch(AceId,GeneNames{TypeL}{ProbeNb1},'exact');
                            if isempty(GenePos1)
                                GenePos1=strmatch(AceId,GeneNames{TypeL}{1},'exact');
                                ProbeNb1=1;
                            end
                        end
                    else
                        if isempty(strmatch('GOP',GeneId{PairL}))
                            GenePos1=strmatch(GeneId{PairL},GeneNames{TypeL}{ProbeNb1},'exact');
                            if isempty(GenePos1)
                                GenePos1=strmatch(GeneId{PairL},GeneNames{TypeL}{1},'exact');
                                ProbeNb1=1;
                            end
                        end
                    end
                end

                if ~isempty(GenePos1)
                    %recover the position of the current Ps1
                    PsPos1=find(PsRanks{TypeL}{ProbeNb1}{GenePos1}==PsPair(PairL,1));
                    %recover information on transcripts and exons targeted by Ps1
                    if ~isempty(PsPos1)
                        Ts1=TargetedTs{TypeL}{ProbeNb1}{GenePos1}{PsPos1};
                        NotTs1=NotTargetedTs{TypeL}{ProbeNb1}{GenePos1}{PsPos1};
                        ExRank1=ExonRanks{TypeL}{ProbeNb1}{GenePos1}{PsPos1};
                        Ex1=ExonProbeNbs{TypeL}{ProbeNb1}{GenePos1}{PsPos1};
                        LastExon1=LastExons{TypeL}{ProbeNb1}(GenePos1);
                        LastGroup1=LastGroups{TypeL}{ProbeNb1}(GenePos1);
                        GrpRank1=GroupRanks{TypeL}{ProbeNb1}{GenePos1}{PsPos1};
                        Grp1=GroupProbeNbs{TypeL}{ProbeNb1}{GenePos1}{PsPos1};
                    else
                        Ts1=[];
                        NotTs1=[];
                        ExRank1=[];
                        Ex1=[];
                        LastExon1=[];
                        LastGroup1=[];
                        GrpRank1=[];
                        Grp1=[];
                    end
                    GenePos2=[];
                    if TypeL==1
                        if ~isempty(strmatch('ENS',GeneId{PairL}))|~isempty(strmatch('GOP',GeneId{PairL}))
                            GenePos2=strmatch(GeneId{PairL},GeneNames{TypeL}{ProbeNb2},'exact');
                            if isempty(GenePos2)
                                GenePos2=strmatch(GeneId{PairL},GeneNames{TypeL}{1},'exact');
                                ProbeNb2=1;
                            end
                        end
                    else
                        if ~isempty(strmatch('ENS',GeneId{PairL}))
                            EnsPos=strmatch(GeneId{PairL},EnsGeneName,'exact');
                            if ~isempty(EnsPos)
                                %use the first AceId if exists several
                                AceId=AceGeneName{EnsPos(1)};
                                GenePos2=strmatch(AceId,GeneNames{TypeL}{ProbeNb2},'exact');
                                if isempty(GenePos2)
                                    GenePos2=strmatch(AceId,GeneNames{TypeL}{1},'exact');
                                    ProbeNb2=1;
                                end
                            end
                        else
                            if isempty(strmatch('GOP',GeneId{PairL}))
                                GenePos2=strmatch(GeneId{PairL},GeneNames{TypeL}{ProbeNb2},'exact');
                                if isempty(GenePos2)
                                    GenePos2=strmatch(GeneId{PairL},GeneNames{TypeL}{1},'exact');
                                    ProbeNb2=1;
                                end
                            end
                        end
                    end
                    if ~isempty(GenePos2)
                        %recover the position of the current Ps1
                        PsPos2=find(PsRanks{TypeL}{ProbeNb2}{GenePos2}==PsPair(PairL,2));
                        %recover information on transcripts and exons targeted by Ps1
                        if ~isempty(PsPos2)
                            Ts2=TargetedTs{TypeL}{ProbeNb2}{GenePos2}{PsPos2};
                            NotTs2=NotTargetedTs{TypeL}{ProbeNb2}{GenePos2}{PsPos2};
                            ExRank2=ExonRanks{TypeL}{ProbeNb2}{GenePos2}{PsPos2};
                            Ex2=ExonProbeNbs{TypeL}{ProbeNb2}{GenePos2}{PsPos2};
                            LastExon2=LastExons{TypeL}{ProbeNb2}(GenePos2);
                            LastGroup2=LastGroups{TypeL}{ProbeNb2}(GenePos2);
                            GrpRank2=GroupRanks{TypeL}{ProbeNb2}{GenePos2}{PsPos2};
                            Grp2=GroupProbeNbs{TypeL}{ProbeNb2}{GenePos2}{PsPos2};
                        else
                            Ts2=[];
                            NotTs2=[];
                            ExRank2=[];
                            Ex2=[];
                            LastExon2=[];
                            LastGroup2=[];
                            GrpRank2=[];
                            Grp2=[];
                        end
                    else
                        Ts2=[];
                        NotTs2=[];
                        ExRank2=[];
                        Ex2=[];
                        LastExon2=[];
                        LastGroup2=[];
                        GrpRank2=[];
                        Grp2=[];
                    end
                else
                    Ts1=[];
                    NotTs1=[];
                    ExRank1=[];
                    Ex1=[];
                    LastExon1=[];
                    LastGroup1=[];
                    GrpRank1=[];
                    Grp1=[];
                    Ts2=[];
                    NotTs2=[];
                    ExRank2=[];
                    Ex2=[];
                    LastExon2=[];
                    LastGroup2=[];
                    GrpRank2=[];
                    Grp2=[];
                end
                if~isempty(GenePos1)&~isempty(GenePos2)
                    %fill TargetInfo
                    %ps ranks
                    TargetInfo(PairL,1)=PsPair(PairL,1);
                    TargetInfo(PairL,2)=PsPair(PairL,2);
                    %class
                    TargetInfo(PairL,3)=PsPair(PairL,8);
                    %total number of transcripts
                    TargetInfo(PairL,4)=length(union(Ts1,union(Ts2,union(NotTs1,NotTs2))));
                    %nb of common trancripts
                    TargetInfo(PairL,5)=length(intersect(Ts1,Ts2));
                    %nb of specific transcript
                    TargetInfo(PairL,6)=length(setdiff(Ts1,Ts2));
                    TargetInfo(PairL,7)=length(setdiff(Ts2,Ts1));
                    %Overlaping index for transcripts
                    if isempty(Ts1)&isempty(Ts2)
                        TargetInfo(PairL,16)=-1;
                    else
                        TargetInfo(PairL,16)=round(TargetInfo(PairL,5)*100/(TargetInfo(PairL,5)+TargetInfo(PairL,6)+TargetInfo(PairL,7)));
                        if isnan(TargetInfo(PairL,16))|TargetInfo(PairL,16)>100
                            'stop'
                        end
                    end
                    %nb of targeted exons
                    TargetInfo(PairL,8)=length(union(ExRank1,ExRank2));
                    %nb of probes in common exons
                    [Temp,ComPos1,ComPos2]=intersect(ExRank1,ExRank2);
                    TargetInfo(PairL,9)=sum(Ex1(ComPos1))+sum(Ex2(ComPos2));
                    %nb of probe specific to uncomon exons
                    [Temp,DiffPos]=setdiff(ExRank1,ExRank2);
                    TargetInfo(PairL,10)=sum(Ex1(DiffPos));
                    [Temp,DiffPos]=setdiff(ExRank2,ExRank1);
                    TargetInfo(PairL,11)=sum(Ex2(DiffPos));
                    %Overlaping index for exons
                    %total number of probes in common and specific exon(s)
                    PNb=TargetInfo(PairL,9)+TargetInfo(PairL,10)+TargetInfo(PairL,11);
                    %list of all targeted exons
                    ExonRank=union(ExRank1,ExRank2);
                    %weighted overlap in each targeted exon
                    Overlap=[];
                    if PNb>0
                        for ExonL=1:length(ExonRank)
                            %number of probes of the first probeset in the current exon
                            Pos1=find(ExRank1==ExonRank(ExonL));
                            PNb1=0;
                            if ~isempty(Pos1)
                                PNb1=Ex1(Pos1);
                            end
                            %number of probes of the second probeset in the current exon
                            Pos2=find(ExRank2==ExonRank(ExonL));
                            PNb2=0;
                            if ~isempty(Pos2)
                                PNb2=Ex2(Pos2);
                            end
                            %weighted overlap for the current exon
                            Overlap(end+1)=((PNb1+PNb2)/PNb)*(min(PNb1,PNb2)/max(PNb1,PNb2));
                        end
                        Overlap=round(100*sum(Overlap)/length(Overlap));
                        TargetInfo(PairL,17)=Overlap;
                    else
                        TargetInfo(PairL,17)=-1;
                    end

                    %nb of groups
                    TargetInfo(PairL,12)=length(union(GrpRank1,GrpRank2));
                    %nb of probes in common groups
                    [Temp,ComPos1,ComPos2]=intersect(GrpRank1,GrpRank2);
                    TargetInfo(PairL,13)=sum(Grp1(ComPos1))+sum(Grp2(ComPos2));
                    %nb of probe specific to uncomon groups
                    [Temp,DiffPos]=setdiff(GrpRank1,GrpRank2);
                    TargetInfo(PairL,14)=sum(Grp1(DiffPos));
                    [Temp,DiffPos]=setdiff(GrpRank2,GrpRank1);
                    TargetInfo(PairL,15)=sum(Grp2(DiffPos));
                    GrpRank=union(GrpRank1,GrpRank2);
                    %Overlaping index for groups
                    PNb=TargetInfo(PairL,13)+TargetInfo(PairL,14)+TargetInfo(PairL,15);
                    Overlap=[];
                    if PNb>0
                        for GrpL=1:length(GrpRank)
                            Pos1=find(GrpRank1==GrpRank(GrpL));
                            PNb1=0;
                            if ~isempty(Pos1)
                                PNb1=Grp1(Pos1);
                            end
                            Pos2=find(GrpRank2==GrpRank(GrpL));
                            PNb2=0;
                            if ~isempty(Pos2)
                                PNb2=Grp2(Pos2);
                            end
                            Overlap(end+1)=((PNb1+PNb2)/PNb)*(min(PNb1,PNb2)/max(PNb1,PNb2));
                        end
                        Overlap=round(100*sum(Overlap)/length(Overlap));
                        TargetInfo(PairL,18)=Overlap;
                    else
                        TargetInfo(PairL,18)=-1;
                    end

                    %percentage of probes in last exon
                    if LastExon1~=LastExon2
                        h=warndlg(sprintf('LastExon1 ~=LastExon2 for pair %u',PairL));
                        waitfor(h)
                        'stop'
                    else
                        [Temp,LastIndex,Temp]=intersect(ExRank1,LastExon1);
                        if ~isempty(LastIndex)
                            %TargetInfo(PairL,19)=round(Ex1(LastIndex)*100/sum(Ex1));
                            TargetInfo(PairL,19)=round(Ex1(LastIndex)*100/PsPair(PairL,9));
                        end
                        [Temp,LastIndex,Temp]=intersect(ExRank2,LastExon2);
                        if ~isempty(LastIndex)
                            %TargetInfo(PairL,20)=round(Ex2(LastIndex)*100/sum(Ex2));
                            TargetInfo(PairL,20)=round(Ex2(LastIndex)*100/PsPair(PairL,12));
                        end
                    end
                    %all probes in one exon
                    if length(ExRank1)==1
                        if Ex1(1)==PsPair(PairL,9)
                            TargetInfo(PairL,21)=1;
                        end
                    end
                    if length(ExRank2)==1
                        if Ex2(1)==PsPair(PairL,12)
                            TargetInfo(PairL,22)=1;
                        end
                    end

                    %percentage of probes in last group
                    if LastGroup1~=LastGroup2
                        h=warndlg(sprintf('LastGroup1 ~=LastGroup2 for pair %u',PairL));
                        waitfor(h)
                        'stop'
                    else
                        [Temp,LastIndex,Temp]=intersect(GrpRank1,LastGroup1);
                        if ~isempty(LastIndex)
                            %TargetInfo(PairL,23)=round(Grp1(LastIndex)*100/sum(Grp1));
                            TargetInfo(PairL,23)=round(Grp1(LastIndex)*100/PsPair(PairL,9));
                        end
                        [Temp,LastIndex,Temp]=intersect(GrpRank2,LastGroup2);
                        if ~isempty(LastIndex)
                            %TargetInfo(PairL,24)=round(Grp2(LastIndex)*100/sum(Grp2));
                            TargetInfo(PairL,24)=round(Grp2(LastIndex)*100/PsPair(PairL,12));
                        end
                    end
                    %all probes in one group
                    if length(GrpRank1)==1
                        if Grp1(1)==PsPair(PairL,9)
                            TargetInfo(PairL,25)=1;
                        end
                    end
                    if length(GrpRank2)==1
                        if Grp2(1)==PsPair(PairL,12)
                            TargetInfo(PairL,26)=1;
                        end
                    end
                end
            end
        end
        PsPair=[PsPair,TargetInfo(:,4:end)];
    end

    %find GeneName
    PsGeneName=cell(size(PsPair,1),1);
    for PsL=1:size(PsPair,1)
        CurrGeneId=GeneId{PsL};
        EnsPos=strmatch(CurrGeneId,EnsId,'exact');
        if ~isempty(EnsPos)
            PsGeneName{PsL}=GeneName{EnsPos};
        else
            PsGeneName{PsL}=CurrGeneId;
        end
    end

    cd(DirName)
    fid=fopen(sprintf('m%u_n%u_netnb%u_probenb%u_pvcorr%u_pspair.txt',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,PvCorrRanks),'w');
    if AceFlag
        Out=[repmat('%s\t',1,4),repmat('%u\t',1,16),'%s',repmat('\t%u',1,12),repmat('\t%d',1,3),repmat('\t%u',1,20),repmat('\t%d',1,3),repmat('\t%u',1,8),'\n'];
        for PsL=1:size(PsPair,1)
            fprintf(fid,Out,GeneId{PsL},PsGeneName{PsL},PsName{PsPair(PsL,1)},PsName{PsPair(PsL,2)},PsPair(PsL,1:16),Mark(PsL),PsPair(PsL,17:62));
        end
    else
        Out=[repmat('%s\t',1,4),repmat('%u\t',1,16),'%s',repmat('\t%u',1,12),repmat('\t%d',1,3),repmat('\t%u',1,8),'\n'];
        for PsL=1:size(PsPair,1)
            fprintf(fid,Out,GeneId{PsL},PsGeneName{PsL},PsName{PsPair(PsL,1)},PsName{PsPair(PsL,2)},PsPair(PsL,1:16),Mark(PsL),PsPair(PsL,17:39));
        end
    end
    fclose(fid)


    % FIG 56
    cd(DirName)
    PsPair=[];    
    if AceFlag
        TypeNb=2;
        [GeneId,PsGeneName,PsName1,PsName2,PsPair(:,1),PsPair(:,2),PsPair(:,3),PsPair(:,4),PsPair(:,5),...
            PsPair(:,6),PsPair(:,7),PsPair(:,8),PsPair(:,9),PsPair(:,10),PsPair(:,11),PsPair(:,12),PsPair(:,13),PsPair(:,14),...
            PsPair(:,15),PsPair(:,16),Mark,PsPair(:,17),PsPair(:,18),PsPair(:,19),PsPair(:,20),...
            PsPair(:,21),PsPair(:,22),PsPair(:,23),PsPair(:,24),PsPair(:,25),PsPair(:,26),PsPair(:,27),...
            PsPair(:,28),PsPair(:,29),PsPair(:,30),PsPair(:,31),PsPair(:,32),PsPair(:,33),PsPair(:,34),...
            PsPair(:,35),PsPair(:,36),PsPair(:,37),PsPair(:,38),PsPair(:,39),PsPair(:,40),PsPair(:,41),...
            PsPair(:,42),PsPair(:,43),PsPair(:,44),PsPair(:,45),PsPair(:,46),PsPair(:,47),PsPair(:,48),...
            PsPair(:,49),PsPair(:,50),PsPair(:,51),PsPair(:,52),PsPair(:,53),PsPair(:,54),PsPair(:,55),...
            PsPair(:,56),PsPair(:,57),PsPair(:,58),PsPair(:,59),PsPair(:,60),PsPair(:,61),PsPair(:,62)]...
            =textread(fullfile('.',sprintf('m%u_n%u_netnb%u_probenb%u_pvcorr%u_pspair.txt',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,PvCorrRanks)),...
            [repmat('%s',1,4),repmat('%u',1,16),'%c',repmat('%u',1,12),repmat('%d',1,3),repmat('%u',1,20),repmat('%d',1,3),repmat('%u',1,8)],'delimiter','\t');
    else
        TypeNb=1;
        [GeneId,PsGeneName,PsName1,PsName2,PsPair(:,1),PsPair(:,2),PsPair(:,3),PsPair(:,4),PsPair(:,5),...
            PsPair(:,6),PsPair(:,7),PsPair(:,8),PsPair(:,9),PsPair(:,10),PsPair(:,11),PsPair(:,12),PsPair(:,13),PsPair(:,14),...
            PsPair(:,15),PsPair(:,16),Mark,PsPair(:,17),PsPair(:,18),PsPair(:,19),PsPair(:,20),...
            PsPair(:,21),PsPair(:,22),PsPair(:,23),PsPair(:,24),PsPair(:,25),PsPair(:,26),PsPair(:,27),...
            PsPair(:,28),PsPair(:,29),PsPair(:,30),PsPair(:,31),PsPair(:,32),PsPair(:,33),PsPair(:,34),...
            PsPair(:,35),PsPair(:,36),PsPair(:,37),PsPair(:,38),PsPair(:,39)]...
            =textread(fullfile('.',sprintf('m%u_n%u_netnb%u_probenb%u_pvcorr%u_pspair.txt',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,PvCorrRanks)),...
            [repmat('%s',1,4),repmat('%u',1,16),'%c',repmat('%u',1,12),repmat('%d',1,3),repmat('%u',1,8)],'delimiter','\t');
    end

    % Filter probe set pairs on selected properties
    %Index of ENSEMBL genes that have a corresponding ACEVIEW GENE
    %Ensembl genes
    a=strmatch('ENS',GeneId);
    if AceFlag
        %Probeset pair that target also AceView genes
        g=find(sum(PsPair(:,40:62),2)>0);
        a=intersect(a,g);
    end
    %index of paired ps that are similar in a fixed number of networks
    %paired ps that are never similar
    b={};
    b{1}=find(PsPair(:,3)==0);
    %Paired ps that are similar in less than x% networks
    for PrcL=1:4
        b{PrcL+1}=find(PsPair(:,2+PrcL)==1&PsPair(:,3+PrcL)==0);
    end
    %Paired ps that are always similar
    b{6}=find(PsPair(:,7)==1);
    %Class MS
    c=find(PsPair(:,8)==3);
    %ps pair that target only one gene in classes CX and HX
    d=find(sum(PsPair(:,[10,11,13,14]),2)==0);
    %ps pair that target gene(s) with 11 probes
    e=find(PsPair(:,9)==11&PsPair(:,12)==11);
    %intersect of previous indexes
    f={};
    for PrcL=1:6
        f{PrcL}=intersect(a,intersect(b{PrcL},intersect(e,union(c,d))));
    end

    FigRank=56;
    Colors=colors(colormap,6);
    for TypeL=1:TypeNb
        h=figure;
        set(gcf,'color',[1,1,1])
        if TypeL==1
            Offset=13;
            set(h,'name',sprintf('FIG%ua: Ensembl statistics on redundant paired probe sets',FigRank))
        else
            Offset=36;
            set(h,'name',sprintf('FIG%ub: AceView statistics on redundant paired probe sets',FigRank))
        end
        for PrcL=1:6
            if ~isempty(f{PrcL})
                % percentage of common transcripts (=overlap index for transcripts)
                subplot(4,3,1)
                hold on               
                T=PsPair(f{PrcL},5+Offset)+PsPair(f{PrcL},6+Offset)+PsPair(f{PrcL},7+Offset);                        
                Pos=find(T>0);
                ComT=round(PsPair(f{PrcL}(Pos),5+Offset)*100./T(Pos));        
                [Value Temp]=histc(ComT,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))                                

                % percentage of uncommon transcripts
                subplot(4,3,4)
                hold on
                DiffT=round((PsPair(f{PrcL}(Pos),6+Offset)+PsPair(f{PrcL}(Pos),7+Offset))*100./T(Pos));
                [Value Temp]=histc(DiffT,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                % percentage of common probes in exons
                subplot(4,3,2)
                hold on
                E=PsPair(f{PrcL},9+Offset)+PsPair(f{PrcL},10+Offset)+PsPair(f{PrcL},11+Offset);
                Pos=find(E>0);
                ComE=round(PsPair(f{PrcL}(Pos),9+Offset)*100./E(Pos));
                [Value Temp]=histc(ComE,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                % percentage of different probes in exons
                subplot(4,3,5)
                hold on
                DiffE=round((PsPair(f{PrcL}(Pos),10+Offset)+PsPair(f{PrcL}(Pos),11+Offset))*100./E(Pos));
                [Value Temp]=histc(DiffE,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                % percentage of common probes in groups
                subplot(4,3,3)
                hold on
                G=PsPair(f{PrcL},13+Offset)+PsPair(f{PrcL},14+Offset)+PsPair(f{PrcL},15+Offset);
                Pos=find(G>0);
                ComG=round(PsPair(f{PrcL}(Pos),13+Offset)*100./G(Pos));
                [Value Temp]=histc(ComG,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))
                % percentage of different probes in groups
                subplot(4,3,6)
                hold on
                DiffG=round((PsPair(f{PrcL}(Pos),14+Offset)+PsPair(f{PrcL}(Pos),15+Offset))*100./G(Pos));
                [Value Temp]=histc(DiffG,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                % nb of targeted exons
                subplot(4,3,7)
                hold on
                T=[PsPair(f{PrcL},8+Offset)];
                [Value Temp]=histc(T,[1:10]);
                plot([1:10],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                %overlap index for exons
                subplot(4,3,8)
                hold on
                Overlap=PsPair(f{PrcL},17+Offset);
                [Value Temp]=histc(Overlap,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                %overlap index for groups
                subplot(4,3,9)
                hold on
                Overlap=PsPair(f{PrcL},18+Offset);
                [Value Temp]=histc(Overlap,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                % percentage of probes in last exon
                subplot(4,3,10)
                hold on
                T=[PsPair(f{PrcL},19+Offset);PsPair(f{PrcL},20+Offset)];
                [Value Temp]=histc(T,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))

                % percentage of probes in last group
                subplot(4,3,11)
                hold on
                T=[PsPair(f{PrcL},23+Offset);PsPair(f{PrcL},24+Offset)];
                [Value Temp]=histc(T,[0:100]);
                plot([0:100],cumsum(Value)/sum(Value),'color',Colors(PrcL,:))                
            end
        end

        for i=1:11
            subplot(4,3,i)           
            set(gca,'box','on')
            set(gca,'ylim',[0,1])

            switch i
                case 1
                    title(sprintf('percentage of transcripts\n targeted in common'))
                case 4
                    title(sprintf('percentage of transcripts\n not targeted in common'))
                case 2
                    title(sprintf('percentage of probes\n in shared exons'))
                case 5
                    title(sprintf('percentage of probes\n not in shared exons'))
                case 3
                    title(sprintf('percentage of probes\n in shared groups'))
                case 6
                    title(sprintf('percentage of probes\n not in shared groups'))
                case 7
                    title(sprintf('number of targeted exons'))
                case 8
                    title('exon overlap index')
                case 9
                    title('group overlap index')
                case 10
                    title(sprintf('percentage of probes\n in last exon'))
                case 11
                    title(sprintf('percentage of probes\n in last group'))
            end

        end
        legend('0%','1%','25%','50%','75%','100%');
        Position=get(gcf,'position');
        Position(3:4)=[1000,1000];
        set(gcf,'position',Position)
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        if TypeL==1
            saveas(h,sprintf('m%u_fig%ua.png',ChipRank,FigRank),'png')
        else
            saveas(h,sprintf('m%u_fig%ub.png',ChipRank,FigRank),'png')
        end
        close(h)
    end
end


%% LOADSIM
% INPUT PARAMETERS
% 1   ChipRank : the rank of chip set model
% 2   NetRanks : a list of network rank
% 3 PvCorrRanks: five series of p-values are calculated (on edges with
%                corr>=0,10,20,30,40,50)
% 4   FileName : the common part of or file name to be loaded

%OUTPUT
% AllVal : all CORR,ANTI and PV (p-values) of the networks
%    Sim : The last loaded Sim
function [AllVal,Sim]=LOADSIM(ChipRank,NetRanks,PvCorrRanks,FileName)

NetNb=length(NetRanks);
AllVal.corr=[];
AllVal.anti=[];
AllVal.pv=[];
for NetL=1:NetNb
    eval(sprintf('load m%u_n%u_%s',ChipRank,NetRanks(NetL),FileName));
    AllVal.corr=[AllVal.corr,Sim.corr];
    AllVal.anti=[AllVal.anti,Sim.anti];
    AllVal.pv=[AllVal.pv,Sim.pv{PvCorrRanks}];
end
%% CONSTRUCT_NEWPS

function [NewPs,PsBy,Genes]=CONSTRUCT_NEWPS(Ps,EnsGenesDict,AceGenesDict)

if isempty(AceGenesDict)
    AceFlag=0;
else
    AceFlag=1;
end
PsNb=length(Ps{1});
NewPs=[];
for PsL=1:PsNb
    NewPs{PsL,1}.geneNames={};
    NewPs{PsL,1}.probeNb={};
    NewPs{PsL,1}.source={};
    NewPs{PsL,1}.target={};
end

Genes.name={};
%Genes.source=[];
PsBy.gene={};

for PsL=1:PsNb
    %recover genes targetted by the current probeset inside exons or splice
    EnsGeneNames=Ps{1}{PsL}.geneNamesSup;
    EnsProbeNb=Ps{1}{PsL}.probeNbSup;
    if AceFlag
        AceGeneNames=Ps{2}{PsL}.geneNamesSup;
        AceProbeNb=Ps{2}{PsL}.probeNbSup;
        CurrAceGeneNames={};
        CurrAceProbeNb=[];
        %search for gene specific to AceView
        if ~isempty(AceGeneNames)
            %process each AceView gene targeted by the current probe set
            for GeneL1=1:length(AceGeneNames)
                %process genes that have AceView names
                DoIt=1;
                if length(AceGeneNames{GeneL1})>=3
                    if isequal(AceGeneNames{GeneL1}(1:3),'ENS')
                        DoIt=0;
                    end
                end
                if DoIt
                    AcePos=strmatch(AceGeneNames{GeneL1},AceGenesDict,'exact');
                    if ~isempty(EnsGenesDict{AcePos})
                        % The current AceView gene exist also in Ensembl
                        if length(EnsGenesDict{AcePos})==1
                            % It is one to one correspondance between
                            % AceView and Ensembl genes
                            %if the gene is not in the current list of targeted genes add it
                            if isempty(strmatch(EnsGenesDict{AcePos}{1},EnsGeneNames,'exact'))
                                EnsGeneNames=[EnsGeneNames,{EnsGenesDict{AcePos}{1}}];
                                EnsProbeNb=[EnsProbeNb,AceProbeNb(GeneL1)];
                            end
                        else
                            %if length >1 it is a composit gene (an AceView gene that overlap two ensembl gene
                            %add an Ensemlb gene only if none are in common
                            %otherwise would systematically add genes close to true Ensembl targets
                            AddFlag=1;
                            for GeneL2=1:length(EnsGenesDict{AcePos})
                                if ~isempty(strmatch(EnsGenesDict{AcePos}{GeneL2},EnsGeneNames,'exact'))
                                    %one of the Ensembl genes overlapped by the AceView gene has been found
                                    %so do not consider the other as true target
                                    AddFlag=0;
                                    break
                                end
                            end
                            if AddFlag
                                for GeneL2=1:length(EnsGenesDict{AcePos})
                                    %if the gene is not in the current list of targeted genes add it
                                    if isempty(strmatch(EnsGenesDict{AcePos}{GeneL2},EnsGeneNames,'exact'))
                                        EnsGeneNames=[EnsGeneNames,{EnsGenesDict{AcePos}{GeneL2}}];
                                        %consider that the nb of probe is the same
                                        EnsProbeNb=[EnsProbeNb,AceProbeNb(GeneL1)];
                                    end
                                end
                            end
                        end

                    else
                        %recover gene names specific to AceView
                        CurrAceGeneNames=[CurrAceGeneNames,{AceGeneNames{GeneL1}}];
                        CurrAceProbeNb=[CurrAceProbeNb,AceProbeNb(GeneL1)];
                    end
                end
            end
        end
    end
    [EnsGeneNames,Pos1,tmp]=unique(EnsGeneNames);
    EnsProbeNb=EnsProbeNb(Pos1);
    [EnsGeneNames,SortIndex]=sort(EnsGeneNames);
    EnsProbeNb=EnsProbeNb(SortIndex);
    if AceFlag
        [CurrAceGeneNames,Pos1,tmp]=unique(CurrAceGeneNames);
        CurrAceProbeNb=CurrAceProbeNb(Pos1);
        [CurrAceGeneNames,SortIndex]=sort(CurrAceGeneNames);
        CurrAceProbeNb=CurrAceProbeNb(SortIndex);
        NewPs{PsL}.geneNames=[EnsGeneNames,CurrAceGeneNames];
        NewPs{PsL}.source=[repmat(1,1,length(EnsGeneNames)),repmat(2,1,length(CurrAceGeneNames))];
    else
        NewPs{PsL}.geneNames=EnsGeneNames;
        NewPs{PsL}.source=[repmat(1,1,length(EnsGeneNames))];
    end
    NewPs{PsL}.target=repmat(1,1,length(NewPs{PsL}.geneNames));


    %correct source for GOP : the corresponding gene if it exists
    %is not referenced in the used database (Ensembl, AceView...)
    GopPos=strmatch('GOP',EnsGeneNames);
    NewPs{PsL}.source(GopPos)=0;
    NewPs{PsL}.target(GopPos)=0;
    if AceFlag
        NewPs{PsL}.probeNb=[EnsProbeNb,CurrAceProbeNb];
    else
        NewPs{PsL}.probeNb=EnsProbeNb;
    end


    %fill probesets by ensembl gene
    if ~isempty(EnsGeneNames)
        for GeneL=1:length(EnsGeneNames)
            Pos=strmatch(EnsGeneNames{GeneL},Genes.name,'exact');
            if isempty(Pos)
                Genes.name{end+1,1}=EnsGeneNames{GeneL};
                %Genes.source(end+1,1)=1;
                Pos=length(Genes.name);
                PsBy.gene{Pos,1}=PsL;
            else
                try
                    PsBy.gene{Pos,1}=unique([PsBy.gene{Pos,1},PsL]);
                catch
                    % only created with out
                    PsBy.gene{Pos,1}=PsL;
                end
            end
%             if length(EnsGeneNames{GeneL})>=3
%                 if isequal('GOP',EnsGeneNames{GeneL}(1:3))
%                     Genes.source=0;
%                 end
%             end
        end
    end


    if AceFlag
    %fill probesets by aceview gene (found only in genes specific to AceView)
    if ~isempty(CurrAceGeneNames)
        for GeneL=1:length(CurrAceGeneNames)
            Pos=strmatch(CurrAceGeneNames{GeneL},Genes.name,'exact');
            if isempty(Pos)
                Genes.name{end+1,1}=CurrAceGeneNames{GeneL};
%                 Genes.source(end+1,1)=2;
                Pos=length(Genes.name);
                PsBy.gene{Pos,1}=PsL;
            else
                try
                    PsBy.gene{Pos,1}=unique([PsBy.gene{Pos,1},PsL]);
                catch
                    % only created with out
                    PsBy.gene{Pos,1}=PsL;
                end
            end
        end
    end
    end

    %recover genes targetted by the current probeset inside up, down and introns
    EnsGeneNames=Ps{1}{PsL}.geneNamesOutSup;
    EnsProbeNb=Ps{1}{PsL}.probeNbOutSup;
    if AceFlag
        AceGeneNames=Ps{2}{PsL}.geneNamesOutSup;
        AceProbeNb=Ps{2}{PsL}.probeNbOutSup;
        CurrAceGeneNames={};
        CurrAceProbeNb=[];
        %search for gene specific to AceView
        if ~isempty(AceGeneNames)
            %process each AceView gene targeted by the current probe set
            for GeneL1=1:length(AceGeneNames)
                %process genes that have AceView names
                if isempty(findstr('ENS',AceGeneNames{GeneL1}))
                    AcePos=strmatch(AceGeneNames{GeneL1},AceGenesDict,'exact');
                    if ~isempty(EnsGenesDict{AcePos})
                        % if length=1, it is a one to one
                        % correspondance between AceVeiw and ensembl
                        % genes
                        if length(EnsGenesDict{AcePos})==1
                            %if the gene is not in the current list of targeted genes add it
                            if isempty(strmatch(EnsGenesDict{AcePos}{1},EnsGeneNames,'exact'))
                                EnsGeneNames=[EnsGeneNames,{EnsGenesDict{AcePos}{1}}];
                                EnsProbeNb=[EnsProbeNb,AceProbeNb(GeneL1)];
                            end
                        else
                            %if length >1 it is a composit gene (an AceView gene that overlap two ensembl gene
                            %add an Ensemlb gene only if none are in common
                            %otherwise would systematically add genes close to true Ensembl targets
                            AddFlag=1;
                            for GeneL2=1:length(EnsGenesDict{AcePos})
                                if ~isempty(strmatch(EnsGenesDict{AcePos}{GeneL2},EnsGeneNames,'exact'))
                                    %one of the Ensembl genes overlapped by the AceView gene has been found
                                    %so do not consider the other as true target
                                    AddFlag=0;
                                    break
                                end
                            end
                            if AddFlag
                                for GeneL2=1:length(EnsGenesDict{AcePos})
                                    %if the gene is not in the current list of targeted genes add it
                                    if isempty(strmatch(EnsGenesDict{AcePos}{GeneL2},EnsGeneNames,'exact'))
                                        EnsGeneNames=[EnsGeneNames,{EnsGenesDict{AcePos}{GeneL2}}];
                                        %consider that the nb of probe is the same
                                        EnsProbeNb=[EnsProbeNb,AceProbeNb(GeneL1)];
                                    end
                                end
                            end
                        end

                    else
                        %recover gene names specific to AceView
                        CurrAceGeneNames=[CurrAceGeneNames,{AceGeneNames{GeneL1}}];
                        CurrAceProbeNb=[CurrAceProbeNb,AceProbeNb(GeneL1)];
                    end
                end
            end
        end
    end

    if ~isempty(EnsGeneNames)
        [EnsGeneNames,Pos1,tmp]=unique(EnsGeneNames);
        EnsProbeNb=EnsProbeNb(Pos1);
        [EnsGeneNames,SortIndex]=sort(EnsGeneNames);
        EnsProbeNb=EnsProbeNb(SortIndex);


        % remove OutGenes if already exists as InGenes
        ClearIndex=[];
        for GeneL=1:length(EnsGeneNames)
            if ~isempty(strmatch(EnsGeneNames{GeneL},NewPs{PsL}.geneNames,'exact'))
                ClearIndex=[ClearIndex,GeneL];
            end
        end
        EnsGeneNames(ClearIndex)=[];
        EnsProbeNb(ClearIndex)=[];
    end
    if AceFlag
        if ~isempty(CurrAceGeneNames)
            [CurrAceGeneNames,Pos1,tmp]=unique(CurrAceGeneNames);
            CurrAceProbeNb=CurrAceProbeNb(Pos1);
            [CurrAceGeneNames,SortIndex]=sort(CurrAceGeneNames);
            CurrAceProbeNb=CurrAceProbeNb(SortIndex);
            % remove OutGenes if already exists as InGenes
            ClearIndex=[];
            for GeneL=1:length(CurrAceGeneNames)
                if ~isempty(strmatch(CurrAceGeneNames{GeneL},NewPs{PsL}.geneNames,'exact'))
                    ClearIndex=[ClearIndex,GeneL];
                end
            end
            CurrAceGeneNames(ClearIndex)=[];
            CurrAceProbeNb(ClearIndex)=[];
        end
    else
        CurrAceGeneNames=[];
    end
    if length([EnsGeneNames,CurrAceGeneNames])>0
        if length(NewPs{PsL}.geneNames)>0
            NewPs{PsL}.geneNames=[NewPs{PsL}.geneNames,[EnsGeneNames,CurrAceGeneNames]];
            NewPs{PsL}.source=[NewPs{PsL}.source,repmat(1,1,length(EnsGeneNames)),repmat(2,1,length(CurrAceGeneNames))];
            NewPs{PsL}.target=[NewPs{PsL}.target,repmat(2,1,length([EnsGeneNames,CurrAceGeneNames]))];
            %correct source for GOP : the corresponding gene if it exists
            %is not referenced in the used database (Ensembl, AceView...)
            GopPos=strmatch('GOP',NewPs{PsL}.geneNames);
            NewPs{PsL}.source(GopPos)=0;
            NewPs{PsL}.target(GopPos)=0;
            if AceFlag
                NewPs{PsL}.probeNb=[NewPs{PsL}.probeNb,EnsProbeNb,CurrAceProbeNb];
            else
                NewPs{PsL}.probeNb=[NewPs{PsL}.probeNb,EnsProbeNb];
            end
        else
            NewPs{PsL}.geneNames=[EnsGeneNames,CurrAceGeneNames];
            NewPs{PsL}.source=[repmat(1,1,length(EnsGeneNames)),repmat(2,1,length(CurrAceGeneNames))];
            NewPs{PsL}.target=repmat(2,1,length([EnsGeneNames,CurrAceGeneNames]));
            GopPos=strmatch('GOP',EnsGeneNames);
            NewPs{PsL}.source(GopPos)=0;
            NewPs{PsL}.target(GopPos)=0;
            if AceFlag
                NewPs{PsL}.probeNb=[EnsProbeNb,CurrAceProbeNb];
            else
                NewPs{PsL}.probeNb=[EnsProbeNb];
            end
        end
    end

    %fill probesets by ensembl gene
    if ~isempty(EnsGeneNames)
        for GeneL=1:length(EnsGeneNames)
            Pos=strmatch(EnsGeneNames{GeneL},Genes.name,'exact');
            if isempty(Pos)
                Genes.name{end+1,1}=EnsGeneNames{GeneL};
%                 Genes.source(end+1,1)=1;
                Pos=length(Genes.name);
                PsBy.gene{Pos,1}=PsL;
            else
                try
                    PsBy.gene{Pos,1}=unique([PsBy.gene{Pos,1},PsL]);
                catch %not yet Ps ranks for Out
                    PsBy.gene{Pos,1}=PsL;
                end
            end
%             if length(EnsGeneNames{GeneL})>=3
%                 if isequal('GOP',EnsGeneNames{GeneL}(1:3))
%                     Genes.source=0;
%                 end
%             end
        end
    end

    if AceFlag
        %fill probesets by aceview gene (found only in genes specific to AceView)
        if ~isempty(CurrAceGeneNames)
            for GeneL=1:length(CurrAceGeneNames)
                Pos=strmatch(CurrAceGeneNames{GeneL},Genes.name,'exact');
                if isempty(Pos)
                    Genes.name{end+1,1}=CurrAceGeneNames{GeneL};
%                    Genes.source(end+1,1)=2;
                    Pos=length(Genes.name);
                    PsBy.gene{Pos,1}=PsL;
                else
                    try
                        PsBy.gene{Pos,1}=unique([PsBy.gene{Pos,1},PsL]);
                    catch  %not yet Ps ranks for Out
                        PsBy.gene{Pos,1}=PsL;
                    end
                end
            end
        end
    end
end

%% FILL_PAIRED
%file PsRanks, Pv and Corr
function [NewPs,NotFound]=FILL_PAIRED(Limit,NewPs,Sim,AllVal,PsBy,Genes,RoundRank)

PsNb=length(NewPs);
NetNb=size(AllVal,2);
if RoundRank==1
    for PsL=1:PsNb
        NewPs{PsL,1}.geneRanks={};
        NewPs{PsL,1}.psRanks={};
        NewPs{PsL,1}.paired={};
        NewPs{PsL,1}.probeNbs={};   
        NewPs{PsL,1}.sources={};
        NewPs{PsL,1}.targets={};
    end
end
NotFound=[];
GoodPosNb=0;
%process each probe set
for PsL=1:PsNb
    GeneList=Genes.name;
    %Genes targeted by the current probe set
    GeneNames=NewPs{PsL}.geneNames;
    %CurrPsBy=PsBy.gene;
    NewPsRanks='NewPs{PsL}.psRanks{GeneL}';
    NewPsPaired='NewPs{PsL}.paired{GeneL}';
    if ~isempty(GeneNames)
        %process each gene
        for GeneL=1:length(GeneNames)
            %recover the rank of probe sets that target the same gene
            Pos=strmatch(GeneNames{GeneL},GeneList,'exact');
            %PsRanks=CurrPsBy{Pos};
            PsRanks=PsBy.gene{Pos};
            if length(PsRanks)>1
                PsPos=find(PsRanks==PsL);
                if isempty(PsPos)
                    h=errordlg(sprintf('%u not in PsRanks',PsL));
                    waitfor(h)
                    error('process canceled')
                end
                %eliminate current probe set from the list
                PsRanks(PsPos)=[];
                PsRanks=unique(PsRanks);
                if RoundRank==1                    
                    %add PsL
                    eval(sprintf('%s=sort([PsL,PsRanks]);',NewPsRanks))
                    Paired=cell(1,length(PsRanks));
                else
                    %recover existing values
                    eval(sprintf('Paired=%s;',NewPsPaired))
                end
                %process all the probe set
                for PsL1=1:length(PsRanks)
                    Pos=find(Sim.firstPsRank==PsL &Sim.sndPsRank==PsRanks(PsL1));
                    if isempty(Pos)
                        Pos=find(Sim.sndPsRank==PsL &Sim.firstPsRank==PsRanks(PsL1));
                    end
                    if isempty(Pos)
                        if RoundRank==1
                            NotFound=[NotFound;sort([PsL,PsRanks(PsL1)])];
                            Paired{PsL1}=ones(1,NetNb)*NaN;
                        end
                    else
                        Paired{PsL1}=uint8(AllVal.corr(Pos(1),:)>0|AllVal.anti(Pos(1),:)>0);
                        GoodPos=AllVal.corr(Pos(1),:)>=Limit.corr&AllVal.anti(Pos(1),:)<=Limit.anti&AllVal.pv(Pos(1),:)<=Limit.pv;
                        if ~isempty(find(GoodPos))
                            Paired{PsL1}(GoodPos)=2;
                            GoodPosNb=GoodPosNb+1;
                        end
                    end
                end
                eval(sprintf('%s=Paired;',NewPsPaired))
            else
                %only one probe set (the current one, PsL) targets the current gene
                if RoundRank==1
                    eval(sprintf('%s=PsRanks;',NewPsRanks))
                    eval(sprintf('%s={};',NewPsPaired))
                else
                    if PsRanks(1)~=PsL
                        errordlg(sprintf('Probe set %u : Gene %s was already targeted by Probe set %u',PsL,GeneNames{GeneL},PsRanks(1)))
                        error('process canceled')
                    end
                end
            end
        end
    end
end
if ~isempty(NotFound)
    NotFound=unique(NotFound,'rows');
    [Temp,SortIndex]=sort(NotFound(:,1));
    NotFound=NotFound(SortIndex,:);
    %transform into cell
    NotFound=mat2cell(NotFound,ones(size(NotFound,1),1),2);
end
if RoundRank==1
    %ADD geneRanks,probeNbs, sources and targets field
    for PsL1=1:PsNb
        if ~isempty(NewPs{PsL1}.geneNames)
            GeneNb=length(NewPs{PsL1}.geneNames);
            NewPs{PsL1}.geneRanks=zeros(1,GeneNb);
            NewPs{PsL1}.probeNbs=cell(1,GeneNb);
            NewPs{PsL1}.sources=cell(1,GeneNb);
            NewPs{PsL1}.targets=cell(1,GeneNb);
            for GeneL=1:GeneNb               
                NewPs{PsL1}.geneRanks(GeneL)=strmatch(NewPs{PsL1}.geneNames{GeneL},Genes.name,'exact');
                CurrPsNb=length(NewPs{PsL1}.psRanks{GeneL});
                NewPs{PsL1}.probeNbs{GeneL}=zeros(1,CurrPsNb);
                NewPs{PsL1}.sources{GeneL}=zeros(1,CurrPsNb);
                NewPs{PsL1}.targets{GeneL}=zeros(1,CurrPsNb);
                for PsL2=1:CurrPsNb
                    CurrPsRank2=NewPs{PsL1}.psRanks{GeneL}(PsL2);
                    if CurrPsRank2==PsL1
                        NewPs{PsL1}.probeNbs{GeneL}(PsL2)=NewPs{PsL1}.probeNb(GeneL);
                        NewPs{PsL1}.sources{GeneL}(PsL2)=NewPs{PsL1}.source(GeneL);
                        NewPs{PsL1}.targets{GeneL}(PsL2)=NewPs{PsL1}.target(GeneL);
                    else
                        GenePos2=strmatch(NewPs{PsL1}.geneNames{GeneL},NewPs{CurrPsRank2}.geneNames,'exact');
                        NewPs{PsL1}.probeNbs{GeneL}(PsL2)=NewPs{CurrPsRank2}.probeNb(GenePos2);
                        NewPs{PsL1}.sources{GeneL}(PsL2)=NewPs{CurrPsRank2}.source(GenePos2);
                        NewPs{PsL1}.targets{GeneL}(PsL2)=NewPs{CurrPsRank2}.target(GenePos2);
                    end
                end
            end
        end
    end
end


%% CONSTRUCT_PSBY
function PsBy=CONSTRUCT_PSBY(NewPs,PsBy,Genes)
%PsBy is initiatied in CONSTRUCT_NEWPS
PsNb=length(NewPs);
PsBy.paired=cell(length(PsBy.gene),1);
for PsL=1:PsNb
    GeneList=Genes.name;
    GeneNames=NewPs{PsL}.geneNames;
    Paired=NewPs{PsL}.paired;
    PsRanks=NewPs{PsL}.psRanks;

    PairedMat=PsBy.paired;

    if ~isempty(GeneNames)
        if ~isempty(PsRanks)
            for GeneL=1:length(GeneNames)
                CurrPsRanks=PsRanks{GeneL};
                if length(CurrPsRanks)>1
                    CurrPaired=Paired{GeneL};
                    %sort all on PsRanks order
                    CurrPsRanks=sort(CurrPsRanks);
                    PsNb=length(CurrPsRanks);
                    %position of the current ps in the matrixes
                    PsPos=find(CurrPsRanks==PsL);
                    if isempty(PsPos)
                        h=errordlg(sprintf('PsL %u absent in PsRanks{gene %u}',PsL,GeneL));
                        waitfor(h)
                        error('process canceled')
                    end
                    %create or load current matrixes
                    GenePos=strmatch(GeneNames{GeneL},GeneList,'exact');
                    if isempty(PairedMat{GenePos})
                        CurrPairedMat=cell(PsNb*(PsNb-1)/2,1);
                    else
                        CurrPairedMat=PairedMat{GenePos};
                    end

                    %insert the CurrPaired values at the rigth position in
                    %in the vectorized upper triangle of the matrix
                    %PsNb*PsNb, respecting i<j
                    if PsPos==1
                        for j=1:PsNb-1
                            CurrPairedMat{j}=CurrPaired{j};
                        end
                    else
                        if PsPos==PsNb
                            j=PsNb;
                            for i=1:PsNb-1
                                %position in the vectorized upper triangle
                                %of the matrix PsNb*PsNb
                                CurrPos=j-PsNb+(i*(2*PsNb-i-1))/2;
                                CurrPairedMat{CurrPos}=CurrPaired{i};
                            end
                        else
                            j=PsPos;
                            for i=1:j-1
                                CurrPos=j-PsNb+(i*(2*PsNb-i-1))/2;
                                CurrPairedMat{CurrPos}=CurrPaired{i};
                            end
                            i=PsPos;
                            for j=PsPos+1:PsNb
                                CurrPos=j+-PsNb+(i*(2*PsNb-i-1))/2;
                                CurrPairedMat{CurrPos}=CurrPaired{j-1};
                            end
                        end
                    end
                    %update PVMat and PairedMat
                    PairedMat{GenePos}=CurrPairedMat;
                end
            end
            PsBy.paired=PairedMat;
        end
    end
end

%% STAT
function [NewPs,Stat,LinkedPs]=STAT(Species,ChipRank,NewPs,PsBy,NetNb,ProbeNb,TestLimit,SumFlag,DisplayFlag)
global K
LinkedPs={};
PsNb=length(NewPs);
for PsL=1:PsNb
    NewPs{PsL,1}.grp=cell(1,length(NewPs{PsL}.geneNames));
end
%The number of genes that are targetted only by the current probe set.
Stat.singleTargNb=zeros(PsNb,1);
%The number of genes that are targetted by the current probe set plus a
%single other one.
Stat.doubleTargNb=zeros(PsNb,1);
%The number of genes that are are targetted by the current probe set plus
%two or more another probe sets.
Stat.multipleTargNb=zeros(PsNb,1);
%The maximum nb of probe sets targetting a gene targeted by the current
%probe set.
Stat.maxPsNb=zeros(PsNb,1);
%The distribution of ps group sizes
Stat.grpSizes=[];
%The type of link between to paired probe sets (0: no corr, 1: corr>0 but
%not similar (don't pass the test), 2: similar (pass the test))
Stat.linkTypes=[];
%properties of pairs that don't pass the test but are in a ps group
%[PsL, GeneL, Node1, Node2, Nb of bad links in the group, Total nb of links in 
%the group, Nb pairs with corr>0 (PairedMat=1), Nb of similar pairs (pairedMat=2)];
Stat.badLinks=[];
%During merging process of triangle, some transitory groups of probe sets
%are split into two new  groups
%single probe sets that are common to these two groups are taken away and
%considered as forming a hubs which is in relation with the two groups.
%during merging process of triangle, some transitory groups of probe sets
%are split into two new  groups.
%single probe sets that are common to these two groups are taken away and
%considered as forming a special group called a hub or a pivot which is in 
%relation with the two groups.
%[hub rank,first group rank,second group rank,size of the first group,size of the snd group]
Stat.hubs=[];
%Information on pairs of probe set that belong to the same group:
%LinkedPs{1}: they target a gene with a number of probe smaller than the number of probes
%             of the current ps targeting the assigned gene
%LinkedPs{2}: they target a gene with a number of probe equal to the number of probes of
%             the current ps targeting the assigned gene
%[Ps1,Ps2,Target type of targeted gene, Nb of probes of Ps1 targeting the gene,...
%         Nb of probes of Ps2 targeting the gene, Source type of the targeted gene,
%         maximum nb of probes of Ps1 that target a gene,
%         maximum nb of probes of Ps1 that target a gene];
%LinkedPs{1}(PairPos,:)=[Ps1,Ps2,NewPs{Ps1}.target(GenePos1),NewPs{Ps1}.probeNb(GenePos1),
%         NewPs{Ps2}.probeNb(GenePos2), NewPs{Ps1}.source(GenePos1),max(NewPs{Ps1}.probeNb),
%         max(NewPs{Ps2}.probeNb)];
GrpSizeNb=1;
if SumFlag
    LinkType=1;
else
    LinkType=2;
end
%process each probe set
for PsL=1:PsNb
    if ~isempty(NewPs{PsL}.geneNames)
        %recover the maximum number of probe sets that can target a gene
        MaxPs=0;
        %process each targeted gene
        for GeneL=1:length(NewPs{PsL}.geneNames)            
            GeneRank=NewPs{PsL}.geneRanks(GeneL);            
            %the current gene is targeted only by the current probe set
            if length(NewPs{PsL}.psRanks{GeneL})==1
                if NewPs{PsL}.psRanks{GeneL}~=PsL
                    h=errordlg('error ~=PsL');
                    waitfor(h)
                    error('cancel process')
                else
                    MaxPs=max(MaxPs,1);
                    NewPs{PsL}.grp{GeneL}{1}=PsL;
                    Stat.singleTargNb(PsL)=Stat.singleTargNb(PsL)+1;
                end
            %the current gene is targeted by the current probe set and one another probe set    
            elseif length(NewPs{PsL}.psRanks{GeneL})==2
                MaxPs=max(MaxPs,2);
                %test if the pair of probe sets targeting the current gene constitutes 
                %a group
                if  length(find((NewPs{PsL}.paired{GeneL}{1}>=LinkType)))>=TestLimit
                    NewPs{PsL}.grp{GeneL}{1}=NewPs{PsL}.psRanks{GeneL};
                else
                    NewPs{PsL}.grp{GeneL}{1}=NewPs{PsL}.psRanks{GeneL}(1);
                    NewPs{PsL}.grp{GeneL}{2}=NewPs{PsL}.psRanks{GeneL}(2);
                end
                Stat.doubleTargNb(PsL)=Stat.doubleTargNb(PsL)+1;
            %the current gene is targeted by the curent probe sets and at least two another
            %probe sets
            else               
                MaxPs=max(MaxPs,length(NewPs{PsL}.psRanks{GeneL}));                
                Stat.multipleTargNb(PsL)=Stat.multipleTargNb(PsL)+1;
                if isempty(NewPs{PsL}.grp{GeneL})
                    PsRanks=NewPs{PsL}.psRanks{GeneL};
                    CurrGeneRank=NewPs{PsL}.geneRanks(GeneL);
                    [NewPs{PsL}.grp{GeneL},GrpSizes,LinkTypes,BadLinks,Hubs]=make_psgroups(PsRanks,PsBy.paired{GeneRank},NetNb,TestLimit,GrpSizeNb,LinkType);
                    if ~isempty(LinkTypes)
                        Stat.linkTypes=[Stat.linkTypes;[repmat(PsL,size(LinkTypes,1),1),repmat(GeneL,size(LinkTypes,1),1),LinkTypes]];
                    end
                    if ~isempty(BadLinks)
                        Stat.badLinks=[Stat.badLinks;[repmat(PsL,size(BadLinks,1),1),repmat(GeneL,size(BadLinks,1),1),BadLinks]];
                    end
                    if length(GrpSizes)+2> size(Stat.grpSizes,2)
                        GrpSizeNb=length(GrpSizes);
                        for i=size(Stat.grpSizes,2)+1:length(GrpSizes)+2
                            Stat.grpSizes=[Stat.grpSizes,zeros(size(Stat.grpSizes,1),1)];
                        end
                    end                   
                    Stat.grpSizes=[Stat.grpSizes;[repmat(PsL,size(GrpSizes,1),1),repmat(GeneL,size(GrpSizes,1),1),GrpSizes]];
                    %Hubs :[PivotRank,FirstGrp,SndGrp,length(Common),length(FirstGrp),length(SndGrp)]
                    %Stat.hubs:[PivotRank,GeneL,PivotPos,FirstGrp,SndGrp,length(Common),length(FirstGr),length(SndGrp)]
                    if ~isempty(Hubs)
                        for HubL=1:size(Hubs,1)
                            %Pivot is a single probe set
                            PivotRank=NewPs{PsL}.grp{GeneL}{Hubs(HubL,1)};
                            CurrGeneL=find(NewPs{PivotRank}.geneRanks==CurrGeneRank);
                            Stat.hubs=[Stat.hubs;[PivotRank,CurrGeneL,Hubs(HubL,:)]];
                        end
                    end                    
                    %Fill otherPs
                    PsRanks=setdiff(PsRanks,PsL);
                    for PsL1=1:length(PsRanks)
                        CurrPsRank=PsRanks(PsL1);
                        GenePos=find(NewPs{CurrPsRank}.geneRanks==GeneRank);
                        NewPs{CurrPsRank}.grp{GenePos}=NewPs{PsL}.grp{GeneL};
                    end                    
                end
            end
        end
        Stat.maxPsNb(PsL)=MaxPs;
    end
end

%% -- CATEGORIES STAT
%s=single
%m=multiple
%first position : ps
%second position : gene
Stat.names={'su';'sm';'mu';'mm';'cx';'hx';'nogene'};
Stat.index=cell(length(Stat.names),1);
Bindex=ones(PsNb,1);

%FIRST CATEGORY
%a single Ps targeting a unique gene
%'su'
a=find(Stat.singleTargNb==1 & Stat.doubleTargNb==0 & Stat.multipleTargNb==0);
if ~isempty(find(Bindex(a)==0)),sprintf '1',end
Stat.index{1}=a;
Bindex(a)=0;
Stat.catcount(1)=length(Stat.index{1});

%SECOND CATEGORY
%a single Ps targeting multiple genes
%'sm'
a=find(Stat.singleTargNb>1 & Stat.doubleTargNb==0 & Stat.multipleTargNb==0);
if ~isempty(find(Bindex(a)==0)),sprintf '2',end
Stat.index{2}=a;
Bindex(a)=0;
Stat.catcount(2)=length(Stat.index{2});

%THIRD CATEGORY
%several Ps targeting  a unique gene
%'mu'
a=find(Stat.singleTargNb==0 & Stat.doubleTargNb+Stat.multipleTargNb==1);
if ~isempty(find(Bindex(a)==0)),sprintf '3',end
Stat.index{3}=a;
Bindex(a)=0;
Stat.catcount(3)=length(Stat.index{3});

%FOURTH CATEGORY
%several Ps targeting several genes
%'mm'
a=find((Stat.singleTargNb==0 & Stat.doubleTargNb+Stat.multipleTargNb>1)|(Stat.singleTargNb>0 & Stat.doubleTargNb+Stat.multipleTargNb>0));
if ~isempty(find(Bindex(a)==0)),sprintf '4',end
Stat.index{4}=a;
Bindex(a)=0;
Stat.catcount(4)=length(Stat.index{4});

%SEVENTH CATEGORY
%'nogene'
Stat.index{7}=find(Bindex==1);
Stat.catcount(7)=length(Stat.index{7});

%CHECK CATEGORIES
Index=Stat.index{7};
Recover=[];
for i=1:length(Index)
    if ~isempty(NewPs{Index(i)}.geneNames)
        Recover=[Recover;Index(i)];
    end
end
if ~isempty(Recover)
    h=errordlg('no gene is false, exist probe set with gene in this category');
    waitfor(h)
    error('canceled')
end
if sum(Stat.catcount)~=PsNb
    h=errordlg(sprintf('%u probe set in categories instead of %u',sum(Stat.catcount),PsNb));
    waitfor(h)
    error('canceled')
end
for CatL1=1:6
    for CatL2=CatL1+1:7
        if length(intersect(Stat.index{CatL1},Stat.index{CatL2}))>0
            h=warndlg(sprintf('%uinter%u: %u',CatL1,CatL2,length(intersect(Stat.index{CatL1},Stat.index{CatL2}))));
            waitfor(h)
        end
    end
end


%% -- TRUE AND COMPLEX MS & MM
%add categories
MemStat=Stat;
%CHECK MS CATEGORY
MuPsNb=MemStat.catcount(3);
MuPsRanks=MemStat.index{3};
TrueMuPsRanks=[];
%prevent from doing the same check several times
MuPsMade=zeros(MuPsNb,1);
MmPsNb=MemStat.catcount(4);
MmPsRanks=MemStat.index{4};
%CHECK MS CATEGORY
for PsL1=1:MuPsNb
    if MuPsMade(PsL1)==0
        PsRank1=MuPsRanks(PsL1);
        %check that there exist only one target
        if length(NewPs{PsRank1}.geneNames)>1
            h=warndlg(sprintf('%u not in MS class',PsRank1));
            waitfor(h)
        else
            %if one of the ps ranks that targets the current unique gene, targets
            % also other genes, it belong to mm class, and the current probe set
            % and all its associated probe sets are not a true MS probe set
            %and are put into ComplexMuPsRanks
            PsRanks2=NewPs{PsRank1}.psRanks{1};
            TrueFlag=1;
            for PsL2=1:length(PsRanks2)
                PsRank2=PsRanks2(PsL2);
                Pos2=find(MuPsRanks==PsRank2);
                if isempty(Pos2)
                    %search it in MM class
                    Pos2=find(MmPsRanks==PsRank2);
                    if isempty(Pos2)
                        h=warndlg(sprintf('%u not in MM class',PsRank2));
                        waitfor(h)
                    end
                    TrueFlag=0;
                else
                    MuPsMade(Pos2)=1;
                end
            end
            if TrueFlag
                for PsL2=1:length(PsRanks2)
                    TrueMuPsRanks=[TrueMuPsRanks;PsRanks2(PsL2)];
                end
            end
            MuPsMade(PsL1)=1;
        end
    end
end
TrueMuPsRanks=unique(TrueMuPsRanks);
ComplexMuPsRanks=setdiff(MuPsRanks,TrueMuPsRanks);

%CHECK MM CATEGORY
MmPsMade=zeros(MmPsNb,1);
TrueMmPsRanks=[];
AssociatedPsRanks={};
for PsL1=1:MmPsNb
    if MmPsMade(PsL1)==0
        PsRank1=MmPsRanks(PsL1);
        GeneNames1=NewPs{PsRank1}.geneNames;
        %check that the current ps targets several genes
        if length(GeneNames1)<2
            h=warndlg(sprintf('%u not in MM class',PsRank1));
            waitfor(h)
        else
            %recover list of all gene names and all ps ranks at depth one
            GeneNames2=GeneNames1;
            PsRanks=PsRank1;
            for GeneL1=1:length(GeneNames1)
                PsRanks2=NewPs{PsRank1}.psRanks{GeneL1};
                PsRanks=union(PsRanks,PsRanks2);
                for PsL2=1:length(PsRanks2)
                    PsRank2=PsRanks2(PsL2);
                    GeneNames2=union(GeneNames2,NewPs{PsRank2}.geneNames);
                    for GeneL2=1:length(NewPs{PsRank2}.geneNames)
                        PsRanks=union(PsRanks,NewPs{PsRank2}.psRanks{GeneL2});
                    end
                end
            end
            %check that all the recovered probe sets target all the recovered genes
            %if it is not the case, recovered probe sets are considered as associated
            TrueFlag=1;
            for PsL2=1:length(PsRanks)
                for GeneL2=1:length(GeneNames2)
                    if isempty(strmatch(GeneNames2{GeneL2},NewPs{PsRanks(PsL2)}.geneNames,'exact'))
                        TrueFlag=0;
                        break
                    end
                end
            end
            if TrueFlag
                for PsL2=1:length(PsRanks)                    
                    TrueMmPsRanks=[TrueMmPsRanks;PsRanks(PsL2)];
                end
            else
                AssociatedPsRanks{end+1,1}=PsRanks;
            end
            for PsL2=1:length(PsRanks)
                PsRank2=PsRanks(PsL2);
                Pos2=find(MmPsRanks==PsRank2);
                if ~isempty(Pos2) %PsRanks2 could refers to a MS class probe set
                    MmPsMade(Pos2)=1;
                end
            end
        end
    end
end
TrueMmPsRanks=unique(TrueMmPsRanks);
ComplexMmPsRanks=setdiff(MmPsRanks,TrueMmPsRanks);
%exist probe sets that have been classified as Complex MS that are in
%reality True MM
ComplexMuPsRanks=setdiff(ComplexMuPsRanks,TrueMmPsRanks);
%RECOVER COMPLEX PS
ComplexPsRanks=unique([ComplexMmPsRanks;ComplexMuPsRanks]);

Miss=setdiff(1:PsNb,unique([MemStat.index{1};MemStat.index{2};MemStat.index{7};TrueMuPsRanks;TrueMmPsRanks;ComplexPsRanks]));
if length(Miss)>0
    h=errordlg('WRONG CLASSIFICATION');
    waitfor(h)
    error('process canceled')
end

%CORRECTION OF MS AND MM CATEGORIES
Stat.index{3}=TrueMuPsRanks;
Stat.catcount(3)=length(Stat.index{3});
Stat.index{4}=TrueMmPsRanks;
Stat.catcount(4)=length(Stat.index{4});

%add PsRanks that are in ComplexPsRanks and not in AssociatedPsRanks

CxPsRanks=[];
for GrpL=1:length(AssociatedPsRanks)
    CxPsRanks=union(CxPsRanks,AssociatedPsRanks{GrpL});
end
MissCxPsRanks=setdiff(ComplexPsRanks,CxPsRanks);

MissAssociatedPsRanks={};
for  PsL=1:length(MissCxPsRanks)
    CurrPsRank=MissCxPsRanks(PsL);
    if length(NewPs{CurrPsRank}.geneNames)>1
        h=errordlg(sprintf('%u not MS',CurrPsRank));
        waitfor(h)
        error('process canceled')
    else
        MissAssociatedPsRanks{end+1,1}=NewPs{CurrPsRank}.psRanks{1};
    end
end

%find group of ps ranks that can be merged (search for cliques) =>
%Hypercomplex class

Continue=1;
while Continue
    CxNb=length(MissAssociatedPsRanks);
    Intersections=[];
    for i=1:CxNb-1
        for j=i+1:CxNb
            CurrIntersect=intersect(MissAssociatedPsRanks{i},MissAssociatedPsRanks{j});
            if ~isempty(CurrIntersect)
                Intersections(end+1,:)=[i,j,length(CurrIntersect)];
            end
        end
    end
    Intersections=unique(Intersections,'rows');
    if isempty(Intersections)
        Continue=0;
    else
        a=unique([Intersections(:,1);Intersections(:,2)]);
        Round=0;
        GrpPos={};
        InterSize=length(Intersections);
        Continue=1;
        while Continue
            Round=Round+1;
            CurrDir=pwd;
            fid=fopen('cliquer.txt','w');
            fprintf(fid,'c rank 1\n');
            fprintf(fid,'p edge %u %u\n',max(a),length(Intersections));
            for GrpL=1:size(Intersections,1)
                fprintf(fid,'e %u %u\n',Intersections(GrpL,1),Intersections(GrpL,2));
            end
            fclose(fid);            
            eval(sprintf('! %s -su %s/cliquer.txt>%s/cliquer_res.txt;',K.dir.cliquer,CurrDir,CurrDir))
            GrpPos{Round}=load('cliquer_res.txt');
            if length(GrpPos{Round})<=1
                Continue=0;
            else
                for GrpL=1:length(GrpPos{Round})
                    Pos=find(Intersections(:,1)==GrpPos{Round}(GrpL));
                    Intersections(Pos,:)=[];
                    Pos=find(Intersections(:,2)==GrpPos{Round}(GrpL));
                    Intersections(Pos,:)=[];
                end
                InterSize(end+1,1)=length(Intersections);
                if length(Intersections)==0
                    Continue=0;
                else
                    a=unique([Intersections(:,1);Intersections(:,2)]);
                end
            end
        end

        Mem=MissAssociatedPsRanks;
        NotMerged=1:size(Mem,1);
        GrpNb=length(GrpPos);
        MissAssociatedPsRanks=cell(GrpNb,1);
        for GrpL1=1:GrpNb
            for GrpL2=1:length(GrpPos{GrpL1})
                MissAssociatedPsRanks{GrpL1,1}=union(MissAssociatedPsRanks{GrpL1,1},Mem{GrpPos{GrpL1}(GrpL2)});
                Pos=find(NotMerged==GrpPos{GrpL1}(GrpL2));
                NotMerged(Pos)=[];
            end
        end
        %add not merged lists of ps ranks
        if ~isempty(NotMerged)
            for i=1:length(NotMerged)
                MissAssociatedPsRanks{end+1,1}=Mem{NotMerged(i)};
            end
        end
    end
end

AssociatedPsRanks=[AssociatedPsRanks;MissAssociatedPsRanks];

CxNb=length(AssociatedPsRanks);
Intersections=[];
for i=1:CxNb-1
    for j=i+1:CxNb
        CurrIntersect=intersect(AssociatedPsRanks{i},AssociatedPsRanks{j});
        if ~isempty(CurrIntersect)
            Intersections(end+1,:)=[i,j,length(CurrIntersect)];
        end
    end
end
Intersections=unique(Intersections,'rows');
HxPos=unique([Intersections(:,1);Intersections(:,2)]);

HxAssociatedPsRanks=AssociatedPsRanks(HxPos);
AssociatedPsRanks(HxPos)=[];

CxPsRanks=[];
for GrpL=1:length(AssociatedPsRanks)
    CxPsRanks=union(CxPsRanks,AssociatedPsRanks{GrpL});
end

ClearPos=[];
%Verify that groups of probe set that collectively target the same set of genes
GrpNb=length(AssociatedPsRanks);
for GrpL=1:GrpNb
    PsRanks=AssociatedPsRanks{GrpL};
    %verify that the PsRankList is closed
    Verif=[];
    for PsL=1:length(PsRanks)
        GeneNb=length(NewPs{PsRanks(PsL)}.geneNames);
        for GeneL=1:GeneNb
            Verif=union(Verif,NewPs{PsRanks(PsL)}.psRanks{GeneL});
        end
    end
    if ~isequal(PsRanks,Verif)
        ClearPos=[ClearPos;GrpL];
    end
end


if ~isempty(ClearPos)
    HxAssociatedPsRanks=[HxAssociatedPsRanks;AssociatedPsRanks(ClearPos)];
    AssociatedPsRanks(ClearPos)=[];
    CxPsRanks=[];
    for GrpL=1:length(AssociatedPsRanks)
        CxPsRanks=union(CxPsRanks,AssociatedPsRanks{GrpL});
    end

end


%verify that CxPsRanks is correct
Verif=[];
GrpNb=length(AssociatedPsRanks);
for GrpL=1:GrpNb
    Verif=[Verif;AssociatedPsRanks{GrpL}'];
end
Verif=unique(Verif);
Verif=sort(Verif);
if ~isequal(CxPsRanks',Verif)
    h=errordlg('differences between Verif and CxPsRanks');
    waitfor(h)
    error('process canceled')
end


% REGISTER CX (COMPLEX) CATEGORY
Stat.index{5}=CxPsRanks';
Stat.catcount(5)=length(Stat.index{5});
Stat.cxPsRanks=AssociatedPsRanks;

%ADD HX (HYPERCOMPLEX) category
HxPsRanks=[];
for GrpL=1:length(HxAssociatedPsRanks)
    HxPsRanks=union(HxPsRanks,HxAssociatedPsRanks{GrpL});
end
Stat.index{6}=HxPsRanks';
Stat.catcount(6)=length(Stat.index{6});
Stat.hxPsRanks=HxAssociatedPsRanks;

%find ps repeated in several HX groups of probe sets
Position=zeros(size(Stat.hxPsRanks,1),length(HxPsRanks));
for HxL=1:size(Stat.hxPsRanks,1)
    [Inter temp Pos]=intersect(Stat.hxPsRanks{HxL},HxPsRanks);
    Position(HxL,Pos)=1;
end
SumPos=sum(Position);
PosNb=zeros(size(Stat.hxPsRanks,1),1);
NoIntersect=zeros(size(Stat.hxPsRanks,1),1);
for HxL=1:size(Stat.hxPsRanks,1)
    PosNb(HxL)=length(find(Position(HxL,:)));
    if sum(SumPos(find(Position(HxL,:))))==PosNb(HxL)
        NoIntersect(HxL)=1;
    end
end

%% -- FIGURES STAT
% Result of triangle merging
% Statistics on added links that do not pass the test show that most of them are similar in a
% large number of networks, which justifies their recruitment in the group.
%figure of results
if DisplayFlag
    FigRank=43;
    h=figure;
    set(h,'name',sprintf('FIG%u_n%02u - m%u: PROPERTIES OF MERGED PS',FigRank,TestLimit,ChipRank))
    %nb genes targeted with 1,2 or more probe sets
    subplot(2,2,1)
    hold on
    Bin=histc(Stat.singleTargNb,[1:max(Stat.singleTargNb)]);
    plot([1:max(Stat.singleTargNb)],Bin,'g')
    %plot([1:max(Stat.singleTargNb)],Bin,'go')
    Bin=histc(Stat.doubleTargNb,[1:max(Stat.doubleTargNb)]);
    plot([1:max(Stat.doubleTargNb)],Bin,'b')
    %plot([1:max(Stat.doubleTargNb)],Bin,'bo')
    Bin=histc(Stat.multipleTargNb,[1:max(Stat.multipleTargNb)]);
    plot([1:max(Stat.multipleTargNb)],Bin,'r')
    %plot([1:max(Stat.multipleTargNb)],Bin,'ro')
    title('nb of ps targeting a gene')
    set(gca,'yscale','log')
    set(gca,'box','on')
    set(gca,'xlim',[1,20])
    xlabel('nb of genes targeted by the current ps')
    ylabel('ps fequency')
    legend({'1 ps','2 ps','>2 ps'},'location','NorthEast')


    %maximal nb of ps targeting a gene targeted by the current ps
    subplot(2,2,2)
    hold on
    Bin=histc(Stat.singleTargNb,[1:max(Stat.maxPsNb)]);
    plot([1:max(Stat.maxPsNb)],Bin,'b')
    plot([1:max(Stat.maxPsNb)],Bin,'bo')
    title({'max nb of ps targeting one of';'the gene of the current ps'})
    xlabel('ps nb')
    set(gca,'box','on')
    ylabel('ps fequency')
    set(gca,'yscale','log')

    %Sizes
    subplot(2,2,3)
    hold on
    Sizes=sum(Stat.grpSizes);
    Sizes=Sizes(3:end);
    plot(1:length(Sizes),Sizes,'g')
    plot(1:length(Sizes),Sizes,'go')
    title('distribution of group sizes')
    set(gca,'yscale','log')
    set(gca,'box','on')
    xlabel('group size')
    ylabel('group fequency')

    if ~isempty(Stat.badLinks)
        subplot(2,2,4)
        hold on
        if SumFlag
            a=histc(Stat.badLinks(:,8)+Stat.badLinks(:,7),0:NetNb);
        else
            a=histc(Stat.badLinks(:,8),0:NetNb);
        end
        plot(0:NetNb,a,'ro')
        plot(0:NetNb,a,'r')
        title('nb of significative pairs in bad links')
        set(gca,'box','on')
        %set(gca,'xlim',[1,30])
        xlabel('nb of significative pairs')
        ylabel('bad link fequency')
        %legend(Legend,'location','NorthEast')

    end
    set(gcf,'color',[1,1,1])    
    Position=get(gcf,'position');
    Position(3:4)=[800,600];
    set(gcf,'position',Position)   
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,FigRank,TestLimit),'png')
    close(h)
end


%information on pairs of probe set targeting the same group of transcripts
%but with a number of probes inferior to the number of probes of the current probe set targeting the
%assigned gene
for RoundL=1:2
    if RoundL==1
        %calculate the number of pairs
        PairNb=0;
    else
        LinkedPs{1}=zeros(PairNb,8);
        PairPos=0;
    end
    for PsL=1:length(NewPs)
        %find the gene which is targeted with the maximal number of probes of the current
        %probe set
        [MaxVal,MaxPos]=max(NewPs{PsL}.probeNb);
        MaxPos=find(NewPs{PsL}.probeNb==MaxVal);
        for GeneL=1:length(NewPs{PsL}.geneNames);
           %process genes not assigned to the current probe set
            if isempty(find(MaxPos==GeneL))
                GeneName=NewPs{PsL}.geneNames{GeneL};
                %process each grp of ps of the current gene
                for GrpL=1:length(NewPs{PsL}.grp{GeneL})
                    %process grp greater than one probe set
                    if length(NewPs{PsL}.grp{GeneL}{GrpL})>1
                        for PsL1=1:length(NewPs{PsL}.grp{GeneL}{GrpL})-1
                            for PsL2=PsL1+1:length(NewPs{PsL}.grp{GeneL}{GrpL})
                                if RoundL==1
                                     %calculate the number of pairs in the first round
                                    PairNb=PairNb+1;
                                else
                                    PairPos=PairPos+1;
                                    Ps1=NewPs{PsL}.grp{GeneL}{GrpL}(PsL1);
                                    Ps2=NewPs{PsL}.grp{GeneL}{GrpL}(PsL2);
                                    GenePos1=strmatch(GeneName,NewPs{Ps1}.geneNames,'exact');
                                    GenePos2=strmatch(GeneName,NewPs{Ps2}.geneNames,'exact');
                                    LinkedPs{1}(PairPos,:)=[Ps1,Ps2,NewPs{Ps1}.target(GenePos1),NewPs{Ps1}.probeNb(GenePos1),NewPs{Ps2}.probeNb(GenePos2),...
                                        NewPs{Ps1}.source(GenePos1),max(NewPs{Ps1}.probeNb),max(NewPs{Ps2}.probeNb)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%information on pairs of probe set targeting the same group of transcripts
%but with a number of probes equal to the number of probes of the current probe set targeting the
%assigned gene
for RoundL=1:2
    if RoundL==1
        %calculate the number of pairs
        PairNb=0;
    else
        LinkedPs{2}=zeros(PairNb,8);
        PairPos=0;
    end
    for PsL=1:length(NewPs)
        [MaxVal,MaxPos]=max(NewPs{PsL}.probeNb);
        MaxPos=find(NewPs{PsL}.probeNb==MaxVal);
        %process all the genes that have the same number of probes than the
        %gene assigned to the current probe set
        for MaxL=1:length(MaxPos)
            GeneL=MaxPos(MaxL);
            GeneName=NewPs{PsL}.geneNames{GeneL};
            for GrpL=1:length(NewPs{PsL}.grp{GeneL})
                if length(NewPs{PsL}.grp{GeneL}{GrpL})>1
                    for PsL1=1:length(NewPs{PsL}.grp{GeneL}{GrpL})-1
                        for PsL2=PsL1+1:length(NewPs{PsL}.grp{GeneL}{GrpL})
                            if RoundL==1
                                PairNb=PairNb+1;
                            else
                                PairPos=PairPos+1;
                                Ps1=NewPs{PsL}.grp{GeneL}{GrpL}(PsL1);
                                Ps2=NewPs{PsL}.grp{GeneL}{GrpL}(PsL2);
                                GenePos1=strmatch(GeneName,NewPs{Ps1}.geneNames,'exact');
                                GenePos2=strmatch(GeneName,NewPs{Ps2}.geneNames,'exact') ;
                                LinkedPs{2}(PairPos,:)=[Ps1,Ps2,NewPs{Ps1}.target(GenePos1),NewPs{Ps1}.probeNb(GenePos1),NewPs{Ps2}.probeNb(GenePos2),...
                                    NewPs{Ps1}.source(GenePos1),max(NewPs{Ps1}.probeNb),max(NewPs{Ps2}.probeNb)];
                            end
                        end
                    end
                end
            end
        end
    end
end




LinkedPs{1}=unique(LinkedPs{1},'rows');
LinkedPs{2}=unique(LinkedPs{2},'rows');

CurrProbeNb=max(max(LinkedPs{1}(:,4:5)));
PairNb(1)=size(LinkedPs{1},1);
PairNb(2)=size(LinkedPs{2},1);

if DisplayFlag
    FigRank=44;
    h=figure;
    set(h,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u_n%02u - m%u: PROPERTIES OF PAIRED PROBE SETS',FigRank,TestLimit,ChipRank))
    for SelL=1:2
        subplot(2,2,(SelL-1)*2+1)
        Bin=histc([LinkedPs{SelL}(:,4);LinkedPs{SelL}(:,5)],1:CurrProbeNb);
        plot(1:CurrProbeNb,Bin,'g')
        hold on

        Val=LinkedPs{SelL}(:,[4,5]);

        Pos=find(LinkedPs{SelL}(:,3)==1);
        Bin=histc([Val(Pos,1);Val(Pos,2)],1:CurrProbeNb);
        plot(1:CurrProbeNb,Bin,'r')

        Pos=find(LinkedPs{SelL}(:,3)==2);
        Bin=histc([Val(Pos,1);Val(Pos,2)],1:CurrProbeNb);
        plot(1:CurrProbeNb,Bin,'m')

        Pos=find(LinkedPs{SelL}(:,6)==1);
        Bin=histc([Val(Pos,1);Val(Pos,2)],1:CurrProbeNb);
        plot(1:CurrProbeNb,Bin,'b')

        Pos=find(LinkedPs{SelL}(:,6)==2);
        if ~isempty(Pos)
            AceFlag=1;
            Bin=histc([Val(Pos,1);Val(Pos,2)],1:CurrProbeNb);
            plot(1:CurrProbeNb,Bin,'c')
        else
            AceFlag=0;
        end

        if SelL==1
            title({'pairs with nb of targeting probes smaller';'than the nb of probes targeting the assigned gene'})
        else
            title({'pairs with nb of targeting probes equal to';'the nb of probes targeting the assigned gene'})
        end
        set(gca,'box','on')
        %set(gca,'xlim',[1,20])
        xlabel('probe nb targeting the common gene')
        ylabel('fequency')
        if AceFlag
            legend({'all','target=1','target=2','Ensembl','AceView'},'location','NorthEast')
        else
            legend({'all','target=1','target=2','Ensembl'},'location','NorthEast')
        end

        subplot(2,2,(SelL-1)*2+2)
        DiffProbeNb=abs(LinkedPs{SelL}(:,4)-LinkedPs{SelL}(:,5));
        MaxProbeNb=max(LinkedPs{SelL}(:,4),LinkedPs{SelL}(:,5));
        [DiffProbeNb,SortIndex]=sort(DiffProbeNb);
        MaxProbeNb=MaxProbeNb(SortIndex);
        [MaxProbeNb SortIndex]=sort(MaxProbeNb);
        DiffProbeNb=DiffProbeNb(SortIndex);
        hold on        
        plot(1:PairNb(SelL),MaxProbeNb,'k.')
        plot(1:PairNb(SelL),-DiffProbeNb,'r.')
        if SelL==1
            title('probe nb difference and max(probe nb) for each pair')
        else
            title('probe nb difference and max(probe nb) for each pair')
        end
        set(gca,'box','on')
        %set(gca,'xlim',[1,20])
        xlabel('ordered ps pairs')
        ylabel('probe nb difference & max(probe nb)')
        legend({'max','diff'},'NorthWest')        
    end
    Position=get(gcf,'position');
    Position(3:4)=[800,600];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,FigRank,TestLimit),'png')
    close(h)



    %FIGURES ON CATEGORIES

    Bindex=zeros(PsNb,1);
    su=Bindex;
    su(Stat.index{1})=1;
    sm=Bindex;
    sm(Stat.index{2})=1;
    mu=Bindex;
    mu(Stat.index{3})=1;
    mm=Bindex;
    mm(Stat.index{4})=1;
    cx=Bindex;
    cx(Stat.index{5})=1;

    %CATEGORIES MS AND MM

    for FigL=1:2
        if FigL==1
            Category=mm;
            CategoryName='MM';
        elseif FigL==2
            Category=cx;
            CategoryName='CX';
        end
        h=figure;
        set(h,'color',[1,1,1])
        set(h,'name',sprintf('FIG%u_n%02u -m%u: PROBE NB DIFFERENCE DISTRIBUTION ACCORDING TO GENE ASSIGNATION FOR %s',44+FigL,TestLimit,ChipRank,CategoryName))
        for SelL=1:2
            for DiffL=1:3
                if DiffL==1
                    subplot(2,3,(SelL-1)*3+1)
                    Diff1=find(LinkedPs{SelL}(:,4)~=LinkedPs{SelL}(:,7)&Category(LinkedPs{SelL}(:,1))&Category(LinkedPs{SelL}(:,2)));
                    Diff2=find(LinkedPs{SelL}(:,5)~=LinkedPs{SelL}(:,8)&Category(LinkedPs{SelL}(:,2)));
                    Diff=intersect(Diff1,Diff2);
                elseif DiffL==2
                    subplot(2,3,(SelL-1)*3+2)
                    Diff1=find(LinkedPs{SelL}(:,4)~=LinkedPs{SelL}(:,7)&Category(LinkedPs{SelL}(:,2)));
                    Diff2=find(LinkedPs{SelL}(:,5)~=LinkedPs{SelL}(:,8)&Category(LinkedPs{SelL}(:,2)));
                    Diff=setxor(Diff1,Diff2);
                else
                    subplot(2,3,(SelL-1)*3+3)
                    Diff1=find(LinkedPs{SelL}(:,4)==LinkedPs{SelL}(:,7)&Category(LinkedPs{SelL}(:,2)));
                    Diff2=find(LinkedPs{SelL}(:,5)==LinkedPs{SelL}(:,8)&Category(LinkedPs{SelL}(:,2)));
                    Diff=intersect(Diff1,Diff2);
                end
                DiffProbeNb=[abs(LinkedPs{SelL}(Diff,4)-LinkedPs{SelL}(Diff,5))];
                MaxProbeNb=[max(LinkedPs{SelL}(Diff,4),LinkedPs{SelL}(Diff,5))];
                [DiffProbeNb,SortIndex]=sort(DiffProbeNb);
                MaxProbeNb=MaxProbeNb(SortIndex);
                [MaxProbeNb SortIndex]=sort(MaxProbeNb);
                DiffProbeNb=DiffProbeNb(SortIndex);
                hold on
                plot(1:length(MaxProbeNb),MaxProbeNb,'k.')
                plot(1:length(MaxProbeNb),-DiffProbeNb,'r.')
                if SelL==1&DiffL==1
                    ylabel({'The gene is not assigned to';'the current probe set';'';'probe nb difference & max(probe nb)'})                    
                end
                if SelL==2&DiffL==1
                    ylabel({'The gene is assigned to';'the current probe set';'';'probe nb difference & max(probe nb)'})                    
                end
                if DiffL==1
                    title({'The gene is not the';'assigned one for both ps'})
                elseif DiffL==2
                    title({'The gene is not the';'assigned one for one ps'})
                else
                    title({'The gene is the';'assigned one for both ps'})
                end
                set(gca,'box','on')
                xlabel('ordered ps pairs')
            end
        end
        Position=get(gcf,'position');
        Position(3:4)=[800,600];
        set(gcf,'position',Position)
        cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
        saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,44+FigL,TestLimit),'png')
        close(h)
    end


    %MS CATEGORY

    FigRank=47;
    h=figure;
    set(h,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u_n%02u - m%u: PROBE NB DIFFERENCE DISTRIBUTION ACCORDING TO GENE ASSIGNATION FOR MS',FigRank,TestLimit,ChipRank))

    Diff1=find(LinkedPs{2}(:,4)==LinkedPs{2}(:,7)&mu(LinkedPs{2}(:,2)));
    Diff2=find(LinkedPs{2}(:,5)==LinkedPs{2}(:,8)&mu(LinkedPs{2}(:,2)));
    Diff=intersect(Diff1,Diff2);
    DiffProbeNb=[abs(LinkedPs{2}(Diff,4)-LinkedPs{2}(Diff,5))];
    MaxProbeNb=[max(LinkedPs{2}(Diff,4),LinkedPs{2}(Diff,5))];
    [DiffProbeNb,SortIndex]=sort(DiffProbeNb);
    MaxProbeNb=MaxProbeNb(SortIndex);
    [MaxProbeNb SortIndex]=sort(MaxProbeNb);
    DiffProbeNb=DiffProbeNb(SortIndex);
    hold on
    plot(1:length(MaxProbeNb),MaxProbeNb,'k.')
    plot(1:length(MaxProbeNb),-DiffProbeNb,'r.')
    title('The gene is the assigned one for both ps')
    set(gca,'box','on')
    %set(gca,'xlim',[1,20])
    xlabel('ordered ps pairs')
    ylabel('probe nb difference & max(probe nb)')
    set(gcf,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[400,350];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,FigRank,TestLimit),'png')
    close(h)





    FigRank=48;
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u_n%02u - m%u: PROBE NB DIFFERENCE DISTRIBUTION ACCORDING TO CATEGORY',FigRank,TestLimit,ChipRank))
    for RoundL=1:3
        if RoundL==1
            Category=find(sm);
        elseif RoundL==2
            Category=find(mm);
        else
            Category=find(cx);
        end
        FirstPNb=zeros(length(Category),1);
        SndPNb=FirstPNb;
        for PsL=1:length(Category)
            if length(NewPs{Category(PsL)}.probeNb)>1
                CurrProbeNb=NewPs{Category(PsL)}.probeNb;
                [FirstPNb(PsL),Pos]=max(CurrProbeNb);
                CurrProbeNb(Pos)=[];
                SndPNb(PsL)=max(CurrProbeNb);
            end
        end
        if RoundL>1
            %exist MS in CX and MM classes
            NullPos=find(FirstPNb==0);
            FirstPNb(NullPos)=[];
            SndPNb(NullPos)=[];
        end
        subplot(1,3,RoundL)
        [SndPNb,SortOrder]=sort(SndPNb,'descend');
        FirstPNb=FirstPNb(SortOrder);
        [FirstPNb,SortOrder]=sort(FirstPNb);
        SndPNb=SndPNb(SortOrder);
        plot(FirstPNb,'k.')
        hold on
        plot(SndPNb-FirstPNb,'r.')
        set(gca,'ylim',[-(ProbeNb-1),ProbeNb])
        if RoundL==1
            title('category SM')
        elseif RoundL==2
            title('category MM')
        else
            title('category CX')
        end
    end
    set(gcf,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[600,300];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,FigRank,TestLimit),'png')
    close(h)

    FigRank=49;
    h=figure;
    set(h,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u_n%02u - m%u: PROBE NB DIFFERENCE DISTRIBUTION ACCORDING TO TARGET TYPE',FigRank,TestLimit,ChipRank))
    for SelL=1:2
        if SelL==1
            subplot(1,2,1)
            Diff=find(LinkedPs{2}(:,3)==1);
        elseif SelL==2
            subplot(1,2,2)
            Diff=find(LinkedPs{2}(:,3)==2);
        end
        DiffProbeNb=[abs(LinkedPs{2}(Diff,4)-LinkedPs{2}(Diff,5))];
        MaxProbeNb=[max(LinkedPs{2}(Diff,4),LinkedPs{2}(Diff,5))];
        [DiffProbeNb,SortIndex]=sort(DiffProbeNb);
        MaxProbeNb=MaxProbeNb(SortIndex);
        [MaxProbeNb SortIndex]=sort(MaxProbeNb);
        DiffProbeNb=DiffProbeNb(SortIndex);
        hold on
        plot(1:length(MaxProbeNb),MaxProbeNb,'k.')
        plot(1:length(MaxProbeNb),-DiffProbeNb,'r.')
        if SelL==1
            title('the two ps are in exons or splice')
        else
            title('the two ps are outside exons or splice')
        end

        set(gca,'box','on')
        %set(gca,'xlim',[1,20])
        xlabel('ordered ps pairs')
        ylabel('probe nb difference & max(probe nb)')

    end
    set(gcf,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[600,300];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,FigRank,TestLimit),'png')
    close(h)


    FigRank=50;
    h=figure;
    set(h,'color',[1,1,1])
    set(h,'name',sprintf('FIG%u_n%02u - m%u: PROBE NB DIFFERENCE DISTRIBUTION ACCORDING TO GENE ASSIGNATION AND TARGET TYPE',FigRank,TestLimit,ChipRank))
    for SelL=1:4
        if SelL==1
            subplot(2,2,1)
            Diff=find((LinkedPs{2}(:,4)~=LinkedPs{2}(:,7)|LinkedPs{2}(:,5)~=LinkedPs{2}(:,8))&LinkedPs{2}(:,3)==1);
        elseif SelL==2
            subplot(2,2,2)
            Diff=find(LinkedPs{2}(:,4)==LinkedPs{2}(:,7)&LinkedPs{2}(:,5)==LinkedPs{2}(:,8)&LinkedPs{2}(:,3)==1);
        elseif SelL==3
            subplot(2,2,3)
            Diff=find((LinkedPs{2}(:,4)~=LinkedPs{2}(:,7)|LinkedPs{2}(:,5)~=LinkedPs{2}(:,8))&LinkedPs{2}(:,3)==2);
        else
            subplot(2,2,4)
            Diff=find(LinkedPs{2}(:,4)==LinkedPs{2}(:,7)&LinkedPs{2}(:,5)==LinkedPs{2}(:,8)&LinkedPs{2}(:,3)==2);
        end
        DiffProbeNb=[abs(LinkedPs{2}(Diff,4)-LinkedPs{2}(Diff,5))];
        MaxProbeNb=[max(LinkedPs{2}(Diff,4),LinkedPs{2}(Diff,5))];
        [DiffProbeNb,SortIndex]=sort(DiffProbeNb);
        MaxProbeNb=MaxProbeNb(SortIndex);
        [MaxProbeNb SortIndex]=sort(MaxProbeNb);
        DiffProbeNb=DiffProbeNb(SortIndex);
        hold on
        plot(1:length(MaxProbeNb),MaxProbeNb,'k.')
        plot(1:length(MaxProbeNb),-DiffProbeNb,'r.')
        if SelL==1
            title({'target=1 & the gene is not';'the assigned one for at least one ps'})
        elseif SelL==2
            title({'target=1 & the gene is';'the assigned one for both ps'})
        elseif SelL==3
            title({'target=2 & the gene is not';'the assigned one for at least one ps'})
        else
            title({'target=2 & the gene is';'the assigned one for both ps'})
        end

        set(gca,'box','on')
        %set(gca,'xlim',[1,20])
        xlabel('ordered ps pairs')
        ylabel('probe nb difference & max(probe nb)')

    end
    set(gcf,'color',[1,1,1])
    Position=get(gcf,'position');
    Position(3:4)=[600,600];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h,sprintf('m%u_fig%u_n%02u.png',ChipRank,FigRank,TestLimit),'png')
    close(h)
end





%find incomplete triangles
CurrLinkedPs{1}=unique(LinkedPs{1}(:,[1:2]),'rows');
CurrLinkedPs{2}=unique(LinkedPs{2}(:,[1:2]),'rows');
for RoundL=1:2
    if RoundL==1
        MissNb=[0,0];
        ExistNb=[0,0];
        PairNb=[0,0];
    else
        MissingPs=cell(2,1);
        MissingPs{1}=zeros(MissNb(1),4);
        MissingPs{2}=zeros(MissNb(2),4);
        ExistingPs=cell(2,1);
        ExistingPs{1}=zeros(ExistNb(1),3);
        ExistingPs{2}=zeros(ExistNb(2),3);
        PairedPs=cell(2,1);
        PairedPs{1}=zeros(PairNb(1),2);
        PairedPs{2}=zeros(PairNb(2),2);
    end
    %SelL=1 => pairs of probe set targeting the same group of transcripts
    %but with a number of probes inferior to the number of probes targeting the assigned gene
    %SelL=2 => pairs of probe set targeting the same group of transcripts
    %but with a number of probes equal to the number of probes targeting the assigned gene
    for SelL=1:2
        MissPos=0;
        ExistPos=0;
        PairPos=0;
        for LinkL=1:size(CurrLinkedPs{SelL},1)
            %Rep indicates the rang of the 'pivot' probe set in the ordered
            %set of three probe sets constituing an incomplete triangle
            %Rep=1 significates that the there exist Ps1-Ps2 & Ps1-Ps3
            %edge but not Ps2-Ps3 edge in the Ps1-Ps2-Ps3 triangle
            Pos=find(CurrLinkedPs{SelL}(:,1)==CurrLinkedPs{SelL}(LinkL,1)&CurrLinkedPs{SelL}(:,2)==CurrLinkedPs{SelL}(LinkL,2));
            Rep1=0;
            if ~isempty(Pos)
                Rep1=1;
                Ps1=CurrLinkedPs{SelL}(LinkL,1);
                Ps2=CurrLinkedPs{SelL}(LinkL,2);
                for PosL=1:length(Pos)
                    Ps3=CurrLinkedPs{SelL}(Pos(PosL),2);
                    if Ps3~=Ps2
                        if  isempty(find(CurrLinkedPs{SelL}(:,1)==Ps2&CurrLinkedPs{SelL}(:,2)==Ps3))
                            if RoundL==1
                                MissNb(SelL)=MissNb(SelL)+1;
                            else
                                MissPos=MissPos+1;
                                MissingPs{SelL}(MissPos,:)=[Ps1,Ps2,Ps3,1];
                            end
                        else
                            if RoundL==1
                                ExistNb(SelL)=ExistNb(SelL)+1;
                            else
                                ExistPos=ExistPos+1;
                                ExistingPs{SelL}(ExistPos,:)=[Ps1,Ps2,Ps3];
                            end
                        end
                    end
                end
            end

            Pos=find(CurrLinkedPs{SelL}(:,2)==CurrLinkedPs{SelL}(LinkL,2)&CurrLinkedPs{SelL}(:,1)~=CurrLinkedPs{SelL}(LinkL,1));
            Rep3=0;
            if ~isempty(Pos)
                Rep3=1;
                Ps1=CurrLinkedPs{SelL}(LinkL,1);
                Ps3=CurrLinkedPs{SelL}(LinkL,2);
                for PosL=1:length(Pos)
                    Ps2=CurrLinkedPs{SelL}(Pos(PosL),1);
                    if Ps2~=Ps1
                        if  isempty(find(CurrLinkedPs{SelL}(:,1)==Ps1&CurrLinkedPs{SelL}(:,2)==Ps2))
                            if RoundL==1
                                MissNb(SelL)=MissNb(SelL)+1;
                            else
                                MissPos=MissPos+1;
                                MissingPs{SelL}(MissPos,:)=[Ps1,Ps2,Ps3,3];
                            end
                        else
                            if RoundL==1
                                ExistNb(SelL)=ExistNb(SelL)+1;
                            else
                                ExistPos=ExistPos+1;
                                ExistingPs{SelL}(ExistPos,:)=[Ps1,Ps2,Ps3];
                            end
                        end
                    end
                end
            end

            Pos=find(CurrLinkedPs{SelL}(:,1)==CurrLinkedPs{SelL}(LinkL,2));
            Rep2=0;
            if ~isempty(Pos)
                Rep2=1;
                Ps1=CurrLinkedPs{SelL}(LinkL,1);
                Ps2=CurrLinkedPs{SelL}(LinkL,2);
                for PosL=1:length(Pos)
                    Ps3=CurrLinkedPs{SelL}(Pos(PosL),2);
                    if Ps3~=Ps1
                        if  isempty(find(CurrLinkedPs{SelL}(:,1)==Ps1&CurrLinkedPs{SelL}(:,2)==Ps3))
                            if RoundL==1
                                MissNb(SelL)=MissNb(SelL)+1;
                            else
                                MissPos=MissPos+1;
                                MissingPs{SelL}(MissPos,:)=[Ps1,Ps2,Ps3,2];
                            end
                        else
                            if RoundL==1
                                ExistNb(SelL)=ExistNb(SelL)+1;
                            else
                                ExistPos=ExistPos+1;
                                ExistingPs{SelL}(ExistPos,:)=[Ps1,Ps2,Ps3];
                            end
                        end
                    end
                end
            end
            if Rep1==0 & Rep2==0 & Rep3==0
                if RoundL==1
                    PairNb(SelL)=PairNb(SelL)+1;
                else
                    PairPos=PairPos+1;
                    Ps1=CurrLinkedPs{SelL}(LinkL,1);
                    Ps2=CurrLinkedPs{SelL}(LinkL,2);
                    PairedPs{SelL}(PairPos,:)=[Ps1,Ps2];
                end
            end
        end
    end
end


MissingPs{1}=unique(MissingPs{1},'rows');
MissingPs{2}=unique(MissingPs{2},'rows');
ExistingPs{1}=unique(ExistingPs{1},'rows');
ExistingPs{2}=unique(ExistingPs{2},'rows');
PairedPs{1}=unique(PairedPs{1},'rows');
PairedPs{2}=unique(PairedPs{2},'rows');

%% TRIANGLE_STAT
function TRIANGLE_STAT(Species,ChipRank,NewPs,PsBy,NetNb)
global K
%STAT ON TRIANGLES

%>=nb of networks where pairs pass the test
TestLimits=1:NetNb;
TestNb=length(TestLimits);
TestStat=cell(TestNb,1);
for PsL4=1:length(NewPs)
    if ~isempty(NewPs{PsL4}.geneNames)
        for GeneL=1:length(NewPs{PsL4}.geneNames)
            if length(NewPs{PsL4}.psRanks{GeneL})>2
                GenePos=NewPs{PsL4}.geneRanks(GeneL);
                PairedMat=PsBy.paired{GenePos};
                %construct a matrix with information on the number of network with
                % existing correlation between probe set which passes (2) or not (1) the test
                NbPaired{1}=zeros(1,length(PairedMat));
                NbPaired{2}=zeros(1,length(PairedMat));
                for PsL=1:size(PairedMat,1)
                    NbPaired{1}(PsL)=length(find(PairedMat{PsL}==1));
                    NbPaired{2}(PsL)=length(find(PairedMat{PsL}==2));
                end                             
                NbPaired{1}=squareform(NbPaired{1});
                NbPaired{2}=squareform(NbPaired{2});
                CurrPsNb=length(NbPaired{2});
                for TestL=1:TestNb-1
                    if ~isempty(find(NbPaired{2}>=TestLimits(TestL) & NbPaired{2}<NetNb))                    
                        %find properties of the third edge
                        for PsL1=1:CurrPsNb-1
                            for PsL2=PsL1+1:CurrPsNb                                
                                if NbPaired{2}(PsL1,PsL2)>=TestLimits(TestL)                                
                                    %search for an existing triangle
                                    if PsL2<CurrPsNb
                                        for PsL3=PsL2+1:CurrPsNb
                                            if NbPaired{2}(PsL1,PsL3)>=TestLimits(TestL)                                            
                                                TestStat{TestL}=[TestStat{TestL};[NbPaired{1}(PsL2,PsL3),NbPaired{2}(PsL2,PsL3)]];
                                            end
                                            if NbPaired{2}(PsL2,PsL3)>=TestLimits(TestL)                                            
                                                TestStat{TestL}=[TestStat{TestL};[NbPaired{1}(PsL1,PsL3),NbPaired{2}(PsL1,PsL3)]];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end                          
                if ~isempty(find(NbPaired{2}==NetNb))
                    %find properties of the third edge
                    for PsL1=1:CurrPsNb-1
                        for PsL2=PsL1+1:CurrPsNb
                            if NbPaired{2}(PsL1,PsL2)==NetNb
                                %search for an existing triangle
                                if PsL2<CurrPsNb
                                    for PsL3=PsL2+1:CurrPsNb
                                        if NbPaired{2}(PsL1,PsL3)==NetNb
                                            TestStat{TestNb}=[TestStat{TestNb};[NbPaired{1}(PsL2,PsL3),NbPaired{2}(PsL2,PsL3)]];                                        
                                        elseif NbPaired{2}(PsL2,PsL3)==NetNb
                                            TestStat{TestNb}=[TestStat{TestNb};[NbPaired{1}(PsL1,PsL3),NbPaired{2}(PsL1,PsL3)]];
                                        end                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

FigRank=41;
Legend={};
for TestL=1:TestNb
    if ~isempty(TestStat{TestL})
        Legend{1,end+1}=sprintf('delta %u (%u)',NetNb-TestLimits(TestL),length(TestStat{TestL}(:,2)));
    end
end
h=figure;
set(h,'name',sprintf('FIGU%u m%u: STAT ON THE THIRD EDGE IN G-G TRIANGLES (>=)',FigRank,ChipRank))
set(gcf,'color',[1,1,1])
Colors=colors(colormap,TestNb);
for PlotL=1:2
    subplot(1,2,PlotL)
    hold on
    for TestL=1:TestNb
        if ~isempty(TestStat{TestL})
            CurrCorrNb=histc(TestStat{TestL}(:,PlotL),0:NetNb);
            plot(0:NetNb,CurrCorrNb/sum(CurrCorrNb),'color',Colors(TestL,:))
        end
    end
    set(gca,'box','on')
    xlabel('number of positive networks')
    ylabel('relative frequency')
    if PlotL==1
        legend(Legend)
        title('nb of network where third edges do not pass the test')
    else
        title('nb of network where third edges pass the test')
    end
end
Position=get(gcf,'position');
Position(3:4)=[1000,450];
set(gcf,'position',Position)
cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
close(h)



%find index (in number of positive index) where the number of type 2 edges
%is greater thant 5% of the total number of tested edges
TestLimits=1:NetNb;
TestNb=length(TestLimits);
TestStat=cell(TestNb,1);
for PsL4=1:length(NewPs)
    if ~isempty(NewPs{PsL4}.geneNames)
        for GeneL=1:length(NewPs{PsL4}.geneNames)
            if length(NewPs{PsL4}.psRanks{GeneL})>2
                GenePos=NewPs{PsL4}.geneRanks(GeneL);
                PairedMat=PsBy.paired{GenePos};
                %construct a matrix with information on the number of network with
                % existing correlation between probe set which passes (2) or not (1) the test
                NbPaired{1}=zeros(1,length(PairedMat));
                NbPaired{2}=zeros(1,length(PairedMat));
                for PsL=1:size(PairedMat,1)
                    NbPaired{1}(PsL)=length(find(PairedMat{PsL}==1));
                    NbPaired{2}(PsL)=length(find(PairedMat{PsL}==2));
                end                             
                NbPaired{1}=squareform(NbPaired{1});
                NbPaired{2}=squareform(NbPaired{2});
                CurrPsNb=length(NbPaired{2});
                for TestL=1:TestNb-1                   
                    if ~isempty(find(NbPaired{2}==TestLimits(TestL)))
                        %find properties of the third edge
                        for PsL1=1:CurrPsNb-1
                            for PsL2=PsL1+1:CurrPsNb                                
                                if NbPaired{2}(PsL1,PsL2)==TestLimits(TestL)
                                    %search for an existing triangle
                                    if PsL2<CurrPsNb
                                        for PsL3=PsL2+1:CurrPsNb
                                            if NbPaired{2}(PsL1,PsL3)==TestLimits(TestL)
                                                TestStat{TestL}=[TestStat{TestL};[NbPaired{1}(PsL2,PsL3),NbPaired{2}(PsL2,PsL3)]];
                                            end
                                            if NbPaired{2}(PsL2,PsL3)==TestLimits(TestL)
                                                TestStat{TestL}=[TestStat{TestL};[NbPaired{1}(PsL1,PsL3),NbPaired{2}(PsL1,PsL3)]];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end                          
                if ~isempty(find(NbPaired{2}==NetNb))
                    %find properties of the third edge
                    for PsL1=1:CurrPsNb-1
                        for PsL2=PsL1+1:CurrPsNb
                            if NbPaired{2}(PsL1,PsL2)==NetNb
                                %search for an existing triangle
                                if PsL2<CurrPsNb
                                    for PsL3=PsL2+1:CurrPsNb
                                        if NbPaired{2}(PsL1,PsL3)==NetNb
                                            TestStat{TestNb}=[TestStat{TestNb};[NbPaired{1}(PsL2,PsL3),NbPaired{2}(PsL2,PsL3)]];                                        
                                        elseif NbPaired{2}(PsL2,PsL3)==NetNb
                                            TestStat{TestNb}=[TestStat{TestNb};[NbPaired{1}(PsL1,PsL3),NbPaired{2}(PsL1,PsL3)]];
                                        end                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


FigRank=42;
Legend={};
for TestL=1:TestNb
    if ~isempty(TestStat{TestL})
        Legend{1,end+1}=sprintf('delta %u (%u)',NetNb-TestLimits(TestL),length(TestStat{TestL}(:,2)));
    end
end
h=figure;
set(h,'name',sprintf('FIGU%u m%u: STAT ON THE THIRD EDGE IN G-G TRIANGLES (==)',FigRank,ChipRank))
set(gcf,'color',[1,1,1])
Colors=colors(colormap,TestNb);
for PlotL=1:2
    subplot(1,2,PlotL)
    hold on
    for TestL=1:TestNb
        if ~isempty(TestStat{TestL})
            CurrCorrNb=histc(TestStat{TestL}(:,PlotL),0:NetNb);
            plot(0:NetNb,CurrCorrNb/sum(CurrCorrNb),'color',Colors(TestL,:))
        end
    end
    set(gca,'box','on')
    xlabel('number of positive networks')
    ylabel('relative frequency')
    if PlotL==1
        legend(Legend)
        title('nb of network where third edges do not pass the test')
    else
        title('nb of network where third edges pass the test')
    end
end
Position=get(gcf,'position');
Position(3:4)=[1000,450];
set(gcf,'position',Position)
cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigRank),'png')
close(h)



%%  SUMMARIZE
function [NewPs,PsMatrix,PsType]=SUMMARIZE(ProbeNb,NewPs,Stat)
global K
PsNb=length(NewPs);
PsMade=zeros(length(NewPs),1);
PsType=zeros(5,6);

%% -- MERGING

%CORRECT SOURCE (GOP => 3)
for PsL=1:length(NewPs)
    for GeneL=1:length(NewPs{PsL}.geneNames)
        if ~isempty(findstr('GOP',NewPs{PsL}.geneNames{GeneL}))
            NewPs{PsL}.source(GeneL)=3;
            NewPs{PsL}.target(GeneL)=3;
        end
    end
end

ProbeNbs=[];
MaxProbeNb=0;
for PsL=1:length(NewPs)
    if ~isempty(NewPs{PsL}.geneNames)
        MaxProbeNb=max(MaxProbeNb,max(NewPs{PsL}.probeNb));
        ProbeNbs(end+1,1)=max(NewPs{PsL}.probeNb);
    end
end
PsMatrix=zeros(PsNb,17+ProbeNb);

Colors='bcrmgk';

%% ***** CLASS SS
CurrPsNb=Stat.catcount(1);
CurrPsRanks=Stat.index{1};


%FILL PsMatrix
CurrProbeNb=zeros(CurrPsNb,1);
for PsL=1:CurrPsNb    
    CurrPsRank=CurrPsRanks(PsL);
    PsMade(CurrPsRank)=1;
    CurrProbeNb(PsL)=NewPs{CurrPsRank}.probeNb;  
    %1st position: rank of assigned gene tot the current probe set
    PsMatrix(CurrPsRank,1)=NewPs{CurrPsRank}.geneRanks(1);
    %2nd position: pos of assigned gene in NewPs
    PsMatrix(CurrPsRank,2)=1;
    %3rd position: target type of assigned gene
    PsMatrix(CurrPsRank,3)=NewPs{CurrPsRank}.target;
    %4th position: source type of assigned gene
    PsMatrix(CurrPsRank,4)=NewPs{CurrPsRank}.source;    
    %5th position: nb of probe targetting the assigned gene
    PsMatrix(CurrPsRank,5)=NewPs{CurrPsRank}.probeNb;
    %8th position: nb of groups of transcripts corresponding to the assigned
    %gene
    PsMatrix(CurrPsRank,8)=1;
    %9th position: rank of the parent probe set
    PsMatrix(CurrPsRank,9)=CurrPsRank;
    %10th position: Rank of the group of transcripts targetted by the current
    %probe set in the assigned gene
    PsMatrix(CurrPsRank,10)=1;
    %17th position ClassRank
    PsMatrix(CurrPsRank,17)=1;
end

PsType=DIST_PROBENB(PsMatrix,CurrPsRanks,PsType,1);



%% ***** CLASS SM
CurrPsNb=Stat.catcount(2);
CurrPsRanks=Stat.index{2};


CurrProbeNb=zeros(CurrPsNb,1);
for PsL=1:CurrPsNb
    CurrProbeNb(PsL)=max(NewPs{CurrPsRanks(PsL)}.probeNb);
end

%STATISTICS ON PROBE NB AND GENE NB and START TO FILL PS

for PsL=1:length(CurrPsRanks)
    CurrPsRank=CurrPsRanks(PsL);    
    ProbeNbs=NewPs{CurrPsRank}.probeNb;
    FirstPNb(PsL)=max(ProbeNbs);
    IdemPos=find(ProbeNbs==FirstPNb(PsL));
    IdemNb=length(IdemPos);    
    %6th position : nb of not assigned genes targetted with the same nb (or greater nb) of
    %probes
    PsMatrix(CurrPsRank,6)=IdemNb-1;        
    %18th position and beyond: nb of  other genes targetted with all possible nb of
    %probes
    for NbL=1:length(ProbeNbs)
        PsMatrix(CurrPsRank,17+min(ProbeNbs(NbL),ProbeNb))=PsMatrix(CurrPsRank,17+min(ProbeNbs(NbL),ProbeNb))+1;
    end
    ProbeNbs(IdemPos)=[];    
    %7th position : nb of not assigned genes targetted with less nb of probes
    PsMatrix(CurrPsRank,7)=length(ProbeNbs);
end

%FILL PS
CurrProbeNb=zeros(CurrPsNb,1);
for PsL=1:length(CurrPsRanks)
    CurrPsRank=CurrPsRanks(PsL);
    PsMade(CurrPsRank)=1;
    %SELECT A GENE
    MaxProbeNb=max(NewPs{CurrPsRank}.probeNb);
    MaxPos=find(NewPs{CurrPsRank}.probeNb==MaxProbeNb);
    SelPos=1;
    if length(MaxPos)>1
        CurrTarget=NewPs{CurrPsRank}.target(MaxPos);
        CurrSource=NewPs{CurrPsRank}.source(MaxPos);
        if ~isempty(find(CurrTarget==1&CurrSource==1))
            SelPos=find(CurrTarget==1&CurrSource==1);
        elseif ~isempty(find(CurrTarget==1&CurrSource==2))
            SelPos=find(CurrTarget==1&CurrSource==2);
        elseif ~isempty(find(CurrTarget==2&CurrSource==1))
            SelPos=find(CurrTarget==2&CurrSource==1);
        elseif ~isempty(find(CurrTarget==2&CurrSource==2))
            SelPos=find(CurrTarget==2&CurrSource==2);
        end
        if length(SelPos)>1
            GeneNames=NewPs{CurrPsRank}.geneNames(MaxPos(SelPos));
            [GeneNames,SortIndex]=sort(GeneNames);
            SelPos=SelPos(SortIndex);
            SelPos=SelPos(1);
        end
    end

    SelPos=MaxPos(SelPos);
    CurrProbeNb(PsL)=NewPs{CurrPsRank}.probeNb(SelPos);
    %1st position: rank of assigned gene tot the current probe set
    PsMatrix(CurrPsRank,1)=NewPs{CurrPsRank}.geneRanks(SelPos);
    %2nd position: pos of assigned gene in NewPs
    PsMatrix(CurrPsRank,2)=SelPos;
    %3rd position: target type of assigned gene
    PsMatrix(CurrPsRank,3)=NewPs{CurrPsRank}.target(SelPos);
    %4th position: source type of assigned gene
    PsMatrix(CurrPsRank,4)=NewPs{CurrPsRank}.source(SelPos);
    %5th positon: nb of probe targetting the assigned gene
    PsMatrix(CurrPsRank,5)=NewPs{CurrPsRank}.probeNb(SelPos);
    %6th position : nb of not assigned genes targetted with the same nb (or greater nb) of probes
    ProbeNbs=NewPs{CurrPsRank}.probeNb;
    IdemPos=find(ProbeNbs>=PsMatrix(CurrPsRank,5));
    IdemNb=length(IdemPos);
    PsMatrix(CurrPsRank,6)=IdemNb-1;    
    if ~isempty(ProbeNbs)
         %18th position and beyond: nb of  other genes targetted with all possible nb of
        %probes
        for NbL=1:length(ProbeNbs)
            PsMatrix(CurrPsRank,17+min(ProbeNbs(NbL),ProbeNb))=PsMatrix(CurrPsRank,17+min(ProbeNbs(NbL),ProbeNb))+1;
        end
    end
    ProbeNbs(IdemPos)=[];
    %7th position : nb of not assigned genes targetted with less nb of probes
    PsMatrix(CurrPsRank,7)=length(ProbeNbs);
    %8th position: nb of groups of transcripts corresponding to the assigned gene
    PsMatrix(CurrPsRank,8)=1;
    %9th position: rank of the parent probe set
    PsMatrix(CurrPsRank,9)=CurrPsRank;
    %10th position: Rank of the group of transcripts targetted by the current
    %probe set in the assigned gene
    PsMatrix(CurrPsRank,10)=1;    
    %17th position ClassRank
    PsMatrix(CurrPsRank,17)=2;
end


% DISTRIBUTION PROBE NB, TARGET TYPE, GENE TYPE
PsType=DIST_PROBENB(PsMatrix,CurrPsRanks,PsType,2);


%% ***** CLASS TRUE MS
CurrPsNb=Stat.catcount(3);
CurrPsRanks=Stat.index{3};
GenePos=ones(CurrPsNb,1);
PsPos={};

CurrProbeNb=zeros(CurrPsNb,1);
for PsL=1:CurrPsNb
    CurrProbeNb(PsL)=max(NewPs{CurrPsRanks(PsL)}.probeNb);
end

%STATISTICS ON PROBE NB AND GENE NB and START TO FILL PS
[PsMatrix,PsMade]=fill_psmatrix(3,ProbeNb,NewPs,PsMatrix,PsMade,CurrPsRanks,Stat,GenePos,PsPos);
% DISTRIBUTION PROBE NB, TARGET TYPE, GENE TYPE
PsType=DIST_PROBENB(PsMatrix,CurrPsRanks,PsType,3);


%% ***** CLASSES TRUE MM, CX and HC are processed similarly
for ClassL=4:6

    CurrPsNb=Stat.catcount(ClassL);
    CurrPsRanks=Stat.index{ClassL};
    CurrPsMade=[];
    %Contruct GenePos
    GenePos=ones(CurrPsNb,1);
    PsPos=cell(CurrPsNb,1);
    for PsL1=1:CurrPsNb
        CurrPsRank1=CurrPsRanks(PsL1);
        if isempty(find(CurrPsMade==CurrPsRank1))
            %FIRST: find gene position
            if length(NewPs{CurrPsRank1}.geneNames)>1
                MaxProbeNb=max(NewPs{CurrPsRank1}.probeNb);
                MaxPos=find(NewPs{CurrPsRank1}.probeNb==MaxProbeNb);
                %privilegiate target type 1
                Targets=NewPs{CurrPsRank1}.target(MaxPos);
                Type1Pos=find(Targets==1);
                if ~isempty(Type1Pos)
                    MaxPos=MaxPos(Type1Pos);
                else
                    %privilegiate target type 2
                    Type2Pos=find(Targets==2);
                    if ~isempty(Type2Pos)
                        MaxPos=MaxPos(Type2Pos);
                    end
                end
                SelPos=1;
                if length(MaxPos)>1
                    %calculate ratio nb of ps/nb of grps
                    PsRanks=NewPs{CurrPsRank1}.psRanks(MaxPos);
                    PsNbs=zeros(length(MaxPos),1);
                    for PsL=1:length(MaxPos)
                        PsNbs(PsL)=length(PsRanks{PsL});
                    end
                    Grps=NewPs{CurrPsRank1}.grp(MaxPos);
                    GrpNbs=zeros(length(MaxPos),1);
                    for PsL=1:length(MaxPos)
                        GrpNbs(PsL)=length(Grps{PsL});
                    end
                    %select max PsNbs/GrpNbs
                    MaxRatio=max(PsNbs./GrpNbs);
                    MaxRatioPos=find(PsNbs./GrpNbs==MaxRatio);
                    if length(MaxRatioPos)>1
                        %select min PsNbs => increases the number of targeted genes
                        PsNbs=PsNbs(MaxRatioPos);
                        MinPsNb=min(PsNbs);
                        MinPsPos=find(PsNbs==MinPsNb);
                        MaxPos=MaxPos(MaxRatioPos(MinPsPos));
                    else
                        MaxPos=MaxPos(MaxRatioPos);
                    end
                    if length(MaxPos)>1
                        CurrTarget=NewPs{CurrPsRank1}.target(MaxPos);
                        CurrSource=NewPs{CurrPsRank1}.source(MaxPos);
                        if ~isempty(find(CurrTarget==1&CurrSource==1))
                            SelPos=find(CurrTarget==1&CurrSource==1);
                        elseif ~isempty(find(CurrTarget==1&CurrSource==2))
                            SelPos=find(CurrTarget==1&CurrSource==2);
                        elseif ~isempty(find(CurrTarget==2&CurrSource==1))
                            SelPos=find(CurrTarget==2&CurrSource==1);
                        elseif ~isempty(find(CurrTarget==2&CurrSource==2))
                            SelPos=find(CurrTarget==2&CurrSource==2);
                        end
                        if length(SelPos)>1
                            %take first in alphabetical order
                            GeneNames=NewPs{CurrPsRank1}.geneNames(MaxPos(SelPos));
                            [GeneNames,SortIndex]=sort(GeneNames);
                            SelPos=SelPos(SortIndex);
                            SelPos=SelPos(1);
                        end
                    end
                end
                GenePos(PsL1)=MaxPos(SelPos);
            else
                GenePos(PsL1)=1;
            end
            CurrGeneName1=NewPs{CurrPsRank1}.geneNames{GenePos(PsL1)};

            %SND: find used probe set position in selected gene
            if length(NewPs{CurrPsRank1}.psRanks{GenePos(PsL1)})>1
                PsPos{PsL1}=[];
                %scan all other probe set targeting the current gene to see if the current gene
                %is the one selected with the same criterium
                for PsL2=1:length(NewPs{CurrPsRank1}.psRanks{GenePos(PsL1)})
                    CurrPsRank2=NewPs{CurrPsRank1}.psRanks{GenePos(PsL1)}(PsL2);
                    if isempty(find(CurrPsMade==CurrPsRank2))
                        if CurrPsRank2==CurrPsRank1
                            PsPos{PsL1}=[PsPos{PsL1},PsL2];
                        else
                            if length(NewPs{CurrPsRank2}.geneNames)>1
                                CurrGenePos=find(NewPs{CurrPsRank2}.geneRanks==NewPs{CurrPsRank1}.geneRanks(GenePos(PsL1)));
                                MaxProbeNb=max(NewPs{CurrPsRank2}.probeNb);
                                MaxPos=find(NewPs{CurrPsRank2}.probeNb==MaxProbeNb);
                                %privilegiate target type 1
                                Targets=NewPs{CurrPsRank2}.target(MaxPos);
                                Type1Pos=find(Targets==1);
                                if ~isempty(Type1Pos)
                                    MaxPos=MaxPos(Type1Pos);
                                else
                                    %privilegiate target type 2
                                    Type2Pos=find(Targets==2);
                                    if ~isempty(Type2Pos)
                                        MaxPos=MaxPos(Type2Pos);
                                    end
                                end
                                if ~isempty(find(MaxPos==CurrGenePos))
                                    if length(MaxPos)>1
                                        %calculate ratio nb of ps/nb of grps
                                        PsRanks=NewPs{CurrPsRank2}.psRanks(MaxPos);
                                        PsNbs=zeros(length(MaxPos),1);
                                        for PsL=1:length(MaxPos)
                                            PsNbs(PsL)=length(PsRanks{PsL});
                                        end
                                        Grps=NewPs{CurrPsRank2}.grp(MaxPos);
                                        GrpNbs=zeros(length(MaxPos),1);
                                        for PsL=1:length(MaxPos)
                                            GrpNbs(PsL)=length(Grps{PsL});
                                        end
                                        %select max PsNbs/GrpNbs
                                        MaxRatio=max(PsNbs./GrpNbs);
                                        MaxRatioPos=find(PsNbs./GrpNbs==MaxRatio);
                                        PsRanks=[];
                                        for i=1:length(MaxRatioPos)
                                            PsRanks=[PsRanks,NewPs{CurrPsRank2}.psRanks{MaxRatioPos(i)}];
                                        end
                                        if ~isempty(find(PsRanks==CurrPsRank2))
                                            if length(MaxRatioPos)>1
                                                %select min PsNbs
                                                PsNbs=PsNbs(MaxRatioPos);
                                                MinPsNb=min(PsNbs);
                                                MinPsPos=find(PsNbs==MinPsNb);
                                                MaxPos=MaxPos(MaxRatioPos(MinPsPos));
                                            else
                                                MaxPos=MaxPos(MaxRatioPos);
                                            end
                                            SelPos=1;
                                            if length(MaxPos)>1
                                                CurrTarget=NewPs{CurrPsRank2}.target(MaxPos);
                                                CurrSource=NewPs{CurrPsRank2}.source(MaxPos);
                                                if ~isempty(find(CurrTarget==1&CurrSource==1))
                                                    SelPos=find(CurrTarget==1&CurrSource==1);
                                                elseif ~isempty(find(CurrTarget==1&CurrSource==2))
                                                    SelPos=find(CurrTarget==1&CurrSource==2);
                                                elseif ~isempty(find(CurrTarget==2&CurrSource==1))
                                                    SelPos=find(CurrTarget==2&CurrSource==1);
                                                elseif ~isempty(find(CurrTarget==2&CurrSource==2))
                                                    SelPos=find(CurrTarget==2&CurrSource==2);
                                                end
                                                if length(SelPos)>1
                                                    GeneNames=NewPs{CurrPsRank2}.geneNames(MaxPos(SelPos));
                                                    [GeneNames,SortIndex]=sort(GeneNames);
                                                    SelPos=SelPos(SortIndex);
                                                    SelPos=SelPos(1);
                                                end
                                            end
                                            CurrGeneName2=NewPs{CurrPsRank2}.geneNames{MaxPos(SelPos)};
                                            %verify if the same gene is assigned
                                            if isequal(CurrGeneName1,CurrGeneName2)
                                                PsPos{PsL1}=[PsPos{PsL1},PsL2];
                                            end
                                        end
                                    else
                                        CurrGeneName2=NewPs{CurrPsRank2}.geneNames{MaxPos};
                                        %verify if the same gene is assigned
                                        if isequal(CurrGeneName1,CurrGeneName2)
                                            PsPos{PsL1}=[PsPos{PsL1},PsL2];
                                        end
                                    end
                                end
                            else
                                %only one targeted gene by PsL2
                                PsPos{PsL1}=[PsPos{PsL1},PsL2];
                            end
                        end
                    end
                end
            else
                PsPos{PsL1}=1;                
            end
            CurrPsMade=union(CurrPsMade,NewPs{CurrPsRank1}.psRanks{GenePos(PsL1)}(PsPos{PsL1}));
        end
    end
    %STATISTICS ON PROBE NB AND GENE NB and START TO FILL PS
    [PsMatrix,PsMade]=fill_psmatrix(ClassL,ProbeNb,NewPs,PsMatrix,PsMade,CurrPsRanks,Stat,GenePos,PsPos);
    % DISTRIBUTION PROBE NB, TARGET TYPE, GENE TYPE
    PsType=DIST_PROBENB(PsMatrix,CurrPsRanks,PsType,ClassL);
end
for ClassL=1:7
    PsMatrix(Stat.index{ClassL},17)=ClassL;
end
%fill parent probe set for class 7
PsMatrix(find(PsMatrix(:,17)==7),9)=find(PsMatrix(:,17)==7);


%% DISPLAY_PS (FIG55a,FIG55b)
function DISPLAY_PS(Species,ChipRank,EnsExonGeneNbs,AceExonGeneNbs,ProbeNbLimit)
global K

EnsGeneNbsSup=sum(EnsExonGeneNbs(:,ProbeNbLimit:end),2);
EnsGeneNbsInf=sum(EnsExonGeneNbs(:,1:ProbeNbLimit-1),2);

if ~isempty(AceExonGeneNbs)
    AceGeneNbsSup=sum(AceExonGeneNbs(:,ProbeNbLimit:end),2);
    AceGeneNbsInf=sum(AceExonGeneNbs(:,1:ProbeNbLimit-1),2);
    AceFlag=1;
    XPlot=2;
else
    AceFlag=0;
    XPlot=1;
end

%histogram of nb of genes with more (or less) ProbeNbLimit probes
if ProbeNbLimit==1
    YPlot=1;
else
    YPlot=2;
end
h(1)=figure;
FigRank=55;
set(h(1),'name',sprintf('FIG%ua - m%u: Nb of targeted genes',FigRank,ChipRank))
set(gcf,'color',[1,1,1])
subplot(YPlot,XPlot,1)
N=histc(EnsGeneNbsSup,0:20);
plot(2:20,N(3:end),'r')
hold on
plot(2:20,N(3:end),'ro')
xlabel({'nb of targeted genes';sprintf('%u ps target no gene; %u ps target one gene',N(1),N(2))})
ylabel('probe set frequency')
title(sprintf('Nb of Ensembl genes \ntargeted by at least %u probes',ProbeNbLimit))

if AceFlag
    subplot(YPlot,XPlot,2)
    N=histc(AceGeneNbsSup,0:20);
    plot(2:20,N(3:end),'m')
    hold on
    plot(2:20,N(3:end),'mo')
    xlabel({'nb of targeted genes';sprintf('%u ps target no gene; %u ps target one gene',N(1),N(2))})
    ylabel('probe set frequency')
    title(sprintf('Nb of AceView genes \ntargeted by at least %u probes',ProbeNbLimit))
end

if YPlot==2
   if AceFlag 
       subplot(YPlot,XPlot,3)
   else
       subplot(YPlot,XPlot,2)
   end
    N=histc(EnsGeneNbsInf,0:40);
    plot(2:40,N(3:end),'r')
    hold on
    plot(2:40,N(3:end),'ro')
    xlabel({'nb of targeted genes';sprintf('%u ps target no gene; %u ps target one gene',N(1),N(2))})
    ylabel('probe set frequency')
    title(sprintf('Nb of Ensembl genes \ntargeted by less than %u  probes',ProbeNbLimit))
    if AceFlag
        subplot(YPlot,XPlot,4)
        N=histc(AceGeneNbsInf,0:40);
        plot(2:40,N(3:end),'m')
        hold on
        plot(2:40,N(3:end),'mo')
        xlabel({'nb of targeted genes';sprintf('%u ps target no gene; %u ps target one gene',N(1),N(2))})
        ylabel('probe set frequency')
        title(sprintf('Nb of AceView genes \ntargeted by less than %u probes',ProbeNbLimit))
    end
end

h(2)=figure;
set(h(2),'name',sprintf('FIG%ub - m%u:  Nb of targeting probes',FigRank,ChipRank))
set(gcf,'color',[1,1,1])
plot(sum(EnsExonGeneNbs),'r');
hold on
plot(sum(EnsExonGeneNbs),'ro');
if AceFlag
    plot(sum(AceExonGeneNbs),'m');
    hold on
    plot(sum(AceExonGeneNbs),'mo');
    ylabel('gene frequency (red: Ensembl, magenta: AceView)')
else
    ylabel('gene frequency (red: Ensembl)')
end
xlabel('nb of targeting probes')


Letters='ab';
for FigL=1:2
    figure(h(FigL))
    Position=get(gcf,'position');
    Position(3:4)=[600,450];
    set(gcf,'position',Position)
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    saveas(h(FigL),sprintf('m%u_fig%u%c.png',ChipRank,FigRank,Letters(FigL)),'png')
    close(h(FigL))
end


%% DIST_PROBENB
function PsType=DIST_PROBENB(PsMatrix,CurrPsRanks,PsType,ClassRank)

% PERCENTAGE OF DIFFERENT TARGET TYPES AND GENE TYPES
for i=1:3
    Pos=find(PsMatrix(CurrPsRanks,3)==i);
    PsType(ClassRank,i)=round(length(Pos)*100/length(CurrPsRanks));
end
PsType(ClassRank,4)=round(length(find(PsMatrix(CurrPsRanks,4)==1))*100/length(CurrPsRanks));
PsType(ClassRank,5)=round(length(find(PsMatrix(CurrPsRanks,4)==2))*100/length(CurrPsRanks));














%==========================%
% FUNCTION IMPORT_TARGETNB %
%==========================%
%
%IMPORT_TARGETNB read a text file (containing either Ensembl or AceView informations)
%  and creates a matrix indicating for each probeset the number of genes that have x
%  probes targeting their exons (with x>=1 and x<=n(max(probe nb)))
%
%INPUT PARAMETERS
%  ChipRank: rank of chip model
%
%EXTERNAL FILES
%  Read file 'm%ChipRank_probesets_ensembl.txt'
%  and eventually 'm%ChipRank_probesets_aceview.txt'
%
%  File format: [Ps Rank, number of targeted genes with n probes, n-1 probes, ...,1 probe]
%               ex : [54, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 2, 0, 0, 2, 0]
%
%  Write 'm%ChipRank_probesets_ensembl.mat' containing variable EnsExonGeneNbs
%  and eventually 'm%ChipRank_probesets_aceview.mat'  containing variable AceExonGeneNbs
%  In output variable, probe nb are in the direct order (1,2,...,n)
%
%FIGURES 16 to 18

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



function import_targetnb(ChipRank,DisplayFlag)
global K

ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
if isempty(ChipPos)
    h=errordlg(sprintf('chip m%u does not exist',ChipRank));
    waitfor(h)
    error('process canceled')
end
Species=K.chip.species{ChipPos};

if ~exist(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)),'dir')
    mkdir((fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank))))
end

cd(fullfile(K.dir.pydata,Species,'txt'))
if exist(fullfile(K.dir.pydata,Species,'txt',sprintf('m%u_probesets_ensembl.txt',ChipRank)),'file')
    Lines=textread(sprintf('m%u_probesets_ensembl.txt',ChipRank),'%s','delimiter','\t');
    FirstLine=eval(Lines{1});    
    GeneNbs=zeros(length(Lines),length(FirstLine));
    for LineL =1:length(Lines)
        GeneNbs(LineL,:)=eval(Lines{LineL});
    end
    %sort according to probe set ranks
    % Probe set ranks start at zero (python output)
    PsRank=GeneNbs(:,1)+1;
    EnsExonGeneNbs=GeneNbs(:,2:end);
    [PsRank SortIndex]=sort(PsRank);
    EnsExonGeneNbs=EnsExonGeneNbs(SortIndex,:);
    %flip to have probe nb in order
    EnsExonGeneNbs=fliplr(EnsExonGeneNbs);
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)))       
    eval(sprintf('save m%u_probesets_ensembl EnsExonGeneNbs',ChipRank))
    if K.chip.probesetNb(ChipPos)~=size(EnsExonGeneNbs,1);
        h=warndlg('%u ps in K.chip and %u ps in PSAWNpy');
        waitfor(h)
    end
    K.chip.maxProbeNb(ChipPos,1)=size(EnsExonGeneNbs,2);
    cd(K.dir.mlprog)
    save global_ps K
    if DisplayFlag
        %FIG 16
        DISPLAY_TARGNB(EnsExonGeneNbs,Species,ChipRank,sprintf('FIG1 - ENSEMBL TARGET NB - CHIP m%u',ChipRank),16)
    end
else    
    h=errordlg(sprintf('m%u_probesets_ensembl.txt does ot exist',ChipRank));
    waitfor(h)
    error('process canceled')
end


cd(fullfile(K.dir.pydata,Species,'txt'))
if exist(fullfile(K.dir.pydata,Species,'txt',sprintf('m%u_probesets_aceview.txt',ChipRank)),'file')
    Lines=textread(sprintf('m%u_probesets_aceview.txt',ChipRank),'%s','delimiter','\t');
    FirstLine=eval(Lines{1});
    GeneNbs=zeros(length(Lines),length(FirstLine));
    for LineL =1:length(Lines)
        GeneNbs(LineL,:)=eval(Lines{LineL});
    end
    %sort according to probeset ranks
    % Probe set ranks start at zero (python output)
    PsRank=GeneNbs(:,1)+1;
    AceExonGeneNbs=GeneNbs(:,2:end);
    [PsRank SortIndex]=sort(PsRank);
    AceExonGeneNbs=AceExonGeneNbs(SortIndex,:);
    %flip to have probe nb in order
    AceExonGeneNbs=fliplr(AceExonGeneNbs);
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank)))
    eval(sprintf('save m%u_probesets_aceview AceExonGeneNbs',ChipRank))
    if DisplayFlag
        %FIG17
        DISPLAY_TARGNB(AceExonGeneNbs,Species,ChipRank,sprintf('FIG2 - ACEVIEW TARGET NB - CHIP m%u',ChipRank),17)
        %FIG18
        DISPLAY_ACE(Species,ChipRank,18)
    end
end

%% function DISPLAY_TARGNB
function DISPLAY_TARGNB(Data,Species,ChipRank,Title,FigNb)
global K

h=figure;
MaxProbeNb=size(Data,2);
set(h,'name',Title)
set(gcf,'color',[1,1,1])
TargNb=sum(Data);
subplot(1,2,1)
plot(1:MaxProbeNb,TargNb,'b-')
hold on
plot(1:MaxProbeNb,TargNb,'bo')
xlabel('exact nb of probes targeting the gene')
ylabel({'total nb of targeted genes (blue)';'nb of genes that are the single target of a probe set (cyan)'})
GeneNb=sum(Data,2);
PsPos=find(GeneNb==1);
SingleData=Data(PsPos,:);
SingleTargNb=sum(SingleData);
plot(1:MaxProbeNb,SingleTargNb,'c-')
plot(1:MaxProbeNb,SingleTargNb,'co')
subplot(1,2,2)
hold on
%use most current probenb for that chip model
[temp,ProbeNb]=max(TargNb);
Colors=colors(colormap,ProbeNb+1);
SingleGeneNb=zeros(1,ProbeNb+1);
for ProbeL=1:ProbeNb
    GeneNb=sum(Data(:,ProbeL:MaxProbeNb),2);
    GeneNbVal=unique(GeneNb);
    ProbeSetNb=histc(GeneNb,GeneNbVal);
    SinglePos=find(GeneNbVal==1);
    if ~isempty(SinglePos)
        SingleGeneNb(ProbeL)=ProbeSetNb(SinglePos);
    end
    plot(GeneNbVal,ProbeSetNb,'color',Colors(ProbeL,:))
end
%group other genes targeted by more than ProbeNb
SingleGeneNb(ProbeNb+1)=0;
for ProbeL=ProbeNb+1:MaxProbeNb;
    GeneNb=sum(Data(:,ProbeL:MaxProbeNb),2);
    GeneNbVal=unique(GeneNb);
    ProbeSetNb=histc(GeneNb,GeneNbVal);
    SinglePos=find(GeneNbVal==1);
    if ~isempty(SinglePos)
        SingleGeneNb(ProbeNb+1)=SingleGeneNb(ProbeNb+1)+ProbeSetNb(SinglePos);
    end
end
plot(GeneNbVal,ProbeSetNb,'color',Colors(ProbeNb+1,:))



set(gca,'box','on')
xlabel('nb of targeted genes')
ylabel('nb of probe sets')
Legend=cell(ProbeNb,1);
for ProbeL=1:ProbeNb+1
    Legend{ProbeL}=sprintf('>=%u probes (%u)',ProbeL,SingleGeneNb(ProbeL));
end
legend(Legend)
set(gca,'xlim',[0,4])
Position=get(gcf,'position');
Position(3:4)=[800,500];
set(gcf,'position',Position)
try
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
catch
    mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))    
end
saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigNb),'png')
plot2svg(sprintf('m%u_fig%u.svg',ChipRank,FigNb))
close(h)
%set(gca,'yscale','log')
    
%% function DISPLAY_ACE
function DISPLAY_ACE(Species,ChipRank,FigNb)

global K
cd(fullfile(K.dir.pydata,Species,'txt'))
if exist(fullfile(K.dir.pydata,Species,'txt',sprintf('%s_ens_by_ace_gene.txt',Species)),'file')
    [Ace,Ens]=textread(sprintf('%s_ens_by_ace_gene.txt',Species),'%s%s','delimiter','\t');
end
%list of AceView genes
AceGenes=unique(Ace);
%remove Ace ='0' corresponding to Ensembl genes without
%match in AceView
AceGenes(strmatch('0',AceGenes))=[];
%nb of Ensembl genes corresponding to a particular AceView gene
EnsGeneNb=zeros(length(AceGenes),1);
%list of Ensembl genes
EnsGenes={};
%nb of AceView genes corresponding to a particular Ensembl gene
AceGeneNb=[];
for AceL=1:length(Ace)
    if isequal(Ace{AceL},'0')
        %process Ensembl gene without match in AceView
        EnsGene=Ens{AceL};
        EnsPos=strmatch(EnsGene,EnsGenes,'exact');
        if ~isempty(EnsPos)     
            if AceGeneNb(EnsPos)>0
                h=errordlg(sprintf('non reciprocity for Ensembl gene %s and AceView gene %s',EnsGene,Ace{EnsPos}));
                waitfor(h)
                error('process canceled')
            end
        else
            EnsGenes{end+1,1}=EnsGene;
            AceGeneNb(end+1,1)=0;
        end
    else    
        Genes=eval(Ens{AceL});
        if ~isempty(Genes)
            %process Ensembl gene with match in AceView
            EnsGeneNb(AceL)=length(Genes);
            for GeneL=1:length(Genes)
                EnsPos=strmatch(Genes{GeneL},EnsGenes,'exact');
                if ~isempty(EnsPos)                
                    AceGeneNb(EnsPos)=AceGeneNb(EnsPos)+1;
                else
                    EnsGenes{end+1,1}=Genes{GeneL};
                    AceGeneNb(end+1,1)=1;
                end
            end
        end
    end
end

AceGeneNbs=unique(AceGeneNb);
EnsGeneNbs=unique(EnsGeneNb);

AceVal=histc(AceGeneNb,AceGeneNbs);
EnsVal=histc(EnsGeneNb,EnsGeneNbs);

h=figure;
set(h,'name','FIG3 - PROBE SETS TARGETING X GENES')
set(gcf,'color',[1,1,1])
plot(AceGeneNbs,AceVal,'bo')
hold on
h1=plot(AceGeneNbs,AceVal,'b-');
plot(EnsGeneNbs,EnsVal,'co')
h2=plot(EnsGeneNbs,EnsVal,'c-');
set(gca,'xlim',[0,10])
xlabel('nb of matching genes')
ylabel('nb of genes with x matching genes')
Legend{1}=sprintf('%u Ensembl genes',length(AceGeneNb));
Legend{2}=sprintf('%u AceView genes',length(EnsGeneNb));
legend([h1,h2],Legend)
Position=get(gcf,'position');
Position(3:4)=[500,400];
set(gcf,'position',Position)
try
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
catch
    mkdir(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))
    cd(fullfile(K.dir.mldata,Species,sprintf('m%u',ChipRank),'png'))    
end
saveas(h,sprintf('m%u_fig%u.png',ChipRank,FigNb),'png')
plot2svg(sprintf('m%u_fig%u.svg',ChipRank,FigNb))
close(h)

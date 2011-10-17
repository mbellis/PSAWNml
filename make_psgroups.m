%========================%
% FUNCTION MAKE_PSGROUPS %
%========================%
%
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

% MAKE_PSGROUPS partitionate a series of probe sets in groups that can be considered
%   as targeting the same transcripts, based on their properties in
%   a particular network or in a set of network

% INPUT
% 1        PsRanks: ranks of the group of probe sets
% 2      PairedMat: set of vectors indicating the type of relation between each
%                   possible couple of probe set in each networks (0: no correlation,
%                   1: not significant correlation, 2: significant correlation)
% 3          NetNb: nb of networks used
% 4      TestLimit: the minimal number of network where a correlation values must exist
%                   (either significative and coded by 2, or not significative and coded by 1
%                   in PairedMat)
% 5 GrpSizeNb:
%OUTPUT
%   PsGrp: groups of PsRanks

%VERSION
%   V03 - 23 07 2010 - Novel approach with count of networks with
%   significative correlation for a given pair of probe set and using
%   agregation of triangle
%   V01 - 2009: only one network
%   V02 - 2010-3-5 = several networks
function [PsGrp,varargout]=make_psgroups(PsRanks,PairedMat,NetNb,TestLimit,GrpSizeNb,LinkType)

if nargout==5
    GrpSizes=zeros(1,GrpSizeNb);
    LinkTypes=[];
    BadLink=[];
    Hubs=[];
end


%construct a matrix of paired probe sets by using information on the number of network
% with existing correlation between the paired probe sets
PairNb=size(PairedMat,1);
%Transform cel into matrix independant of link type (=> significatvie links)
NetMat=zeros(PairNb,NetNb);
for PairL=1:PairNb
    NetMat(PairL,:)=PairedMat{PairL};
end
NetMat(find(NetMat<LinkType))=0;
NetMat(find(NetMat))=1;

%find one of the largest set of pairs that have >= TestLimit significative 
%links in the same networks
MemPairs=[];
for PairL1=1:PairNb
    CurrPair=NetMat(PairL1,:);
    Pairs=[];
    for PairL2=1:PairNb
        if length(find(CurrPair&NetMat(PairL2,:)))>=TestLimit
            Pairs=[Pairs;PairL2];
        end
    end
    if length(Pairs)>length(MemPairs)
        MemPairs=Pairs;
    end
end
Paired=zeros(1,PairNb);
Paired(MemPairs)=1;

% for PsL=1:PairNb
%     if length(find(PairedMat{PsL}>=LinkType))>=TestLimit
%         Paired(PsL)=1;
%     end
% end
%transform this list into a symetric matrix (pairs are arranged in PairedMat
%in a way taht allows this transformation
Paired=squareform(Paired);
PsNb=length(Paired);

if nargout==5
    %recover the type of the pairs
    %either 2, if correlation values are all significative in all the networks
    %or 1 if the number of networks with a positive correlation is greater or equal
    %to TestLimit
    PairedType=zeros(1,PairNb);
    for PairL=1:PairNb
        if sum(PairedMat{PairL})==NetNb*2
            PairedType(PairL)=2;
        %elseif length(find(PairedMat{PairL}>=LinkType))>=TestLimit
        elseif ~isempty(find(MemPairs==PairL))
            PairedType(PairL)=1;
        end
    end
    PairedType=squareform(PairedType);
end

%ALGORITHM FOR SEARCHING TRIANGLES
%construct a dictionnary of neighbhoors
Nodes=cell(PsNb,1);
for PsL=1:PsNb
    Nodes{PsL}=find(Paired(PsL,:));
end

%find existing triangle
Triangle=[];
%nodes must have been already visited before being
%eventually searched for the presence of a triangle
%=>prevents from writting doublons, and reduces search time
Visited=[];
for FirstNode=1:PsNb
    if ~isempty(Nodes{FirstNode})
        for PsL2=1:length(Nodes{FirstNode})
            SndNode=Nodes{FirstNode}(PsL2);
            if ~isempty(find(Visited==SndNode))
                for PsL3=1:length(Nodes{SndNode})
                    ThirdNode=Nodes{SndNode}(PsL3);
                    if ~isempty(find(Visited==ThirdNode))
                        if ~isempty(find(Nodes{ThirdNode}==FirstNode))
                            Triangle=[Triangle;sort([FirstNode,SndNode,ThirdNode])];
                        end
                    end
                end
            end
        end
        Visited=[Visited,FirstNode];
    end
end
Triangle=unique(Triangle,'rows');


PsGrp={};
Control=[];
PsGrpNb1=0;


TriangleNb=size(Triangle,1);
if TriangleNb>0
    if TriangleNb==1
        PsGrp{1}=Triangle(1,:);
    else
        Continue=1;

        Processed=zeros(TriangleNb,1);
        TriL1=1;
        TriL2=2;

        while Continue
            if Processed(TriL1)==0|Processed(TriL2)==0
                %process current pair of triangles
                %                 TriL1
                %                 TriL2
                if length(intersect(Triangle(TriL1,:),Triangle(TriL2,:)))==2
                    %the two triangles have an edge in common => merge them
                    %point specific to triangle 1
                    %!!! Point1=setdiff(Triangle(TriL1,:),Triangle(TriL2,:));
                    %point specific to triangle 2
                    %!!! Point2=setdiff(Triangle(TriL2,:),Triangle(TriL1,:));
                    %common points
                    %!!! Point3=intersect(Triangle(TriL1,:),Triangle(TriL2,:));

                    %!!!NewTriangle(1,:)=sort([Point1,Point2,Point3(1)]);
                    %!!!NewTriangle(2,:)=sort([Point1,Point2,Point3(2)]);

                    %!!!
                    %PsGrpPos=length(PsGrp)+1;
                    %if ~isempty(Processed(TriL1))!!!
                    %    PsGrp{PsGrpPos,1}=unique([Triangle(TriL1,:),Triangle(TriL2,:)]);
                    %else
                    %    PsGrp{PsGrpPos,1}=unique([PsGrp{PsGrpPos},Triangle(TriL2,:)]);
                    %end
                    %!!!

                    %PsGrpPos indicates the ps group which contains Triangle TriL1
                    if Processed(TriL1)==0
                        PsGrpPos=length(PsGrp)+1;
                        PsGrp{PsGrpPos,1}=unique([Triangle(TriL1,:),Triangle(TriL2,:)]);
                    else
                        PsGrp{PsGrpPos,1}=unique([PsGrp{PsGrpPos},Triangle(TriL2,:)]);
                    end

                    Processed(TriL1)=1;
                    Processed(TriL2)=1;

                    %add eventually Triangle to current ps group if they have an edge in common
                    AddFlag=1;
                    TriIndex=1:TriangleNb;
                    %!!!TriIndex1=TriIndex;
                    IsAdded=0;
                    while AddFlag
                        AddFlag=0;
                        %!!!for TriL3=length(TriIndex1):-1:1
                        for TriL3=length(TriIndex):-1:1
                            %!!!if Processed(TriIndex(TriIndex1(TriL3)))==0
                            if Processed(TriIndex(TriL3))==0
                                %!!!if length(intersect(Triangle(TriIndex(TriIndex1(TriL3)),:),PsGrp{PsGrpPos}))>=2
                                if length(intersect(Triangle(TriIndex(TriL3),:),PsGrp{PsGrpPos}))>=2
                                    %!!!!!!PsGrp{PsGrpPos,1}=unique([PsGrp{PsGrpPos},Triangle(TriL3,:)]);
                                    PsGrp{PsGrpPos,1}=unique([PsGrp{PsGrpPos},Triangle(TriIndex(TriL3),:)]);
                                    IsAdded=1;
                                    %!!!!!!Processed(TriL3)=1;
                                    Processed(TriIndex(TriL3))=1;
                                    %suppress the processed triangle from TriIndex
                                    TriIndex(TriL3)=[];

                                    %process existing ps groups
                                    if ~isempty(PsGrp)
                                        [PsGrp,PsGrpPos]=MERGE_PSGROUPS(PsGrp,PsGrpPos);
                                    end

                                    if ~isempty(TriIndex)
                                        AddFlag=1;
                                    end
                                    %each time a triangle is added, start again from the
                                    %begining (merging of triangles or of ps groups may create
                                    %new possibilities of merging previously tested triangles
                                    break
                                end
                            end
                        end
                    end

                    if IsAdded==0
                        %process ps groups
                        if ~isempty(PsGrp)
                            [PsGrp,PsGrpPos]=MERGE_PSGROUPS(PsGrp,PsGrpPos);
                        end
                    end
                end %if found a common edge between two current triangles
            end %if current pair of triangles is to be processed

            %Test if subsequent TriL2 triangles have been processed
            if TriL2<TriangleNb
                if sum(Processed(TriL2+1:end))==TriangleNb-TriL2
                    TriL2=TriangleNb;
                end
            end
            if length(find(Processed))>=TriangleNb|(TriL1==TriangleNb-1&TriL2==TriangleNb)
                Continue=0;
            else
                %skip Triangle that have been already processed
                %increment TriL1
                if TriL2==TriangleNb
                    TriL1=TriL1+1;
                    %search next 0 position for TriL1
                    if Processed(TriL1)&TriL1<TriangleNb-1
                        for TriPos=TriL1+1:TriangleNb-1
                            if Processed(TriPos)==0
                                TriL1=TriPos;
                                break
                            end
                        end
                        if TriPos==TriangleNb-1&Processed(TriPos)==1
                            TriL1=TriangleNb-1;
                        end
                    elseif Processed(TriL1)&TriL1==TriangleNb-1
                        TriL1=TriangleNb-1;
                    end
                    %restart TriL2
                    %!!!if Continue
                    %!!!if Processed(1)
                    % exist at leat one zero position for TriL2
                    for TriPos=1:TriangleNb
                        if Processed(TriPos)==0
                            TriL2=TriPos;
                            break
                        end
                    end
                    %!!!else
                    %!!!    TriL2=1;
                    %!!!end
                    %end
                else
                    %%%TriL2=TriL2+1;
                    %!!!if Processed(TriL2)
                    %!!!for TriPos=TriL2+1:TriangleNb
                    for TriPos=TriL2+1:TriangleNb
                        % exist at leat one zero position for TriL2
                        % otherwise TriL2 should be set to TriangleNb
                        if Processed(TriPos)==0
                            TriL2=TriPos;
                            break
                        end
                    end
                    %end;
                end
            end
        end %while Continue

        %process not used triangles
        NotUsed=find(Processed==0);
        if ~isempty(NotUsed)
            for i=1:length(NotUsed)
                PsGrp{end+1,1}=Triangle(NotUsed(i),:);
            end
            %process ps groups
            if ~isempty(PsGrp)
                %PsGrpPos not used (might be not defined at this step)
                [PsGrp,Temp]=MERGE_PSGROUPS(PsGrp,0);
            end
        end
    end %construction ps groups from triangle

    PsGrpNb1=length(PsGrp);
    if PsGrpNb1>1
        %split if there is only one point in common
        %between two ps groups (two points in common
        %is not possible, because in this case, ps groups have been merged)
        PsGrpNb1=length(PsGrp);
        MemPsGrp1=PsGrp;
        SplittedFlag=1;
        while SplittedFlag
            SplittedFlag=0;
            for i=1:PsGrpNb1-1               
                for j=i+1:PsGrpNb1
                    Common=intersect(PsGrp{i},PsGrp{j});
                    if length(Common)==1
                        SplittedFlag=1;
                        %eliminate the common ps from the two intersecting groups
                        MakeNew=1;
                        if length(PsGrp{i})>1
                            PsGrp{i}=setdiff(PsGrp{i},Common);
                        else 
                            PsGrpNb1=i;
                            PsGrpNb2=j;
                            MakeNew=0;
                        end
                        if length(PsGrp{j})>1
                            PsGrp{j}=setdiff(PsGrp{j},Common);
                        else
                            PsGrpNb1=j;
                            PsGrpNb2=i;
                            MakeNew=0;
                        end
                        %create a new single ps grp => hub ps cannot be in the two first grps
                        if MakeNew                            
                            PsGrp{end+1,1}=Common;                            
                            PsGrpNb1=PsGrpNb1+1;
                            Hubs=[Hubs;[PsGrpNb1,i,j,0,0]];
                            EndPos=size(Hubs,1);
                            %add new relation if a group have been splitted
                            if EndPos>1
                                for k=[i,j]
                                    for m=[2,3]
                                        Pos=find(Hubs(1:EndPos-1,m)==k);
                                        if ~isempty(Pos)
                                            for PosL=1:length(Pos)
                                                Hubs=[Hubs;[Hubs(Pos(PosL),1),min(k,PsGrpNb1),max(k,PsGrpNb1),0,0]];
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            EndPos=size(Hubs,1);
                            if EndPos>0
                                FoundPos=[];
                                for m=[2,3]
                                    Pos=find(Hubs(1:EndPos,m)==PsGrpNb1);
                                    if ~isempty(Pos)
                                        for PosL=1:length(Pos)
                                            if isempty(find(FoundPos==Hubs(Pos(PosL),1)))
                                                FoundPos=[FoundPos,Hubs(Pos(PosL),1)];
                                                Hubs=[Hubs;[PsGrpNb1,min(Hubs(Pos(PosL),1),PsGrpNb2),max(Hubs(Pos(PosL),1),PsGrpNb2),0,0]];
                                            end
                                        end
                                    end
                                end
                            else
                                Hubs=[PsGrpNb1,min(PsGrpNb1,PsGrpNb2),max(PsGrpNb1,PsGrpNb2),0,0];
                            end
                        end
                        break
                    end
                end
            end
        end
        %add grp sizes
        if ~isempty(Hubs)
            for HubL=1:size(Hubs,1)
                Hubs(HubL,4)=length(PsGrp{Hubs(HubL,2)});
                Hubs(HubL,5)=length(PsGrp{Hubs(HubL,3)});
            end
        end
        %some ps groups might be empty due to splitting process
        Changed=[];
        IsCleared=0;
        MemPsGrp2=PsGrp;
        for i=length(PsGrp):-1:1
            if isempty(PsGrp{i})
                %clear empty Ps Grp
                PsGrp(i)=[];
                IsCleared=1;
            else
                Changed=[Changed;i];
            end
        end
        Changed=[Changed,[length(Changed):-1:1]'];
        if ~isempty(Hubs)&IsCleared
            %update indices of grps          
            for HubL=1:size(Hubs,1)
                for i=1:3
                    Pos=find(Changed(:,1)==Hubs(HubL,i));                    
                    Hubs(HubL,i)=Changed(Pos,2);
                end
            end
        end
    end

    %control absence of doublon
    PsGrpNb1=length(PsGrp);
    for GrpL=1:PsGrpNb1
        Control=[Control,PsGrp{GrpL}];
    end
    if length(Control)~=length(unique(Control))
        h=errordlg('Current ps groups are not a partition');
        waitfor(h)
        error('process canceled')
    end
end

% eventual extra probe set can't form triangle
ToGrp=setdiff(1:PsNb,Control);
if ~isempty(ToGrp)
    %search couples in the remaining points
    Pairs=[];
    if length(ToGrp)>1
        for PsL1=1:length(ToGrp)-1
            for PsL2=PsL1+1:length(ToGrp)
                if Paired(ToGrp(PsL1),ToGrp(PsL2))
                    Pairs=[Pairs;[ToGrp(PsL1),ToGrp(PsL2)]];
                    Control=[Control,ToGrp(PsL1),ToGrp(PsL2)];
                end
            end
        end
    end %of grouping pairs

    %split overlapping ps groups
    if size(Pairs,1)>1
        %keep unique pairs
        for PairL=size(Pairs,1):-1:1
            if length(find(Pairs==Pairs(PairL,1)))==1&length(find(Pairs==Pairs(PairL,2)))==1
                PsGrp{end+1,1}=Pairs(PairL,:);
                Pairs(PairL,:)=[];
            end
        end
        PsGrpNb1=length(PsGrp);
        %proceed hubs
        if ~isempty(Pairs)
            %all remaining pairs will be splitted because they contains hubs
            UsedPs=unique(Pairs);
            %forces to a raw vector format
            if size(UsedPs,2)==1
                UsedPs=UsedPs';
            end
            %recover hubs
            HubPs=[];
            for PsL=1:length(UsedPs)
                %recover all ps paired with the current ps
                CurrPaired=[];
                for PairL=1:size(Pairs,1)
                    if ~isempty(find(Pairs(PairL,:)==UsedPs(PsL)))
                        CurrPaired=[CurrPaired,Pairs(PairL,:)];
                    end
                end
                CurrPaired=setdiff(unique(CurrPaired),UsedPs(PsL));
                if length(CurrPaired)>1
                    HubPs=[HubPs,UsedPs(PsL)];
                end
            end
            NotHubPs=setdiff(UsedPs,HubPs);
            %put hubs at the end
            UsedPs=[NotHubPs,HubPs];
            for PsL=1:length(UsedPs)
                PsGrp{end+1,1}=UsedPs(PsL);
                if ~isempty(HubPs==UsedPs(PsL))            
                    CurrPaired=[];
                    for PairL=1:size(Pairs,1)
                        if ~isempty(find(Pairs(PairL,:)==UsedPs(PsL)))
                            CurrPaired=[CurrPaired,Pairs(PairL,:)];
                        end
                    end
                    CurrPaired=setdiff(unique(CurrPaired),UsedPs(PsL));
                    for PsL1=1:length(CurrPaired)-1
                        for PsL2=PsL1+1:length(CurrPaired)
                            Hubs=[Hubs;[length(PsGrp),PsGrpNb1+find(UsedPs==CurrPaired(PsL1)),PsGrpNb1+find(UsedPs==CurrPaired(PsL2)),1,1]];
                        end
                    end
                end
            end
        end
    elseif size(Pairs,1)==1
        PsGrp{end+1,1}=Pairs;
    end
end

%group orphan probe sets
ToGrp=setdiff(1:PsNb,Control);
if ~isempty(ToGrp)
    for PsL=1:length(ToGrp)
        PsGrp{end+1,1}=ToGrp(PsL);
    end %of grouping orphan
end

if ~isempty(PsGrp)
    GrpNb=length(PsGrp);
    GrpSize=zeros(GrpNb,1);
    for GrpL=1:GrpNb
        GrpSize(GrpL)=length(PsGrp{GrpL});
    end
    if ~isempty(Hubs)       
        %reorder grp
        for HubL=1:size(Hubs,1)
            if Hubs(HubL,3)<Hubs(HubL,2)
                Hubs(HubL,[2,3])=sort(Hubs(HubL,[2,3]));
                Hubs(HubL,[4,5])=fliplr(Hubs(HubL,[4,5]));

            end
        end
    end

    %recover statistiques on ps groups

    if nargout==5
        for GrpL=1:GrpNb
            %group size
            CurrSize=GrpSize(GrpL);
            if length(GrpSizes)>=CurrSize
                GrpSizes(CurrSize)=GrpSizes(CurrSize)+1;              
            else
                GrpSizes(CurrSize)=1;
            end

            if CurrSize>1
                % statistics on all the possible pairs in the current group
                %[number of pairs that don't pass the test (corr>0 in a nb of network<limit),...
                % number of pairs that pass the test (corr>0 in a nb of network>=limit),...
                % number of pairs that are significative in all the network,...
                % group size]
                CurrLinkType=zeros(1,4);
                CurrLinkType(4)=CurrSize;
                BadNb=0;
                for i=1:CurrSize-1
                    for j=i+1:CurrSize
                        Node1=PsGrp{GrpL}(i);
                        Node2=PsGrp{GrpL}(j);
                        CurrType=PairedType(Node1,Node2);
                        CurrLinkType(CurrType+1)=CurrLinkType(CurrType+1)+1;
                        if CurrType==0
                            BadNb=BadNb+1;
                            Pos=(Node1*(2*PsNb-Node1-1)/2)+Node2-PsNb;
                            BadLink=[BadLink;[Node1,Node2,0,CurrSize*(CurrSize-1)/2,length(find(PairedMat{Pos}==1)),length(find(PairedMat{Pos}==2))]];
                        end
                    end
                end
                if BadNb>0
                    for BadL=1:BadNb
                        BadLink(BadL,3)=BadNb;
                    end
                end
                LinkTypes=[LinkTypes;CurrLinkType];
            end
        end
    end

    %convert ps groups of Ps pos to ps groups of Ps ranks
    for i=1:GrpNb
        PsGrp{i}=PsRanks(PsGrp{i});
    end
end

if nargout==5
    varargout{1}=GrpSizes;
    varargout{2}=LinkTypes;
    varargout{3}=BadLink;
    varargout{4}=Hubs;
end

%% MERGE_PSGROUPS
function [PsGrp,PsGrpPos]=MERGE_PSGROUPS(PsGrp,PsGrpPos)
PsgrpNb=length(PsGrp);
ClearedGrp=[];
%merge ps groups which have a common edge
for GrpL1=1:PsgrpNb-1
    for GrpL2=GrpL1+1:PsgrpNb
        %tested PsGrp may have been merged in a previous step
        %and could be empty
        if ~isempty(PsGrp{GrpL2})
            if length(intersect(PsGrp{GrpL1},PsGrp{GrpL2}))>=2
                PsGrp{GrpL1}=unique([PsGrp{GrpL1},PsGrp{GrpL2}]);
                PsGrp{GrpL2}=[];
                ClearedGrp=[ClearedGrp;GrpL2];
            end
        end
        if PsGrpPos==GrpL2
            %the current ps group containing triangle(TriL1) has been merged in a
            %prevous ps group => update the current ps group
            PsGrpPos=GrpL1;
        end
    end
end
if ~isempty(ClearedGrp)
    %clear merged ps groups
    PsGrp(ClearedGrp)=[];
    %correct the current PsGrpPos
    GrpRanks=1:PsgrpNb;
    GrpRanks(ClearedGrp)=[];
    PsGrpPos=find(GrpRanks==PsGrpPos);
end
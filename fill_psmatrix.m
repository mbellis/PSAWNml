%========================%
% FUNCTION FILL_PSMATRIX %
%========================%
%
%FILL_PSMATRIX constructs PsMatrix which summarize information on probe sets in a 
% numeric table
%
%PsMatrix columns:
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


function [PsMatrix,PsMade]=fill_psmatrix(ClassRank,ProbeNb,NewPs,PsMatrix,PsMade,CurrPsRanks,Stat,GenePos,PsPos)

for PsL1=1:length(CurrPsRanks)
    CurrPsRank=CurrPsRanks(PsL1);
    if PsMade(CurrPsRank)==0
        PsMade(CurrPsRank)=1;

        % position of the gene targeted with the maximal number of probes
        CurrGenePos=GenePos(PsL1);       
        
        %other ps considered as targeting the same gene
        if ~isempty(PsPos)
            UsedPsPos=PsPos{PsL1};
        else
            UsedPsPos=1:length(NewPs{CurrPsRank}.psRanks{CurrGenePos});
        end
        
        %1st position: rank of the assigned gene to the current probe set
        PsMatrix(CurrPsRank,1)=NewPs{CurrPsRank}.geneRanks(CurrGenePos);

        %2nd position: position of the assigned gene in NewPs.geneNames
        PsMatrix(CurrPsRank,2)=CurrGenePos;

        %3rd position: target type of assigned gene
        PsMatrix(CurrPsRank,3)=NewPs{CurrPsRank}.target(CurrGenePos);

        %4th position: source type of assigned gene
        PsMatrix(CurrPsRank,4)=NewPs{CurrPsRank}.source(CurrGenePos);

        %5th positon: nb of probe targetting the assigned gene
        PsMatrix(CurrPsRank,5)=NewPs{CurrPsRank}.probeNb(CurrGenePos);

        %6th position : nb of not assigned genes targetted with the same nb of probes
        ProbeNbs=NewPs{CurrPsRank}.probeNb;
        IdemPos=find(ProbeNbs==PsMatrix(CurrPsRank,5));
        IdemNb=length(IdemPos);
        PsMatrix(CurrPsRank,6)=IdemNb-1;        
        if ~isempty(ProbeNbs)           
            %17th position and beyond: nb of  other genes targetted with all possible nb of
            %probes
            for NbL=1:length(ProbeNbs)            
                PsMatrix(CurrPsRank,17+min(ProbeNbs(NbL),ProbeNb))=PsMatrix(CurrPsRank,17+min(ProbeNbs(NbL),ProbeNb))+1;
            end
        end
        ProbeNbs(IdemPos)=[];
        %7th position : nb of not assigned genes targetted with less nb of probes
        PsMatrix(CurrPsRank,7)=length(ProbeNbs);
        %Ranks of all probe sets that are assigned to the current gene
        AllPsRanks=NewPs{CurrPsRank}.psRanks{CurrGenePos}(UsedPsPos);
        %idem minus the current probe set rank
        OtherPsRanks=setdiff(AllPsRanks,CurrPsRank);
        
        CurrGeneRank=NewPs{CurrPsRank}.geneRanks(CurrGenePos);

        %9th position: rank of the parent probe set
        PsMatrix(CurrPsRank,9)=CurrPsRank;
        
        %Process groups of transcript(s)
        Grps=NewPs{CurrPsRank}.grp{CurrGenePos};
        %all grps are identical and in the same order in all the PsMatrix
        %assign Got
        Got=zeros(length(AllPsRanks),2);
        GotRank=0;
        for GrpL=1:length(Grps)
            for PsL2=1:length(Grps{GrpL})
                if PsL2==1
                    GotRank=GotRank+1;
                end
                CurrPsPos=find(AllPsRanks==Grps{GrpL}(PsL2));
                Got(CurrPsPos,1)=GotRank;
            end
        end
        if length(AllPsRanks)>1
            if~isempty(Stat.hubs)
                for GrpL=1:length(Grps)
                    for PsL2=1:length(Grps{GrpL})
                        %process only if the ps target the current gene
                        if ~isempty(find(AllPsRanks==Grps{GrpL}(PsL2)))
                            %recover position of the current hub
                            Pos=find(Stat.hubs(:,1)==Grps{GrpL}(PsL2)&Stat.hubs(:,2)==find(NewPs{Grps{GrpL}(PsL2)}.geneRanks==CurrGeneRank));
                            if ~isempty(Pos)
                                % 13th position: [0,1] indicates that the probe set is paired with a pivot                                
                                GotRanks=[];
                                for HubL=1:length(Pos)
                                    RecoverGot=1;
                                    for PsL3=1:length(Grps{Stat.hubs(Pos(HubL),4)})
                                        if ~isempty(find(AllPsRanks==Grps{Stat.hubs(Pos(HubL),4)}(PsL3)))
                                            PsMatrix(Grps{Stat.hubs(Pos(HubL),4)}(PsL3),13)=1;
                                            if RecoverGot                                    
                                                RecoverGot=0;
                                                CurrPsPos=find(AllPsRanks==Grps{Stat.hubs(Pos(HubL),4)}(PsL3));
                                                GotRanks(end+1)=Got(CurrPsPos,1);                                                
                                            end
                                        end
                                    end

                                    RecoverGot=1;
                                    for PsL3=1:length(Grps{Stat.hubs(Pos(HubL),5)})
                                        if ~isempty(find(AllPsRanks==Grps{Stat.hubs(Pos(HubL),5)}(PsL3)))
                                            PsMatrix(Grps{Stat.hubs(Pos(HubL),5)}(PsL3),13)=1;
                                            if RecoverGot                                    
                                                RecoverGot=0;
                                                CurrPsPos=find(AllPsRanks==Grps{Stat.hubs(Pos(HubL),5)}(PsL3));
                                                GotRanks(end+1)=Got(CurrPsPos,1);                                                
                                            end
                                        end
                                    end
                                end
                                if ~isempty(GotRanks)
                                    %12th position: [0,1] indicates if the current
                                    %probe set is a pivot
                                    PsMatrix(Grps{GrpL}(PsL2),12)=1;
                                    %construct binary string
                                    if ~isempty(GotRanks)
                                        BinStr=repmat('0',1,max(GotRanks));
                                        for GotL=1:length(GotRanks)
                                            BinStr(GotRanks(GotL))='1';
                                        end
                                        BinStr=fliplr(BinStr);
                                        CurrPsPos=find(AllPsRanks==Grps{GrpL}(PsL2));
                                        Got(CurrPsPos,2)=bin2dec(BinStr);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        %8th position: nb of groups of transcripts corresponding to the assigned gene
        PsMatrix(CurrPsRank,8)=max(Got(:,1));
        
        %fill Group of Transcripts
        %10th position: Rank of the group of transcripts targetted by the current
        %               probe set in the assigned gene
        %11th position: Rank of the group(s) of transcripts targetted by the current
        %     probe set, if it is a pivot
        for PsL2=1:length(AllPsRanks)
            PsMatrix(AllPsRanks(PsL2),10)=Got(PsL2,1);
            PsMatrix(AllPsRanks(PsL2),11)=Got(PsL2,2);
        end                        
        
        %14th position: nb of probe sets that do not target the assigned gene but target
        % a common gene with the current probe set
        %find all probe set that target a gene targeted by the current probe set
        FullPsRanks=[];
        for GeneL=1:length(NewPs{CurrPsRank}.geneNames)
            FullPsRanks=union(FullPsRanks,NewPs{CurrPsRank}.psRanks{GeneL});
        end
        DiffPsRanks=setdiff(FullPsRanks,AllPsRanks);
        PsMatrix(CurrPsRank,14)=length(DiffPsRanks);
        
        %15th position: nb of genes that are targetted by probe sets that do not target the
        %assigned gene
        if PsMatrix(CurrPsRank,14)>0
            GeneNb=0;
            for GeneL=1:length(NewPs{CurrPsRank}.geneNames)
                for PsL2=1:length(NewPs{CurrPsRank}.psRanks{GeneL})
                    if ~isempty(find(DiffPsRanks==NewPs{CurrPsRank}.psRanks{GeneL}(PsL2)))
                        GeneNb=GeneNb+1;
                        break
                    end
                end
            end
            PsMatrix(CurrPsRank,15)=GeneNb;
        end
        
        %16th position: nb of genes that are targetted by other probe set with a nb
        %               of probes higher than the number of probes of the current probe set that
        %               target the assigned gene
        PsMatrix(CurrPsRank,16)=length(find(NewPs{CurrPsRank}.probeNb>NewPs{CurrPsRank}.probeNb(CurrGenePos)));
        
        %17th position ClassRank
        PsMatrix(CurrPsRank,17)=ClassRank;

        %fill PsMatrix for other probe sets
        %mark as made all other probe sets targetting the same gene
        for PsL2=1:length(OtherPsRanks)
            CurrPsRank2=OtherPsRanks(PsL2);
            if PsMade(CurrPsRank2)==0
                PsMade(CurrPsRank2)=1;
                %1st position: rank of assigned gene tot the current probe set
                PsMatrix(CurrPsRank2,1)=PsMatrix(CurrPsRank,1);
                %2nd position: pos of assigned gene in NewPs
                LocalGenePos=find(NewPs{OtherPsRanks(PsL2)}.geneRanks==CurrGeneRank);
                PsMatrix(CurrPsRank2,2)=LocalGenePos;
                %3rd position: target type of assigned gene
                PsMatrix(CurrPsRank2,3)=PsMatrix(CurrPsRank,3);
                %4th position: source type of assigned gene
                PsMatrix(CurrPsRank2,4)=PsMatrix(CurrPsRank,4);              
                %5th positon: nb of probe targetting the assigned gene
                PsMatrix(CurrPsRank2,5)=NewPs{CurrPsRank2}.probeNb(LocalGenePos);
                %6th position : nb of not assigned genes targetted with the same nb (or greater nb) of probes
                ProbeNbs=NewPs{CurrPsRank2}.probeNb;
                IdemPos=find(ProbeNbs>=PsMatrix(CurrPsRank2,5));
                IdemNb=length(IdemPos);
                PsMatrix(CurrPsRank2,6)=IdemNb-1;                
                if ~isempty(ProbeNbs)                    
                    %18th position and beyond: nb of  other genes targetted with all possible nb of
                    %probes
                    for NbL=1:length(ProbeNbs)
                        PsMatrix(CurrPsRank2,17+min(ProbeNbs(NbL),ProbeNb))=PsMatrix(CurrPsRank2,17+min(ProbeNbs(NbL),ProbeNb))+1;
                    end
                end
                ProbeNbs(IdemPos)=[];
                %7th position : nb of not assigned genes targetted with less nb of probes
                PsMatrix(CurrPsRank2,7)=length(ProbeNbs);
                
                %8th position: nb of groups of transcripts corresponding to the assigned
                %gene
                PsMatrix(CurrPsRank2,8)=PsMatrix(CurrPsRank,8);
                %9th position: rank of the parent probe set
                PsMatrix(CurrPsRank2,9)=PsMatrix(CurrPsRank,9);

                %14th position: nb of probe sets that do not target the assigned gene but target
                % a common gene with the current probe set
                %find all probe set that target a gene targeted by the current probe set
                FullPsRanks=[];
                for GeneL=1:length(NewPs{CurrPsRank2}.geneNames)
                    FullPsRanks=union(FullPsRanks,NewPs{CurrPsRank2}.psRanks{GeneL});
                end
                DiffPsRanks=setdiff(FullPsRanks,AllPsRanks);
                PsMatrix(CurrPsRank2,14)=length(DiffPsRanks);

                %15th position: nb of genes that are targetted by probe sets that do not target the
                %assigned gene
                if PsMatrix(CurrPsRank2,14)>0
                    GeneNb=0;
                    for GeneL=1:length(NewPs{CurrPsRank2}.geneNames)
                        for PsL3=1:length(NewPs{CurrPsRank2}.psRanks{GeneL})
                            if ~isempty(find(DiffPsRanks==NewPs{CurrPsRank2}.psRanks{GeneL}(PsL3)))
                                GeneNb=GeneNb+1;
                                break
                            end
                        end
                    end

                    PsMatrix(CurrPsRank2,15)=GeneNb;
                end

                %16th position: nb of genes that are targetted by other probe set with a nb
                %               of probes higher than the number of probes of the current probe set that
                %               target the assigned gene
                PsMatrix(CurrPsRank2,16)=length(find(NewPs{CurrPsRank2}.probeNb>NewPs{CurrPsRank2}.probeNb(LocalGenePos)));

                %17th position ClassRank
                PsMatrix(CurrPsRank2,17)=ClassRank;           
            end
        end
        if length(find(PsMatrix(:,10)>PsMatrix(:,8)))>0
            h=errordlg('PsMatrix(x,10)>PsMatrix(x,8)');
            waitfor(h)
            error('process canceled')
        end
        if length(find(PsMatrix(:,12)>0&PsMatrix(:,11)==0))>0
            h=errordlg('PsMatrix(x,12)>0&PsMatrix(x,11)==0');
            waitfor(h)
            error('process canceled')
        end       

    end % if PsMade(PsL1)
end %for PsL1


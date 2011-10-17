%==================
% FUNCTION PS2NET %
%==================

% PS2NET constructs the corespondance between probe set order in probeset file
% and probe set order in networks (CVM)

%INPUT PARAMETERS
% 1 ChipRank : chip rank

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

function ps2net(ChipRank)
global K
cd(K.dir.rawdata)
PsOrder=textread(sprintf('m%u_probeset.txt',ChipRank),'%s');
NetOrder=textread(sprintf('m%u_net.txt',ChipRank),'%s');
NotFound=[];
Ps2Net=zeros(size(PsOrder));
for PsL=1:length(PsOrder)
    NetPos=strmatch(PsOrder{PsL},NetOrder,'exact');
    if~isempty(NetPos)
        Ps2Net(PsL)=NetPos;
    else
        NotFound=[NotFound;PsL];
    end
end
if~isempty(NotFound)
    h=warndlg(sprintf('%u probe sets do not exist in network',length(NotFound)));
    waitfor(h)
    PsOrder(NotFound)
end
ChipPos=strmatch(sprintf('m%u',ChipRank),K.chip.myName,'exact');
SaveDir=fullfile(K.dir.mldata,K.chip.species{ChipPos},sprintf('m%u',ChipRank));
mkdir(SaveDir)
cd(SaveDir)
eval(sprintf('save m%u_ps2net Ps2Net',ChipRank)) 
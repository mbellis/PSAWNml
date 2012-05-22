%==================
% FUNCTION DEMO_PS
%==================
% 
%DEMO runs a demonstration
% demo must be run from inside the main directory:
% 'cd .../psawnML/'
% 'demo_ps'
%
%GLOBAL VARIABLES
% K.dir contains directory paths
% K.chip contains information on chips


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

global K

%control that demo is run from the right directory
ProgDir=pwd;
cd(ProgDir)
DirContent=dir;
Ok=0;
for i=1:length(DirContent)
    if isequal(DirContent(i).name,'import_targetnb.m')
        Ok=1;
        break
    end
end
if Ok==0
    h=errordlg('you must run demo from inside psawnML directory');
    waitfor(h)
    error('process canceled')
end

%inform directory names
K.dir.mlprog=pwd;
%if global_ps does not exist, directory names are informed automatically
%assuming that program has been installed from zipped file, and that data folders
%has been created properly under the main directory (e.g. psawn)
%if data have been installed elsewhere, global_ps must be loaded tochange K.dir where
%necessary

if ~exist('global_ps.mat','file')

    ProgPos=findstr('prog',K.dir.mlprog);
    if ~isempty(ProgPos)
        RootDir=ProgDir(1:ProgPos-2);
        K.dir.mldata=fullfile(RootDir,'data','mldata');
        K.dir.pydata=fullfile(RootDir,'data','pydata');
        K.dir.rawdata=fullfile(RootDir,'data','rawdata');
        K.dir.net=fullfile(RootDir,'data','net');
    else
        h=errordlg('create variable global_ps containing K.dir.mldata, K.dir.pydata, K.dir.rawdata and K.dir.net, save it as global_ps.mat, and start again demo_ps');
        waitfor(h)
    end
    cd(K.dir.mlprog)
    IsK=0;
else
    load global_ps
    MemK=K;
    IsK=1;
end

%create structures
cd(K.dir.rawdata)
[K.chip.myName,K.chip.name,K.chip.shortName,K.chip.species,K.chip.probesetNb,K.chip.probeNb,K.chip.compName,K.chip.ens47Name,K.chip.ens48Name,K.chip.geoName]=textread('chip.txt','%s%s%s%s%u%u%s%s%s%s','delimiter','\t');
[K.chip.ref.rank,K.chip.ref.signal]=textread('signal_vs_rank.txt','%f%f','delimiter','\t');
if IsK==0
    K.chip.maxProbeNb=zeros(length(K.chip.myName),1);
else
    K.chip.maxProbeNb=MemK.chip.maxProbeNb;    
    clear MemK
end

%save global variable
cd(K.dir.mlprog)
save global_ps K
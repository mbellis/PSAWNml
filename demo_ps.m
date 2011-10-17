%==================
% FUNCTION DEMO_PS
%==================
% 
% DEMO runs a demonstration
% demo must be run from inside the main directory:
% 'cd .../psawnML/'
% 'demo_ps'
%
% GLOBAL VARIABLES
% K.dir contains directory paths


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %                               
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

global K

%control that demo is run from the right directory
RootDir=pwd;
cd(RootDir)
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
K.dir.mldata='/home/mbellis/sosma/data/psawn/mldata';
K.dir.pydata='/home/mbellis/sosma/data/psawn/pydata';
K.dir.rawdata='/home/mbellis/sosma/data/psawn/rawdata';
K.dir.net='/home/mbellis/net';
K.dir.table='/usr/data/net';

cd(K.dir.mlprog)
IsK=0;
if exist('global_ps.mat','file')
    load global_ps
    MemK=K;
    IsK=1;
end

%create structures
cd(K.dir.rawdata)
[K.chip.myName,K.chip.name,K.chip.shortName,K.chip.species,K.chip.probesetNb,K.chip.probeNb,K.chip.compName,K.chip.ens47Name,K.chip.ens48Name,K.chip.geoName]=textread('chip.txt','%s%s%s%s%u%u%s%s%s%s','delimiter','\t');
[K.chip.ref.rank,K.chip.ref.signal]=textread('signal_vs_rank.txt','%f%f','delimiter','\t');
if IsK==0
    K.chip.probesetNb=zeros(length(K.chip.myName),1);
    K.chip.maxProbeNb=zeros(length(K.chip.myName),1);
else
    K.chip.probesetNb=MemK.chip.probesetNb;
    K.chip.maxProbeNb=MemK.chip.maxProbeNb;
    
    clear MemK
end

%save global variable
cd(K.dir.mlprog)
save global_ps K
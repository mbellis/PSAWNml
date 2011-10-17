%==========================
% Function HYPERGEOMETRIC %
%==========================

% HYPERGEOMETRIC calculates probability p(X<=k) if k is less than the expected number
% or p(X<=k) if k is more than the expected number where X is an hypergeometric
% random variable and k the number of elements of the extracted sample
% which have the property

%INPUT
% 1 InSampleNb: the number of elements of the sample which 
%               have the property
% 2 SampleSize: the number of elements which are randomly extracted
% 3   InPopuNb: the number of elements of the population which 
%               have the property
% 4   PopuSize: the total number of elements of the population

%OUTPUT
% 1 PVal: probability

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

function PVal=hypergeometric(InSampleNb,SampleSize,InPopuNb,PopuSize)
% put asside situation with p-value~=0 because O link in each


if InSampleNb>SampleSize
    InSampleNb=SampleSize;
end
if InSampleNb>InPopuNb
    InSampleNb=InPopuNb;
end

if InSampleNb>SampleSize*InPopuNb/PopuSize
    PVal = sum(hypgeo([InSampleNb:SampleSize],SampleSize,InPopuNb,PopuSize));
else
    PVal= sum(hypgeo([0:InSampleNb],SampleSize,InPopuNb,PopuSize));
end
PVal=log10(PVal);

if InSampleNb<SampleSize*InPopuNb/PopuSize
    PVal=-PVal;
end

 
function  p = hypgeo(InSampleNb,SampleSize,InPopuNb,PopuSize)
%INPUT
% 1 InSampleNb: the number of elements of the sample which 
%               have the studied property
% 2 SampleSize: the number of elements which are randomly extracted
% 3   InPopuNb: the number of elements of the population which 
%               have the studied property
% 4   PopuSize: the total number of elements of the population

%OUTPUT
% p: probability Prob(X=InSampleNb) where X is an hypergeometric
%    random variable with parameters SampleSize,InPopuNb,PopuSize

% Anders Holtsberg, 18-11-93
% Copyright (c) Anders Holtsberg

p = gammaln(InPopuNb+1) - gammaln(InSampleNb+1) - gammaln(InPopuNb-InSampleNb+1);
p = p + gammaln(PopuSize-InPopuNb+1) - gammaln(SampleSize-InSampleNb+1) - gammaln(PopuSize-InPopuNb-SampleSize+InSampleNb+1);
p = p - gammaln(PopuSize+1) + gammaln(SampleSize+1) + gammaln(PopuSize-SampleSize+1);
p = exp(p);


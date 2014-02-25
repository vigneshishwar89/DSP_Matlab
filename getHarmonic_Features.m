function R = getHarmonic_Features(queryFolder,ext)
addpath '/Data2/Masters_UPF/Code/'
addpath '/Data2/Data/Code_Genmax'
wavfiles=dir(fullfile(queryFolder,ext))
for i=1:length(wavfiles)
    tic
    actrName=[wavfiles(i).name '.allctrs'];
    pctrName=[wavfiles(i).name '.pitch'];
    barName=[wavfiles(i).name '.bars.txt'];
    barnvName=[wavfiles(i).name '.barsnv.txt'];
    actr=load(actrName);
    pctr=load(pctrName);
    pctr=pctr(:,2);
    R1 = getHarmWeights_aggr_mod(wavfiles(i).name,barName,barnvName,pctr,actr,hann(2047), 2048, 30, -80, 0.2)
    a=R1.VocInd;
    b=R1.NVocInd;
    c=R1.OConInd;
    vocFeat = R1.VocFeat;
    %vocFeat_locs = R1.VocFeatTime;
    nvocFeat = R1.NVocFeat;
    %nvocFeat_locs = R1.NVocFeatTime;
    otherFeat = R1.otherFeat;
     outFvoc = [wavfiles(i).name '.vocHfeat_aggr'];
%     outFvoc_locs= [wavfiles(i).name '.vocHfeat_locs'];
    outFnvoc = [wavfiles(i).name '.nvocHfeat_aggr'];
    outFother = [wavfiles(i).name '.otherFeat_aggr'];
    %outFnvoc_locs = [wavfiles(i).name '.nvocHfeat_locs_1ctr'];
    dlmwrite(outFvoc,vocFeat,'delimiter','\t')
%     dlmwrite(outFvoc_locs,vocFeat_locs,'delimiter','\t')
    dlmwrite(outFnvoc,nvocFeat,'delimiter','\t')
    dlmwrite(outFother,otherFeat,'delimiter','\t')
    %dlmwrite(outFnvoc_locs,nvocFeat_locs,'delimiter','\t')
    toc
end
R.V=a;
R.NV=b;
R.O=c;
R.Ratio = R1.Ratio;
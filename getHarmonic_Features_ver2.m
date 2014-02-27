function R = getHarmonic_Features(queryFolder,ext)
addpath '/Data2/Masters_UPF/Code/'
addpath '/Data2/Data/Code_Genmax'
wavfiles=dir(fullfile(queryFolder,ext))
for i=1:length(wavfiles)
    tic
    actrName=[wavfiles(i).name '.allctrs'];
    pctrName=[wavfiles(i).name '.pitch']
    barName=[wavfiles(i).name '.bars.txt'];
    barnvName=[wavfiles(i).name '.barsnv.txt'];
    actr=load(actrName);
    pctr=load(pctrName);
    pctr=pctr(:,2);
    R1 = getHarmWeights_aggr_mod(wavfiles(i).name,barName, barnvName, pctr,actr,hann(2047), 2048, 30, -80, 0.2)
    vocFeat = R1.VocFeat;
    harm2Feat = R1.harm2Feat;
    harm3Feat = R1.harm3Feat;
    subharmFeat = R1.subharmFeat;
    viogrFeat = R1.viogrFeat;
    vioctrFeat = R1.vioctrFeat;
    
    
    outFvoc = [wavfiles(i).name '.vocHfeat_aggr'];
    outFharm2 = [wavfiles(i).name '.harm2_aggr'];
    outFharm3 = [wavfiles(i).name '.harm3_aggr'];
    outFsubharm = [wavfiles(i).name '.subharm_aggr'];
    outFviogr = [wavfiles(i).name '.viogr_aggr'];
    outFvioctr = [wavfiles(i).name '.vioctr_aggr'];
    
    dlmwrite(outFvoc,vocFeat,'delimiter','\t')
    dlmwrite(outFharm2,harm2Feat,'delimiter','\t')
    dlmwrite(outFharm3,harm3Feat,'delimiter','\t')
    dlmwrite(outFsubharm,subharmFeat,'delimiter','\t')
    dlmwrite(outFviogr,viogrFeat,'delimiter','\t')
    dlmwrite(outFvioctr,vioctrFeat,'delimiter','\t')
    toc
end
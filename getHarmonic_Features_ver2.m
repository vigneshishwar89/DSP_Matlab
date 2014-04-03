function getHarmonic_Features_ver2(queryFolder,ext)
addpath '/Data2/Masters_UPF/Code/'
addpath '/Data2/Data/Code_Genmax'
wavfiles=dir(fullfile(queryFolder,ext))
for i=1:length(wavfiles)
    %vocalslp = []; harm2slp = []; harm3slp = []; subharmslp = []; viogrslp = []; violrslp = []; vioctrslp = [];
    tic
    actrName=[wavfiles(i).name '.allctrs'];
    pctrName=[wavfiles(i).name '.pitch']
    barName=[wavfiles(i).name '.bars_corr.txt'];
    %barnvName=[wavfiles(i).name '.barsnv.txt'];
    actr=load(actrName);
    pctr=load(pctrName);
    pctr=pctr(:,2);
    
    R1 = getHarmWeights_aggr_mod_ver3(wavfiles(i).name,barName,pctr,actr,hann(2047), 2048, 30, -80, 0.2,30);
    
    outfName = [wavfiles(i).name(1:end-4) '.mat'];
    save(outfName, '-struct', 'R1');
    toc
end
disp('Finally over!!!!')


function getHarmonic_Features_ver2(queryFolder,ext)
addpath '/Data2/Masters_UPF/Code/'
addpath '/Data2/Data/Code_Genmax'
wavfiles=dir(fullfile(queryFolder,ext))
for i=1:length(wavfiles)
    vocalslp = []; harm2slp = []; harm3slp = []; subharmslp = []; viogrslp = []; violrslp = []; vioctrslp = [];
    tic
    actrName=[wavfiles(i).name '.allctrs'];
    pctrName=[wavfiles(i).name '.pitch']
    barName=[wavfiles(i).name '.bars.txt'];
    %barnvName=[wavfiles(i).name '.barsnv.txt'];
    actr=load(actrName);
    pctr=load(pctrName);
    pctr=pctr(:,2);
    R1 = getHarmWeights_aggr_mod_ver3(wavfiles(i).name,barName,pctr,actr,hann(2047), 2048, 30, -80, 0.2,20)
    vocFeat = R1.VocFeat;
    vocSlopes = R1.vocFrameSlopes;
    harm2Feat = R1.harm2Feat;
    harm2Slopes = R1.harm2FrameSlopes;    
    harm3Feat = R1.harm3Feat;
    harm3Slopes = R1.harm3FrameSlopes;    
    subharmFeat = R1.subharmFeat;
    subharmSlopes = R1.subharmFrameSlopes;    
    viogrFeat = R1.viogrFeat;
    viogrSlopes = R1.viogrFrameSlopes;    
    vioctrFeat = R1.vioctrFeat;
    vioctrSlopes = R1.vioctrFrameSlopes;    
    violrFeat = R1.violrFeat;
    violrSlopes = R1.violrFrameSlopes;   
    
    for voc = 1:length(vocSlopes)
    vocalslp = cat(1,vocalslp,vocSlopes(voc).frameslopes);
    end
    vocalslp(:,1) = vocalslp(:,1)*(128/44100);
    for h2 = 1:length(harm2Slopes)
        harm2slp = cat(1,harm2slp,harm2Slopes(h2).frameslopes);
    end
    harm2slp(:,1) = harm2slp(:,1)*(128/44100);
    for h3 = 1:length(harm3Slopes)
        harm3slp = cat(1,harm3slp,harm3Slopes(h3).frameslopes);
    end
    harm3slp(:,1) = harm3slp(:,1)*(128/44100);
    for sh = 1:length(subharmSlopes)
        subharmslp = cat(1,subharmslp,subharmSlopes(sh).frameslopes);
    end
    subharmslp(:,1) = subharmslp(:,1)*(128/44100);
    for vg = 1:length(viogrSlopes)
        viogrslp = cat(1,viogrslp,viogrSlopes(vg).frameslopes);
    end
    viogrslp(:,1) = viogrslp(:,1)*(128/44100);
    for vl = 1:length(violrSlopes)
        violrslp = cat(1,violrslp,violrSlopes(vl).frameslopes);
    end
    violrslp(:,1) = violrslp(:,1)*(128/44100);
    for vc = 1:length(vioctrSlopes)
        vioctrslp = cat(1,vioctrslp,vioctrSlopes(vc).frameslopes);
        vc
    end
    vioctrslp(:,1) = vioctrslp(:,1)*(128/44100);
    
    
    outFvoc = [wavfiles(i).name '.vocHfeat_aggr'];
    outFvocslp = [wavfiles(i).name '.vocHfeat_fslp'];
    outFharm2 = [wavfiles(i).name '.harm2_aggr'];
    outFharm2slp = [wavfiles(i).name '.harm2_fslp'];
    outFharm3 = [wavfiles(i).name '.harm3_aggr'];
    outFharm3slp = [wavfiles(i).name '.harm3_fslp'];
    outFsubharm = [wavfiles(i).name '.subharm_aggr'];
    outFsubharmslp = [wavfiles(i).name '.subharm_fslp'];
    outFviogr = [wavfiles(i).name '.viogr_aggr'];
    outFviogrslp = [wavfiles(i).name '.viogr_fslp'];
    outFvioctr = [wavfiles(i).name '.vioctr_aggr'];
    outFvioctrslp = [wavfiles(i).name '.vioctr_fslp'];
    outFviolr = [wavfiles(i).name '.violr_aggr'];
    outFviolrslp = [wavfiles(i).name '.violr_fslp'];
    
    % write to file harmonic weight features
    dlmwrite(outFvoc,vocFeat,'delimiter','\t')
    dlmwrite(outFharm2,harm2Feat,'delimiter','\t')
    dlmwrite(outFharm3,harm3Feat,'delimiter','\t')
    dlmwrite(outFsubharm,subharmFeat,'delimiter','\t')
    dlmwrite(outFviogr,viogrFeat,'delimiter','\t')
    dlmwrite(outFvioctr,vioctrFeat,'delimiter','\t')
    dlmwrite(outFviolr,violrFeat,'delimiter','\t')
    % write to file  slopes
    dlmwrite(outFvocslp,vocalslp,'delimiter','\t')
    dlmwrite(outFharm2slp,harm2slp,'delimiter','\t')
    dlmwrite(outFharm3slp,harm3slp,'delimiter','\t')
    dlmwrite(outFsubharmslp,subharmslp,'delimiter','\t')
    dlmwrite(outFviogrslp,viogrslp,'delimiter','\t')
    dlmwrite(outFvioctrslp,vioctrslp,'delimiter','\t')
    dlmwrite(outFviolrslp,violrslp,'delimiter','\t')
    toc
end
disp('Finally over!!!!')


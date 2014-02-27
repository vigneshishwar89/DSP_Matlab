function R = getHarmWeights_aggr_mod(wavFile,barFilev,barFilenv,pctr,actr,win, NFFT, nHarm, thsld, maxhd)
% Gives Vocal and NON Vocal sections based on vocal and non-vocal
% annotations. Other contours that are not annotated are contours which are
% either vocal, violin or harmonic errors. They have to be filtered. Also
% gives the ratio of the contours with the original contour.
[wave1,Fs]=wavread(wavFile);
%Fs = 44100;
hfeat_voc=[];
hfeat_nvoc=[];
hfeat_locs_voc=[];
hfeat_locs_nvoc=[];
%aggrMag_Voc=[];
%aggrMag_NVoc=[];
[m n]=size(actr);
hop = 128;
cnt1 = 1;
cnt2 = 1;
cnt3 = 1;
cnt4 = 1;
cnt5 = 1;
cnt6 = 1;
bars=load(barFilev);
barsnv = load(barFilenv);
pbuff1=zeros(length(pctr),1);
pbuff = zeros(length(pctr),1);
[g h]=size(bars);
[u v]=size(barsnv);
for i = 1:g
    stbars = bars(i,1);
    durbars = bars(i,3);
    endbars = stbars+durbars;
    stbarSamp = round((stbars*Fs)/hop);
    endbarSamp = round((endbars*Fs)/hop);
    if stbarSamp == 0
        stbarSamp = stbarSamp +1;
    end
    pbuff(stbarSamp:endbarSamp)=1;
end

for l=1:m
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the All Contours File(each line is a contour) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ac1=actr(l,:);
    ac1=ac1(ac1~=inf);
    ac2=ac1(3:end);
    st=ac1(1,2);
    stSamp=round((st*Fs)/hop);
    if (stSamp == 0)
        stSamp = stSamp + 1;
    end
    lenCtr=length(ac2);
    endSamp=stSamp+lenCtr-1;
    xCtr=[stSamp:endSamp];
    xCtr1=xCtr+1;
    stSamp_wf=xCtr1(1)*128;
    endSamp_wf=xCtr1(end)*128;
    wave1_clip=wave1(stSamp_wf:endSamp_wf);
    yCtr=55*(2.^(((10.*(ac2)))./1200));
    ppYctr = pctr(xCtr1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Checking with GT and getting Voiced getHarmWeights.mand Not Voice Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    err1=sum(abs(yCtr-ppYctr'))/length(yCtr);
    flag = ~sum(pbuff(stSamp:endSamp)~=1);
    if flag
        if (err1<0.01)
            vocY=yCtr';
            vocX=xCtr1';
            vocPitch = [vocX vocY];
            pitcharr=vocPitch;
            [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
            hmags = hmags';
            linMag = 10.^(hmags./20);
            temp = sum(linMag);
            aggrMag_Voc(cnt1,:) = temp./sum(temp);
            vocconInd(cnt1) = l;
            cnt1 = cnt1 + 1;
            else
                doublePth = ppYctr.*2;
                compPth = yCtr;
                buff1 = zeros(length(ppYctr),1);
                buff2 = zeros(length(yCtr),1);
                buff1(ppYctr~=0)=1;
                buff2(yCtr~=0)=1;
                com = and(buff1,buff2);
                comInd = find(com);
                arr1 = doublePth(comInd);
                arr2 = yCtr(comInd);
                d1 = sum(abs(arr1-arr2'))/length(arr1);
                if (d1<4)
                    harm2X = xCtr1';
                    harm2Y = yCtr';
                    harm2Pitch = [harm2X harm2Y];
                    pitcharr=harm2Pitch;
                    [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp = sum(linMag);
                    aggrMag_harm2(cnt2,:) = temp./sum(temp);
                    harm2Ind(cnt2) = l;
                    cnt2 = cnt2 + 1;
                else
                    triplePth = ppYctr.*3;
                    compPth = yCtr;
                    buff1 = zeros(length(ppYctr),1);
                    buff2 = zeros(length(yCtr),1);
                    buff1(ppYctr~=0)=1;
                    buff2(yCtr~=0)=1;
                    com = and(buff1,buff2);
                    comInd = find(com);
                    arr1 = triplePth(comInd);
                    arr2 = yCtr(comInd);
                    d1 = sum(abs(arr1-arr2'))/length(arr1);
                    if (d1<4)
                        harm3X = xCtr1';
                        harm3Y = yCtr';
                        harm3Pitch = [harm3X harm3Y];
                        pitcharr=harm3Pitch;
                        [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                        hmags = hmags';
                        linMag = 10.^(hmags./20);
                        temp = sum(linMag);
                        aggrMag_harm3(cnt3,:) = temp./sum(temp);
                        harm3Ind(cnt3) = l;
                        cnt3 = cnt3 + 1;
                    else
                        halfPth = ppYctr.*0.5;
                        compPth = yCtr;
                        buff1 = zeros(length(ppYctr),1);
                        buff2 = zeros(length(yCtr),1);
                        buff1(ppYctr~=0)=1;
                        buff2(yCtr~=0)=1;
                        com = and(buff1,buff2);
                        comInd = find(com);
                        arr1 = halfPth(comInd);
                        arr2 = yCtr(comInd);
                        d1 = sum(abs(arr1-arr2'))/length(arr1);
                        if (d1<4)
                            subharmX = xCtr1';
                            subharmY = yCtr';
                            subharmPitch = [subharmX subharmY];
                            pitcharr=subharmPitch;
                            [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                            hmags = hmags';
                            linMag = 10.^(hmags./20);
                            temp = sum(linMag);
                            aggrMag_subharm(cnt4,:) = temp./sum(temp);
                            harm3Ind(cnt4) = l;
                            cnt4 = cnt4 + 1;                        
                        elseif (median(yCtr)>median(ppYctr))
                            vioCtrGrRgX = xCtr1;
                            vioCtrGrRgY = yCtr;
                            vioGrPitch = [vioCtrGrRgX vioCtrGrRgY]; 
                            pitcharr = vioGrPitch;
                            [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                            hmags = hmags';
                            linMag = 10.^(hmags./20);
                            temp = sum(linMag);
                            aggrMag_viogr(cnt5,:) = temp./sum(temp);
                            viogrInd(cnt5) = l;
                            cnt5 = cnt5 + 1;               
                        else
                            vioCtrX = xCtr1;
                            vioCtrY = yCtr;
                            vioctr = [vioCtrX vioCtrY];
                            pitcharr = vioctr
                            [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                            hmags = hmags';
                            linMag = 10.^(hmags./20);
                            temp = sum(linMag);
                            aggrMag_vioctr(cnt6,:) = temp./sum(temp);
                            vioctrInd(cnt6) = l;
                            cnt6 = cnt6 + 1;               
                        end
                    end
                end
        end
    end
end
    R.VocFeat = aggrMag_Voc;
    %R.VocFeatTime=hfeat_locs_voc;
    if exist('aggrMag_harm2')==1
        R.harm2Feat = aggrMag_harm2;
    else
        R.harm2Feat = 0;
    end
    if exist('aggrMag_harm3')
        R.harm3Feat = aggrMag_harm3;
    else
        R.harm3Feat = 0;
    end
    if exist('aggrMag_subharm')
        R.subharmFeat = aggrMag_subharm;
    else
        R.subharmFeat = 0;
    end
    if exist('aggrMag_viogr')
        R.viogrFeat = aggrMag_viogr;
    else
        R.viogrFeat = 0;
    end
    if exist('aggrMag_vioctr')
        R.vioctrFeat = aggrMag_vioctr;
    else
        R.vioctrFeat = 0;
    end
    
    
    
    
%     R.VocInd = vocconInd;
%     if exist('nvocconInd')
%         R.NVocInd = nvocconInd;
%     else
%         R.NVocInd = NaN;
%     end
%     if exist('ocontInd')
%         R.OConInd = ocontInd;
%     else
%         R.OConInd = NaN;
%     end
%     R.Ratio = rat1
    %R.NVocFeatTime=hfeat_locs_nvoc;



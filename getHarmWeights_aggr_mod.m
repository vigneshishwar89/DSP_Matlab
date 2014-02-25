function R = getHarmWeights_aggr(wavFile,barFilev,barFilenv,pctr,actr,win, NFFT, nHarm, thsld, maxhd)
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
for i = 1:u
    stbarsnv = barsnv(i,1);
    durbarsnv = barsnv(i,3);
    endbarsnv = stbarsnv+durbarsnv;
    stbarSampnv = round((stbarsnv*Fs)/hop);
    endbarSampnv = round((endbarsnv*Fs)/hop);
    if stbarSampnv == 0
        stbarSampnv = stbarSampnv +1;
    end
    pbuff1(stbarSampnv:endbarSampnv)=1;
end

for j=1:m
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the All Contours File(each line is a contour) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ac1=actr(j,:);
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
    err=sum(abs(yCtr-ppYctr'))/length(yCtr);
    rat = ppYctr'./yCtr;
    rat1(j) = mean(round(rat));
    if (pbuff(stSamp:endSamp)==1 & err<0.01)
        vocY=yCtr';
        vocX=xCtr1';
        vocPitch = [vocX vocY];
        pitcharr=vocPitch;
%         pin=stSamp_wf;
%         pend=endSamp_wf;
        [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
        %if (~isempty(hmags))
        hmags = hmags';
        linMag = 10.^(hmags./20);
        temp = sum(linMag);
        aggrMag_Voc(cnt1,:) = temp./sum(temp);
        vocconInd(cnt1) = j;
        cnt1 = cnt1 + 1;
        %end    
%         hfeat_voc=cat(1,hfeat_voc,hmags');
%         hfeat_locs_voc=cat(1,hfeat_locs_voc,hmagsLoc');
    elseif(pbuff1(stSamp:endSamp)==1)
        nvocY=yCtr';
        nvocX=xCtr1';
        nvocPitch = [nvocX nvocY];
        pitcharr=nvocPitch;
%         pin=stSamp_wf(1);
%         pend=endSamp_wf(end);
        [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
        hmags = hmags';
        %if (~isempty(hmags))
        linMag = 10.^(hmags./20);
        temp = sum(linMag);
        aggrMag_NVoc(cnt2,:) = temp./sum(temp);
        nvocconInd(cnt2) = j;
        cnt2 = cnt2 + 1;
    else
        otherY=yCtr';
        otherX=xCtr1';
        otherPitch = [otherX otherY];
        pitcharr = otherPitch;
        [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
        hmags = hmags';
        %if (~isempty(hmags))
        linMag = 10.^(hmags./20);
        temp = sum(linMag);
        aggrMag_other(cnt3,:) = temp./sum(temp);
        ocontInd(cnt3) = j;
        cnt3 = cnt3 + 1;
        
    end
end
R.VocFeat = aggrMag_Voc;
%R.VocFeatTime=hfeat_locs_voc;
if exist('aggrMag_NVoc')==1
    R.NVocFeat = aggrMag_NVoc;
else
    R.NVocFeat = 0;
end
if exist('aggrMag_other')
    R.otherFeat = aggrMag_other;
else
    R.otherFeat = 0;
end
R.VocInd = vocconInd;
if exist('nvocconInd')
    R.NVocInd = nvocconInd;
else
    R.NVocInd = NaN;
end
if exist('ocontInd')
    R.OConInd = ocontInd;
else
    R.OConInd = NaN;
end
R.Ratio = rat1
%R.NVocFeatTime=hfeat_locs_nvoc;
    
    

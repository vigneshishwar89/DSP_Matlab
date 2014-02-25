function R = getHarmWeights(wavFile,pctr,actr,win, NFFT, nHarm, thsld, maxhd)
[wave1,Fs]=wavread(wavFile);
hfeat_voc=[];
hfeat_nvoc=[];
hfeat_locs_voc=[];
hfeat_locs_nvoc=[];
[m n]=size(actr);
hop = 128;
for j=1:m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the All Contours File(each line is a contour) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ac1=actr(j,:);
    ac1=ac1(ac1~=inf);
    ac2=ac1(3:end);
    st=ac1(1,2);
    stSamp=round((st*Fs)/hop);
    lenCtr=length(ac2);
    endSamp=stSamp+lenCtr-1;
    xCtr=[stSamp:endSamp];
    xCtr1=xCtr+1;
    stSamp_wf=xCtr1(1)*128;
    endSamp_wf=xCtr1(end)*128;
    yCtr=55*(2.^(((10.*(ac2)))./1200));
    ppYctr = pctr(xCtr1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Checking with GT and getting Voiced getHarmWeights.mand Not Voice Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    err=sum(abs(yCtr-ppYctr'))/length(yCtr);
    if (err<0.01)
        vocY=yCtr';
        vocX=xCtr1';
        vocPitch = [vocX vocY];
        pitcharr=vocPitch;
        pin=stSamp_wf;
        pend=endSamp_wf;
        [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wavFile, win, NFFT, nHarm, thsld, maxhd, pin,pend);
        hfeat_voc=cat(1,hfeat_voc,hmags');
        hfeat_locs_voc=cat(1,hfeat_locs_voc,hmagsLoc');
    elseif (err>100)
        nvocY=yCtr';
        nvocX=xCtr1';
        nvocPitch = [nvocX nvocY];
        pitcharr=nvocPitch;
        pin=stSamp_wf(1);
        pend=endSamp_wf(end);
        [hmags,hmagsLoc]=HarmonicSubtraction_New_v2(pitcharr, wavFile, win, NFFT, nHarm, thsld, maxhd, pin,pend);
        hfeat_nvoc=cat(1,hfeat_nvoc,hmags');
        hfeat_locs_nvoc=cat(1,hfeat_locs_nvoc,hmagsLoc');
    end
end
R.VocFeat = hfeat_voc;
R.VocFeatTime=hfeat_locs_voc;
R.NVocFeat = hfeat_nvoc;
R.NVocFeatTime=hfeat_locs_nvoc;
    
    

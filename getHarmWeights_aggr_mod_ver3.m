function R = getHarmWeights_aggr_mod_ver3(wavFile,barFilev,pctr,actr,win, NFFT, nHarm, thsld, maxhd,slpnum)
% Gives Vocal and NON Vocal sections based on vocal and non-vocal
% annotations. Other contours that are not annotated are contours which are
% either vocal, violin or harmonic errors. They have to be filtered. Also
% gives the ratio of the contours with the original contour.
% wavFile = waveFile barFilev = Barfile, pctr = pitchcontour, actr = all
% contour , win = window function, NFFT = NFFT size, nHarm = numHarmonics,
% thsld = threshold of peak picking of harmonics, maxhd: max. relative
% deviation in harmonic detection (ex: .2), slpnum = number of harmonics
% considered for slope calculation.

%% Read Wave
[wave1,Fs]=wavread(wavFile);

%% Initialize variables
hfeat_voc=[];
hfeat_nvoc=[];
hfeat_locs_voc=[];
hfeat_locs_nvoc=[];

idmat = -1*ones(5,length(pctr));
actrmat = -1*ones(5,length(pctr));

[m n]=size(actr);

hop = 128;

cnt1 = 1;
cnt2 = 1;
cnt3 = 1;
cnt4 = 1;
cnt5 = 1;
cnt6 = 1;
cnt7 = 1;
bars=load(barFilev);

pbuff1=zeros(length(pctr),1);
pbuff = zeros(length(pctr),1);
[g h]=size(bars);


%% Initialize buffer of size of pitch contour with ones at annotation locations

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

%% Process All contours to initialize idmatrix with the contour number and actr matrix with the contour values 

% Assume 5 contours can be there at max at the same time
% If the annotation region is true and if one of the contours lies within a
% certain starttime and end time and another contour starts before this
% contoour, then the next line of idmat and actrmat gets filled with the id number and
% contour value. This is to analyse the contours in blocks of time since
% one contour does not pertain to one source and the predominant final
% output considers start time to be accurate and truncates the previous
% contour if the start of the next contour overlaps. 

for l=1:m    
    ac1=actr(l,:);
    ac1=ac1(ac1~=inf);
    ac2=ac1(3:end);
    st=ac1(1,2);  
    stSamp(l)=round((st*Fs)/hop); 
    if (stSamp(l) == 0)
        stSamp(l) = stSamp(l) + 1;
    end
    endSamp(l) = stSamp(l) + length(ac2) -1;
    xCtr = [stSamp(l):endSamp(l)]; xCtr = xCtr+1;
    yCtr = 55*(2.^(((10.*(ac2)))./1200));
    flag1 = ~sum(idmat(1,xCtr)~=-1);
    flag2 = ~sum(idmat(2,xCtr)~=-1);
    flag3 = ~sum(idmat(3,xCtr)~=-1);
    flag4 = ~sum(idmat(4,xCtr)~=-1);
    flag5 = ~sum(idmat(5,xCtr)~=-1);
    id = l;
    if flag1
        idmat(1,xCtr) = id ;
        actrmat(1,xCtr) = yCtr;
    elseif flag2
        idmat(2,xCtr) = id;
        actrmat(2,xCtr) = yCtr;
    elseif flag3
        idmat(3,xCtr) = id ;
        actrmat(3,xCtr) = yCtr;
    elseif flag4
        idmat(4,xCtr) = id;
        actrmat(4,xCtr) = yCtr;
    elseif flag5
        idmat(5,xCtr) = id ;
        actrmat(5,xCtr) = yCtr;
    end
end

%% Decide based on some heuristics the classes of contours and calculate framewise slope

% 1) Vocal 2) 2nd harmonic 3) Third Harmonic 4) Subharmonic 5) Other contours in the same pitch range or below f0  6) Other contours higher than f0


sId = ones(length(stSamp),1);
eId = -1*ones(length(endSamp),1);
tId = [sId;eId];
timeStamp = [stSamp endSamp];
[TS1,indic1] = sort(timeStamp);
sorTID = tId(indic1);

TS2 = TS1.*sorTID';
[TS3,indic2,indic3]=unique(TS2,'first');
uniqid = sorTID(indic2);
TS3 = TS3.*uniqid';
[TS3,indic4] = sort(TS3);
uniqid = uniqid(indic4);

disp('Begin Processing blocks...')
preflag = 0;

for k = 1:length(TS3)-1
%     tic    
    if abs(TS3(k)-TS3(k+1))==1
        P1 = 0; P2 = 0;
    else
        P1 = (1-uniqid(k))/2;
        P2 = (-1-uniqid(k+1))/2;
    end
    
    if (TS3(k)==TS3(k+1) & uniqid(k)==-1)
        P1 = 0; P2 = 0;
    end
    
    st_chunk = TS3(k)+P1;
    et_chunk = TS3(k+1)+P2;
    
    if (preflag == 1)
        temp = [st_chunk+1:et_chunk];
        temp = temp+1;
    else
        temp = [st_chunk:et_chunk];
        temp = temp+1;
    end
%     if isempty(temp)
%         disp('hurray')
%     end
    %temp = [timeStamp(k):timeStamp(k+1)]; temp = temp+1;
    flag = ~sum(pbuff(temp)~=1); % check whether the region of the block is within the annotation region
    stSamp_wf = temp(1)*128; endSamp_wf = temp(end)*128;
    wave1_clip=wave1((stSamp_wf-1024):(endSamp_wf+1025));
    if flag
        id1 = idmat((1:5),temp(1:end)); 
        u = unique(id1); % Check the contour 
        u = u(u~=-1);    % numbers in this region 
        ctrId = u;
        for f = 1:length(ctrId)
            [r c] = find(id1==ctrId(f));
            rowid = unique(r);
            x_chunk = temp(1:end); y_chunk = actrmat(rowid,x_chunk); % at one time chunk pick the contours that are present and process them
            
            negInd = find(y_chunk~=-1); % Check if the chunk contains a contour or not in that row.
            
            if ~isempty(negInd)
                x_chunk = x_chunk(negInd); y_chunk = y_chunk(negInd);
                %%%%%%%%%%%%%%%%%%%% Calculate the match of the contour %%%%%%%%%%%%%%%%%%%%
                err1=sum(abs(y_chunk'-pctr(x_chunk)))/length(y_chunk);
                err2=sum(abs(y_chunk'-2*pctr(x_chunk)))/length(y_chunk);
                err3=sum(abs(y_chunk'-3*pctr(x_chunk)))/length(y_chunk);
                errhalf=sum(abs(y_chunk'-0.5*pctr(x_chunk)))/length(y_chunk);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if err1 < 0.01
                    voc(cnt1).xcont = x_chunk; voc(cnt1).ycont = y_chunk;
                    vocPitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(vocPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd); % get Harmonic mags and locs 
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                    [w z] = size(temp1);

                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        vocContour(cnt1).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        vocContour(cnt1).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    vocContour(cnt1).aggslopes = R2.AggregateSlope;
                    
                    if z == 1
                        aggrMag_Voc(cnt1,:) = linMag./sum(linMag);
                    else
                        aggrMag_Voc(cnt1,:) = temp1./sum(temp1);
                    end
                    vocconInd(cnt1) = ctrId(f);
                    cnt1 = cnt1 + 1;
                    %plot(x_chunk,y_chunk,'.r')
                elseif err2 < 4
                    harm2(cnt2).xcont = x_chunk; harm2(cnt2).ycont = y_chunk;
                    harm2Pitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(harm2Pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                    [w z]=size(temp1);
                    
                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        harm2Contour(cnt2).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        harm2Contour(cnt2).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    harm2Contour(cnt2).aggslopes = R2.AggregateSlope;
                    
                    if z == 1
                        aggrMag_harm2(cnt2,:) = linMag./sum(linMag);
                    else
                        aggrMag_harm2(cnt2,:) = temp1./sum(temp1);
                    end
                    harm2Ind(cnt2) = ctrId(f);
                    cnt2 = cnt2 + 1;
                    %plot(x_chunk,y_chunk,'.g')
                elseif err3 < 4
                    harm3(cnt6).xcont = x_chunk; harm3(cnt6).ycont = y_chunk;
                    harm3Pitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(harm3Pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                       [w z]=size(temp1);
                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        harm3Contour(cnt6).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        harm3Contour(cnt6).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    harm3Contour(cnt6).aggslopes = R2.AggregateSlope;
                     
                    if z == 1
                        aggrMag_harm3(cnt6,:) = linMag./sum(linMag);
                    else
                        aggrMag_harm3(cnt6,:) = temp1./sum(temp1);
                    end
                    harm3Ind(cnt6) = ctrId(f);
                    cnt6 = cnt6 + 1;
                    %plot(x_chunk,y_chunk,'*r')
                elseif errhalf < 4
                    subharm(cnt3).xcont = x_chunk; subharm(cnt3).ycont = y_chunk;
                    subharmPitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(subharmPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                    [w z] =size(temp1);
                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        subharmContour(cnt3).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        subharmContour(cnt3).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    subharmContour(cnt3).aggslopes = R2.AggregateSlope;

                    if z==1
                        aggrMag_subharm(cnt3,:) = linMag./sum(linMag);
                    else
                        aggrMag_subharm(cnt3,:) = temp1./sum(temp1);
                    end
                    
                    subharmInd(cnt3) = ctrId(f);
                    cnt3 = cnt3 + 1;
                    %plot(x_chunk,y_chunk,'.k')

                elseif (median(y_chunk)>median(pctr(x_chunk)))

                    viogr(cnt4).xcont = x_chunk; harm2(cnt4).ycont = y_chunk;
                    viogrPitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(viogrPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                    [w z] = size(temp1);

                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        viogrContour(cnt4).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        viogrContour(cnt4).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    viogrContour(cnt4).aggslopes = R2.AggregateSlope;
                    if z == 1
                        aggrMag_viogr(cnt4,:) = linMag./sum(linMag);
                    else
                        aggrMag_viogr(cnt4,:) = temp1./sum(temp1);
                    end
                    viogrInd(cnt4) = ctrId(f);
                    cnt4 = cnt4 + 1;
                    %plot(x_chunk,y_chunk,'.m')

                elseif (median(y_chunk)<median(pctr(x_chunk)))
                    
                    violr(cnt5).xcont = x_chunk; violr(cnt5).ycont = y_chunk;
                    violrPitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(violrPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                    [w z] = size(temp1);

                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        violrContour(cnt5).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        violrContour(cnt5).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    violrContour(cnt5).aggslopes = R2.AggregateSlope;

                    if z == 1
                        aggrMag_violr(cnt5,:) = linMag./sum(linMag);
                    else
                        aggrMag_violr(cnt5,:) = temp1./sum(temp1);
                    end

                    violrInd(cnt5) = ctrId(f);
                    cnt5 = cnt5 + 1;
                    %plot(x_chunk(negInd),y_chunk(negInd),'*k')

                else
                    
                    vioctr(cnt7).xcont = x_chunk; vioctr(cnt5).ycont = y_chunk;
                    vioctrPitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(vioctrPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                    temp1 = sum(linMag);
                    [w z] = size(temp1);
                    R2 = genSlope(hmags,slpnum); % Calculate Slope!!!!!!!!

                    if iscolumn(R2.SLOPE)
                        vioctrContour(cnt7).frameslopes = [x_chunk' R2.SLOPE];
                    else
                        vioctrContour(cnt7).frameslopes = [x_chunk' R2.SLOPE'];
                    end

                    vioctrContour(cnt7).aggslopes = R2.AggregateSlope;

                    if z == 1
                        aggrMag_vioctr(cnt7,:) = linMag./sum(linMag);
                    else
                        aggrMag_vioctr(cnt7,:) = temp1./sum(temp1);
                    end

                    vioctrInd(cnt7) = ctrId(f);
                    cnt7 = cnt7 + 1;

                end
            end
        end
    end
    if TS3(k)==TS3(k+1)
        preflag = 1;
    else
        preflag = 0;
    end
%     toc
end

%% Return Values computed

R.VocFeat = aggrMag_Voc;
%R.VocFeatTime=hfeat_locs_voc;
if exist('aggrMag_harm2')==1
    R.harm2Feat = aggrMag_harm2;
else
    R.harm2Feat = zeros(1,nHarm);
end
if exist('aggrMag_harm3')
    R.harm3Feat = aggrMag_harm3;
else
    R.harm3Feat = zeros(1,nHarm);
end
if exist('aggrMag_subharm')
    R.subharmFeat = aggrMag_subharm;
else
    R.subharmFeat = zeros(1,nHarm);
end
if exist('aggrMag_viogr')
    R.viogrFeat = aggrMag_viogr;
else
    R.viogrFeat = zeros(1,nHarm);
end
if exist('aggrMag_vioctr')
    R.vioctrFeat = aggrMag_vioctr;
else
    R.vioctrFeat = zeros(1,nHarm);
end
if exist('aggrMag_violr')
    R.violrFeat = aggrMag_violr;
else
    R.violrFeat = zeros(1,nHarm);
end
if exist('vocContour')
    R.vocFrameSlopes = vocContour;
else
    vocDummy(1).frameslopes = [0 inf];
    R.vocFrameSlopes = vocDummy;
end
if exist('harm2Contour')
    R.harm2FrameSlopes = harm2Contour;
else
    harm2Dummy(1).frameslopes = [0 inf];
    R.harm2FrameSlopes = harm2Dummy;
end
if exist('harm3Contour')
    R.harm3FrameSlopes = harm3Contour;
else
    harm3Dummy(1).frameslopes = [0 inf];
    R.harm3FrameSlopes = harm3Dummy;
end
if exist('subharmContour')
    R.subharmFrameSlopes = subharmContour;
else
    subharmDummy(1).frameslopes = [0 inf];
    R.subharmFrameSlopes = subharmDummy;
end
if exist('viogrContour')
    R.viogrFrameSlopes = viogrContour;
else
    viogrDummy(1).frameslopes = [0 inf];
    R.viogrFrameSlopes = viogrDummy;
end
if exist('R.vioctrContour')
    R.vioctrFrameSlopes = vioctrContour;
else
    vioctrDummy(1).frameslopes = [0 inf];
    R.vioctrFrameSlopes = vioctrDummy;
end
if exist('violrContour')
    R.violrFrameSlopes =violrContour;
else
    violrDummy(1).frameslopes = [0 inf];
    R.violrFrameSlopes = violrDummy;
end


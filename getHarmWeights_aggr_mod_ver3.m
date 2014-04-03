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


sId = ones(length(stSamp),1); % 1 if the frame under consideration is a start of the contour
eId = -1*ones(length(endSamp),1); % -1 if the frame under consideration is a end of the contour
tId = [sId;eId]; % append all the ids 
timeStamp = [stSamp endSamp]; % append start times and end times
[TS1,indic1] = sort(timeStamp); % sort the times to get the segments 
sorTID = tId(indic1); % sort indices according to sorted time stamps
%%%%%% Processing the timestamps to give unique timestamps %%%%
TS2 = TS1.*sorTID'; % storing the timestamps with the sign according to start or end if start it is 1 if end its -2
[TS3,indic2,indic3]=unique(TS2,'first'); %uniquing with sign
uniqid = sorTID(indic2);
TS3 = TS3.*uniqid'; % undoing the sign to the timestamps and storing the signs in the uniqid array separately
[TS3,indic4] = sort(TS3); 
uniqid = uniqid(indic4);

disp('Begin Processing blocks...')
preflag = 0;

for k = 1:length(TS3)-1

    % applying the conditions of start timestamp and endtimestamp to see
    % where to start the segment and where to end the segment for
    % processing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     if isempty(temp)
%         disp('hurray')
%     end
    %temp = [timeStamp(k):timeStamp(k+1)]; temp = temp+1;
    flag = ~sum(pbuff(temp)~=1); % check whether the region of the block is within the annotation region
    stSamp_wf = temp(1)*128; endSamp_wf = temp(end)*128; % index in the waveform
    wave1_clip=wave1((stSamp_wf-1024):min((endSamp_wf+1024),length(wave1))); % clip the waveform corresponding to the segment
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
                
                %%%%%%%%%%%%%%%%%%%% Calculate the error of the contour with each class %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                err1=sum(abs(y_chunk'-pctr(x_chunk)))/length(y_chunk); % error with fundamental ie the predominant final o/p
                err2=sum(abs(y_chunk'-2*pctr(x_chunk)))/length(y_chunk); % error with 2nd harmonic of predominant final o/p
                err3=sum(abs(y_chunk'-3*pctr(x_chunk)))/length(y_chunk); % error with 3rd harmonic of predominant final o/p
                errhalf=sum(abs(y_chunk'-0.5*pctr(x_chunk)))/length(y_chunk); % error with subharmonic of predominant final o/p
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if err1 < 0.01
                    
                    vocPitcharr = [x_chunk' y_chunk']; % time freq of the chunk
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(vocPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd); % get Harmonic mags and locs 
                    hmags = hmags';
                    linMag = 10.^(hmags./20); % convert to linear scale
%                     temp1 = sum(linMag);
%                     [w z] = size(temp1);

                    R2 = genSlope(hmags,hmagsLoc,slpnum); % Calculate Slope!!!!!!!!
                    vocframeslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_Voc(cnt1,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_Voc(cnt1,:) = temp1./sum(temp1);
%                     end
%                   
                    % To store all attributes in a structure
                    vocmags = linMag; vocLocs = hmagsLoc';
                    vpthout = [(x_chunk'*(128/44100)) y_chunk'];
                    vocalContour(cnt1).freq = vpthout;
                    vocalContour(cnt1).id = ctrId(f);
                    vocalContour(cnt1).harmMag = vocmags;
                    vocalContour(cnt1).harmLocs = vocLocs;
                    vocalContour(cnt1).frameslopes = vocframeslopes;
                    
                    cnt1 = cnt1 + 1; % increment count
                    
                elseif err2 < 4
                    
                    harm2Pitcharr = [x_chunk' y_chunk']; 
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(harm2Pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
%                     temp1 = sum(linMag);
%                     [w z]=size(temp1);
                    
                    R2 = genSlope(hmags,hmagsLoc,slpnum); % Calculate Slope!!!!!!!!
                    harm2frameslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_harm2(cnt2,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_harm2(cnt2,:) = temp1./sum(temp1);
%                     end
                    harm2mags = linMag; harm2Locs = hmagsLoc';
                    harm2pthout = [(x_chunk'*(128/44100)) y_chunk'];
                    harm2Contour(cnt2).freq = harm2pthout;
                    harm2Contour(cnt2).id = ctrId(f);
                    harm2Contour(cnt2).harmMag = harm2mags;
                    harm2Contour(cnt2).harmLocs = harm2Locs;
                    harm2Contour(cnt2).frameslopes = harm2frameslopes;
                    
                    cnt2 = cnt2 + 1;
                    
                    %plot(x_chunk,y_chunk,'.g')
                elseif err3 < 4
                    
                    harm3Pitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(harm3Pitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
                     
%                     temp1 = sum(linMag);
%                     [w z]=size(temp1);
                    R2 = genSlope(hmags, hmagsLoc, slpnum); % Calculate Slope!!!!!!!!
                    harm3frameslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_harm3(cnt3,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_harm3(cnt3,:) = temp1./sum(temp1);
%                     end
                    harm3mags = linMag; harm3Locs = hmagsLoc';
                    harm3pthout = [(x_chunk'*(128/44100)) y_chunk'];
                    harm3Contour(cnt3).freq = harm3pthout;
                    harm3Contour(cnt3).id = ctrId(f);
                    harm3Contour(cnt3).harmMag = harm3mags;
                    harm3Contour(cnt3).harmLocs = harm3Locs;
                    harm3Contour(cnt3).frameslopes = harm3frameslopes;
                    
                    cnt3 = cnt3 + 1;
                    
                    %plot(x_chunk,y_chunk,'*r')
                elseif errhalf < 4
                    
                    subharmPitcharr = [x_chunk' y_chunk']; 
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(subharmPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
%                     temp1 = sum(linMag);
%                     [w z]=size(temp1);
                    
                    R2 = genSlope(hmags, hmagsLoc, slpnum); % Calculate Slope!!!!!!!!
                    subharmframeslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_subharm(cnt4,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_subharm(cnt4,:) = temp1./sum(temp1);
%                     end
                    subharmmags = linMag; subharmLocs = hmagsLoc';
                    subharmpthout = [(x_chunk'*(128/44100)) y_chunk'];
                    subharmContour(cnt4).freq = subharmpthout;
                    subharmContour(cnt4).id = ctrId(f);
                    subharmContour(cnt4).harmMag = subharmmags;
                    subharmContour(cnt4).harmLocs = subharmLocs;
                    subharmContour(cnt4).frameslopes = subharmframeslopes;
                    
                    cnt4 = cnt4 + 1;
                    
                    %plot(x_chunk,y_chunk,'.k')

                elseif (median(y_chunk)>median(pctr(x_chunk)))
                    
                    viogrPitcharr = [x_chunk' y_chunk']; 
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(viogrPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
%                     temp1 = sum(linMag);
%                     [w z]=size(temp1);
                    
                    R2 = genSlope(hmags,hmagsLoc,slpnum); % Calculate Slope!!!!!!!!
                    viogrframeslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_viogr(cnt5,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_viogr(cnt5,:) = temp1./sum(temp1);
%                     end

                    viogrmags = linMag; viogrLocs = hmagsLoc';
                    viogrpthout = [(x_chunk'*(128/44100)) y_chunk'];
                    viogrContour(cnt5).freq = viogrpthout;
                    viogrContour(cnt5).id = ctrId(f);
                    viogrContour(cnt5).harmMag = viogrmags;
                    viogrContour(cnt5).harmLocs = viogrLocs;
                    viogrContour(cnt5).frameslopes = viogrframeslopes;
                    
                    cnt5 = cnt5 + 1;



                elseif (median(y_chunk)<median(pctr(x_chunk)))
                    
                    violrPitcharr = [x_chunk' y_chunk']; 
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(violrPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
%                     temp1 = sum(linMag);
%                     [w z]=size(temp1);
                    
                    R2 = genSlope(hmags,hmagsLoc, slpnum); % Calculate Slope!!!!!!!!
                    violrframeslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_violr(cnt6,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_violr(cnt6,:) = temp1./sum(temp1);
%                     end
                    
                    violrmags = linMag; violrLocs = hmagsLoc';
                    violrpthout = [(x_chunk'*(128/44100)) y_chunk'];
                    violrContour(cnt6).freq = violrpthout;
                    violrContour(cnt6).id = ctrId(f);
                    violrContour(cnt6).harmMag = violrmags;
                    violrContour(cnt6).harmLocs = violrLocs;
                    violrContour(cnt6).frameslopes = violrframeslopes;
                    
                    cnt6 = cnt6 + 1;

                else
                    
                    vioctrPitcharr = [x_chunk' y_chunk'];
                    [hmags,hmagsLoc]=HarmonicSubtraction_ver3(vioctrPitcharr, wave1_clip, win, NFFT, nHarm, thsld, maxhd);
                    hmags = hmags';
                    linMag = 10.^(hmags./20);
%                     temp1 = sum(linMag);
%                     [w z]=size(temp1);
                    
                    R2 = genSlope(hmags,hmagsLoc,slpnum); % Calculate Slope!!!!!!!!
                    vioctrframeslopes = R2.SLOPE;
                    
                    % Check if theres only one frame or more. If one then
                    % dont sum if more then one then take aggregate.
                    
%                     if z == 1
%                         aggrMag_vioctr(cnt7,:) = linMag./sum(linMag);
%                     else
%                         aggrMag_vioctr(cnt7,:) = temp1./sum(temp1);
%                     end
                    
                    vioctrmags = linMag; vioctrLocs = hmagsLoc';
                    vioctrpthout = [(x_chunk'*(128/44100)) y_chunk'];
                    vioctrContour(cnt7).freq = vioctrpthout;
                    vioctrContour(cnt7).id = ctrId(f);
                    vioctrContour(cnt7).harmMag = vioctrmags;
                    vioctrContour(cnt7).harmLocs = vioctrLocs;
                    vioctrContour(cnt7).frameslopes = vioctrframeslopes;
                    
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

if exist('vocalContour')==1
    R.VocFeat = vocalContour;
else
    vocaldummy(1).freq = inf;
    vocaldummy(1).id = inf;
    vocaldummy(1).harmMag = inf*ones(1,nHarm);
    vocaldummy(1).harmLocs = inf*ones(1,nHarm);
    vocaldummy(1).frameslopes = inf;
    R.VocFeat = vocaldummy;
end

if exist('harm2Contour')==1
    R.harm2Feat = harm2Contour;
else
    harm2dummy(1).freq = inf;
    harm2dummy(1).id = inf;
    harm2dummy(1).harmMag = inf*ones(1,nHarm);
    harm2dummy(1).harmLocs = inf*ones(1,nHarm);
    harm2dummy(1).frameslopes = inf;
    R.harm2Feat = harm2dummy;
end

if exist('harm3Contour')==1
    R.harm3Feat = harm3Contour;
else
    harm3dummy(1).freq = inf;
    harm3dummy(1).id = inf;
    harm3dummy(1).harmMag = inf*ones(1,nHarm);
    harm3dummy(1).harmLocs = inf*ones(1,nHarm);
    harm3dummy(1).frameslopes = inf;
    R.harm3Feat = harm3dummy;
end

if exist('subharmContour')==1
    R.subharmFeat = subharmContour;
else
    subharmdummy(1).freq = inf;
    subharmdummy(1).id = inf;
    subharmdummy(1).harmMag = inf*ones(1,nHarm);
    subharmdummy(1).harmLocs = inf*ones(1,nHarm);
    subharmdummy(1).frameslopes = inf;
    R.subharmFeat = subharmdummy;
end

if exist('viogrContour')==1
    R.viogrFeat = viogrContour;
else
    viogrdummy(1).freq = inf;
    viogrdummy(1).id = inf;
    viogrdummy(1).harmMag = inf*ones(1,nHarm);
    viogrdummy(1).harmLocs = inf*ones(1,nHarm);
    viogrdummy(1).frameslopes = inf;
    R.viogrFeat = harm3dummy;
end

if exist('violrContour')==1
    R.violrFeat = violrContour;
else
    violrdummy(1).freq = inf;
    violrdummy(1).id = inf;
    violrdummy(1).harmMag = inf*ones(1,nHarm);
    violrdummy(1).harmLocs = inf*ones(1,nHarm);
    violrdummy(1).frameslopes = inf;
    R.violrFeat = violrdummy;
end

if exist('vioctrContour')==1
    R.vioctrFeat = vioctrContour;
else
    vioctrdummy(1).freq = inf;
    vioctrdummy(1).id = inf;
    vioctrdummy(1).harmMag = inf*ones(1,nHarm);
    vioctrdummy(1).harmLocs = inf*ones(1,nHarm);
    vioctrdummy(1).frameslopes = inf;
    R.vioctrFeat = vioctrdummy;
end



a=load('1_02_Rama_Dayajudave.wav.allctrs');
%fp=load('Kapi_alapana_KVN_clip.pitch_2');
%fp=fp(:,2);
pc=load('1_02_Rama_Dayajudave.wav.pitch');
pc=pc(:,2);
plot(pc,'.')
Fs=44100;
hop=128;
[m n]=size(a);
for i=1:length(ind)
    a11=a(ind(i),:);
    a11=a11(a11~=inf);
    contour=a11(3:end);
    st=a11(1,2);
    stSamp=round((st*Fs)/hop);
    lenCtr=length(contour);
    endSamp=stSamp+lenCtr;
    xCtr=[stSamp:endSamp-1];
    yCtr=55*(2.^(((10.*contour)-1)./1200));
    
%     if length(xCtr)<length(yCtr)
%         yCtr=yCtr(1:length(xCtr));
%     elseif length(yCtr)<length(xCtr)
%         xCtr=xCtr(1:length(yCtr));
%     end
    hold on
    plot(xCtr,yCtr,'.r')
     pause
    clear xCtr
    clear yCtr
end

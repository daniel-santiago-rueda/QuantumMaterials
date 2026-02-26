function arrfourier=fourierdatav2(arrsubs,minH)

arrsubs(:,1)=1./arrsubs(:,1);

[sizey,sizex]=size(arrsubs);


array=[];
for i=1:1:sizey
    if arrsubs(i,1)<=1/minH 
        array=[array; arrsubs(i,1) arrsubs(i,2)];
    end
end

%arrsubs=subtractbackground(arr(:,1),arr(:,2), 3, Xmin,Xmax)
minfieldtot=min(array(:,1));
maxfieldtot=max(array(:,1));

xi=transpose(minfieldtot:(maxfieldtot-minfieldtot)/(sizey*10):maxfieldtot);
[x1p, idxp]=unique(array(:,1));
yi1=interp1(array(idxp,1),array(idxp,2),xi);

array=[];
array(:,1)=xi;
array(:,2)=yi1;

%array=inverseBarray(arrsubs,points);
%array=arrsubs;
[sizey,sizex]=size(array);
fs=1/(array(2,1)-array(1,1));

win=hann(sizey); %Hanning window. If other window, use hamming(n), turkeywin(n), kaiser(n)
%win=kaiser(sizey); %Hanning window. If other window, use hamming(n), turkeywin(n), kaiser(n)
xw = win(:).*array(:,2); %Applies Hanning window to data
%xw = array(:,2); %If no window applied

NFFT = 2^nextpow2(sizey*10);
%NFFT = 2^nextpow2(sizey*100);
fftarray=fft(xw,NFFT);
%f = array(sizey,1)*linspace(0,1,NFFT/2+1);
f = (fs/2)*linspace(0,1,NFFT/2+1);

arrfourier(:,1)=f;
arrfourier(:,2)=2*abs(fftarray(1:NFFT/2+1));

%x2=array(:,1);
%y2=array(:,2);
%subplot(1,2,1); plot(x2,y2,'-b','LineWidth',2)
%subplot(1,2,2); plot(f,2*abs(fftarray(1:NFFT/2+1)), '-r','LineWidth',2);
%xlim([0,500]);
%xlabel('frequency (T)')
%ylabel('FFT')
%title('Fourier Transform');


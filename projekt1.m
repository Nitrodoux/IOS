%Identyfikacja obiekt�w projekt
%wczytanie sygna��w
load('T1425.MAT')%sygna� wej�ciowy
load('T1426.MAT')%sygna� eksperymentalny
load('T1427.MAT')%sygna� wyj�ciowy
load('arx80901')
%% 1. Przebiegi czasowe
figure(1)
plot(c1x,c1)
ylabel('Napi�cie [V]');
grid on;xlabel('Czas[s]');title('Sygna� wymuszenia u(t)');
figure(2)
plot(c2x,c2)
grid on;title('Sygna� pomiarowy y(t)');
xlabel('Czas[s]');ylabel('Napi�cie [V]');

%% 2. Transformacja sygna��w z t do f
U=fft(c1);%fft sygna�u wej�ciowego
Y=fft(c2);%fft sygna�u wyj�ciowego
G=Y./U;

%% 3. Eksperymentalny wykres Bodego 
L=20*log(abs(G)); %transformacja sygna�u z dziedziny t do f
Q=atan(imag(G)./real(G));

fs=2048; %cz�stotliwo�� pr�bkowania
n=length(c1); %d�. sygna�u
f=0:1:n-1; 
kx=(f/n)*fs;
fi = atan(imag(G)./real(G));

figure(3)
subplot(2,1,1);
semilogx(kx,L)
xlabel('Czestotliwo�� [Hz]');
ylabel('Amplituda [dB]');
title('Eksperymentalny wykres Bodego- charakterystyka amplitudowa')
subplot(2,1,2);
semilogx(kx,fi)
xlabel('Czestotliwo�� [Hz]');
ylabel('Amplituda [rad]');
title('Eksperymentalny wykres Bodego- charakterystyka fazowa')

figure(4)
semilogx(f, 20*log(abs(G)));
grid on;
title('Wykres Bodego uzyskany z transmitancji G');
xlabel('Cz�stotliwo�� [Hz]');
ylabel('Amplituda [dB]');
axis([10 800 -200 50]);

% transmitancja sygna�u eksperymentalnego
G2=o2i1./o2i1x;
 
L2=20*log10(abs(G2)); % Trasnformacja sygna��w z dziedziny czasu do dziedziny cz�stotliwo�ci
Q2=atan(imag(G2)./real(G2));
 
fs=201;
n=length(o2i1x);
f=0:1:n-1;
kx2=(f/n)*fs;
% 4 Estymacja parametr�w w IT
%% Por�wnanie dopasowania modelu estymowanego modelu obiektu
[licz,mian]=tfdata(arx80901,'v');
[A,B,C,D]=tf2ss(licz,mian);
sys=ss(A,B,C,D,1/2048);
[mag,phase,wout] = bode(sys);
kat=wout/6.28;
figure(5)
mag=squeeze(mag);
semilogx(kat,20*log10(mag))
grid on;
title('Wykres Bodego uzyskany z modelu arx80901 ');
xlabel('Cz�stotliwo�� [Hz]');ylabel('Amplituda [dB]');
axis([10 800 -100 40]);

figure(6)
mag=squeeze(mag);
semilogx(kat,20*log10(mag))
hold on
semilogx(o2i1x,20*log10((abs(o2i1))));
title('Charakterystyka Bodego');
grid on;
title('Charakterystyka Bodego - por�wnanie modeli ');
xlabel('Cz�stotliwo�� [Hz]');ylabel('Amplituda [dB]');
legend('estymowany model obiektu','model obiektu','location','best');
axis([10 800 -100 40]);
%% 6. Redukcja rz�du modelu
[sysb,g]=balreal(sys)
sys_red_new=modred(sysb,g<1/2048)
%% 7. Porownanie zredukowanego rz�du modelu obiektu z danymi eksperymentalnymi

[mag2,phase2,wout2]=bode(sys_red_new)
kat2=wout2/6.28;
mag2=squeeze(mag2);
figure(7)
semilogx(kat2,20*log10(mag2));
hold on
semilogx(kat,20*log10(mag));
hold on
semilogx(o2i1x,20*log10((abs(o2i1))),'-k');
title('Por�wnanie charakterystyk amplitudowych','FontSize',12,'FontWeight','Bold')
xlabel('Cz�stotliwo�� [Hz]')
ylabel('Amplituda [dB]')
legend('Model niezredukowany 90 rz�du','model zredukowany 47 rz�du','model obiektu','location','best')
axis([10 800 -100 40]);


[licz1,mian1]=tfdata(sys_red_new,'v');
[Ar,Br,Cr,Dr]=tf2ss(licz1,mian1);
Gr=ss(Ar,Br,Cr,Dr,0.0001);
[z,p,k]=ss2zp(Ar,Br,Cr,Dr)

%% 8. Dalsza redukcja rz�du modelu obiektu o wydzielenie pasma cz�stotliwo�ci zawieraj�cych
%pierwsze 3 cz�stotliwo�ci rezonansowe i anytreznansowe.
%(bieguny i zera modelu obiektu poza tym zakresem nale�y odrzuci�)

dty=0.00048;
 
for i=1:length(z)
    omegaA(i,1)=abs(log(z(i,1)))/dty;
end
 
for i=1:length(p)
    omegaB(i,1)=abs(log(p(i,1)))/dty;
end
 
bieg=[p,omegaB/6.28];
zera=[z,omegaA/6.28];
 
posortbieg = sortrows(bieg,2);
posortzera = sortrows(zera,2);
Gr2=zpk(oA1,oB1,0.01,dty);
[mag3,phase3,wout3]=bode(Gr2);
kat3=wout3/6.28;

%% 9. Por�wnania modelu zredukowanego_2 z danymi eksperymentlanymi. (Bode)
figure(8)
semilogx(kat3,20*log10(mag3(:,:)));
hold on
semilogx(kat2,20*log10(mag2));
hold on
semilogx(kat,20*log10(mag));
hold on
semilogx(o2i1x,20*log10((abs(o2i1))),'-k');
title('Por�wnanie charakterystyk amplitudowych','FontSize',12,'FontWeight','Bold')
xlabel('Cz�stotliwo�� [Hz]','FontSize',12)
ylabel('Amplituda [dB]','FontSize',12)
legend('model zredukowany 23 rz�du','model niezredukowany 87 rz�du','model zredukowany 45 rz�du','model obiektu','location','best')
axis([10 800 -100 40]);


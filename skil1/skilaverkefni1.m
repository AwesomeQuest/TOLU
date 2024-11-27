clc;clear all;
% part 1
lambda = 0.8; % g/mm^3 eða 
EI = 1.09*1e+10; % Pa*mm^2
L = 30; % mm
%k = ((w^2.*lambda)/(EI))^(1/4);
%cos(kL)*cosh(kL) = -1;
x = linspace(0, 10, 100);
plot(x, f(x))
yline(0)
ylim([-30,30])
xlim([0, 10])

%% part 2
a=1;b=3;tol=0.00001;
kL1 = bisect(a, b, tol);
k1 = kL1/L;
w1 = k1^2*sqrt(EI/lambda)

%% part 3
deltan = 10^(-5); % óvissan
n = log2((b-a)/(2*deltan))

%% part 4
x0=5;tol=0.0001;
kL2 = newton(x0,tol);
k2 = kL2/L;
w2 = k2^2*sqrt(EI/lambda)

%% part 5
tol=1e-15;
xn=2;
newton_roots=[];
for i=1:20
    root=newton(xn,tol);
    xn=xn+3.1;
    newton_roots=[newton_roots root];

end
%disp(newton_roots)
insert=f(newton_roots);
double(insert)
fprintf(" %0.5f", newton_roots)

%% part 6
L=30;lambda=0.8;EI=1.09*10^(10);
k=0.1565;
w=2.8578e+3;
delta=-3.5;

T = 2*pi/w; % Ein lota

% Búa til glugga til að teikna hreyfimyndinna í
% Við felum hann á meðan við teiknum.
hfig = figure('visible', 'off');

FPS = 30; % Fjöldi ramma á sek í hreyfimyndinni.
          % Til að hægja á hreyfimyndinni er hægt að auka FPS
          % og halda FPS_PLAY föstu.

FPS_PLAY = 30; % Fjöldi ramma á sek þegar við spilum hreyfimyndinna.
               % Best að hafa það 30.

t_start = 0.0; % Byrjunar tími í sek
t_end = 2*T; % Enda tími í sek (Tvær lotur)

% Reikna út fjölda ramma sem þarf og stærð eins tíma skrefs.
n = ceil(FPS*(t_end - t_start));
dt = (t_end - t_start)/n;

F(1:n) = getframe(hfig);
t = zeros(1, n);

% Teikna hvert tíma skref
for i = 1:n
    % Hreinsa allt úr glugganum
    cla
    hold on
    
    % Stilla stærð ása.
    % Best er að þeir séu alltaf jafn stórir í hverjum ramma.
    axis([0.0, L, -5.0, 5.0]) %[x_min, x_max, y_min, y_max]
    %axis equal
    
    % Reikna tímann
    t(i) = t_start + i*dt;
    
    % Teikna dót - þessu þarf að breyta og finna (x,y) hnit
    x=linspace(0,L,100);
    plot(x, A(x,k,L)*B(t(i),k,delta,L,w))
    
    % Ef við viljum sýna tímann á grafinu
    text(5.0, 5.0, strcat(num2str(t(i)), ' s'))
    
    grid on % Grid on eða off
    
    % Merkja ásanna og setja titil
    xlabel('Staðsetning á bitanum (mm)')
    ylabel('Hliðrun á bitanum (mm)')
    title('Sveiflan á bitanum')
    
    % Geyma rammann
    F(i) = getframe(hfig);
end

% Sýna hreyfimyndina einu sinni
close all
movie(gcf, F, 1, FPS_PLAY)

% Skrifa hreyfimyndinni niður a diskinn
v = VideoWriter('Video_1.mp4', 'MPEG-4');
v.FrameRate = FPS_PLAY;
open(v);
writeVideo(v, F);
close(v);
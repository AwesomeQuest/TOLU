length=30;
lambda=0.8; %g/m
EI=1.09e10 ;
%% part 1

x=linspace(0,100,100)
plot(x,f(x))
%xlabel(x)
%ylabel(f(x))
xlim([0,10])
ylim([-20,20])
yline(0)

%%
% PART 2 helmingunaraðferð
tol=1e-5 %
% fyrsta rót
disp("--fyrsta rót")
bisect(1,2,tol)
% Önnur rót
disp("--önnur rót")
bisect(4,6,tol)
%þriðja rót
disp("--þriðja rót")
bisect(7,8,tol)
%%
%PART 3
%counter var bætt inn og prentar, sjá part 2

%%
%PART 4
%fyrsta rót
fprintf("fyrsta rót %0.5f",newton(2,tol))
x1=newton(2,tol);
fprintf("1 tíðni %0.5f", freq(x1/30,EI,lambda))

%önnur rót
x2=newton(5,tol);
fprintf("önnur rót %0.5f",x2)
fprintf("tíðni 2 = %0.5f",freq(x2/length,EI,lambda))

%þriðja rót
fprintf("þriðja rót %0.5f", newton(9,tol))



%%

%Part 5 Newtonian
tol=1e-15;
xn=2; 
newton_roots=[];
for i=1:20
    root=newton(xn,tol);
    xn=xn+pi;
    newton_roots=[newton_roots root];

end
%disp(newton_roots)
insert=f(newton_roots);
double(insert);
fprintf(" %0.5f", newton_roots)
frequencies=freq(newton_roots/length,EI,lambda);
fprintf("%1.0f   ", frequencies)

%%
%plotting 
close all
x=linspace(0,100,100);
plot(x,f(x))
xlabel("x")
ylabel("f(x)")
xlim([0,62])
%ylim([-20])
%yline(0)
hold on
scatter((newton_roots),0,'k', 'filled')
%%
%% part 6
L=30;lambda=0.8;EI=1.09*10^(10);
k=0.1565;
w=2.8578e+3;
delta=-3.5;

T = 2*pi/w; % Ein lota

% Búa til glugga til að teikna hreyfimyndinna í
% Við felum hann á meðan við teiknum.
hfig = figure('visible', 'off');

FPS = 30000; % Fjöldi ramma á sek í hreyfimyndinni.
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
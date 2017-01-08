%% TP Modèles de Markov cachés

clc
clear all
close all



K = 2;
load('etat_cache_1.mat','x','tc');
n = size(x,2);
fprintf('instant de changement vrai %i\n',tc);
load('observation_1.mat','y');

for i=1:K
    axe = subplot(K+1,1,i);
    plot((x==i),'+r');
    axis([1 n -0.5 1.5]);
    set(axe,'YTick',[0 1]);
    title(['etat cache # ',num2str(i)]);
end

subplot(K+1,1,K+1);
plot(y);
title('observation')

%% Modèle Nominal

vo = [.2 .8]
Pio = [ .99 .01; .01 .99];
ho = [-2 +2];
Ro = [5 5];

%%Modele  alternatif
Pi=Pio;
h = [-10 10];
R = Ro;

%% Log vraisemblence

[logvrai pk vk] = markov(y ,vo,Pio,ho, Ro)














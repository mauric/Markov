%% TP Mod�les de Markov caches

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
    plot((x==i),'+r','LineWidth',1);
    axis([1 n -0.5 1.5]);
    set(axe,'YTick',[0 1]);
    title(['etat cache # ',num2str(i)]);
end

subplot(K+1,1,K+1);
plot(y,'LineWidth',1);
title('observation')

%Nominal model

vo = [.2 .8]
Pio = [ .99 .01; .01 .99];
ho = [-2 +2];
Ro = [5 5];

%%Modele  alternatif
Pi=Pio;
h = [-10 10];
R = Ro;

%% Log vraisemblenc
% alloc
n = size(y,2);
pk = zeros(n,2);
vk = zeros(2,n);

%Initialisation de pk
gio1 = (1/sqrt(2*pi*R(1)))*exp(-1/2*(y(1)-ho(1))^2/(R(1)));
gio2 = (1/sqrt(2*pi*R(2)))*exp(-1/2*(y(1)-ho(2))^2/(R(2)));

pk(1,1) = vo(1)*gio1/(vo(1)*gio1+vo(2)*gio2);
pk(1,2) = vo(2)*gio2/(vo(1)*gio1+vo(2)*gio2);

%Initialisation de vk
vk(:,n) = [1 ; 1];


for k = 2:n

    %forward
    gi1 = (1/sqrt(2*pi*R(1)))*exp(-1/2*(y(k)-ho(1))^2/(R(1)));
    gi2 = (1/sqrt(2*pi*R(2)))*exp(-1/2*(y(k)-ho(2))^2/(R(2)));
    g =[gi1 ; gi2];

    pk(k,:) =(pk(k-1,:)*Pio).*g'; %calcul de pk
    ck(k) = pk(k,:)*ones(2,1); %
    pk(k,:)= pk(k,:)/ck(k); %normalisation

end

 for k = n:-1:2

   gi1 = (1/sqrt(2*pi*R(1)))*exp(-1/2*(y(k)-ho(1))^2/(R(1)));
   gi2 = (1/sqrt(2*pi*R(2)))*exp(-1/2*(y(k)-ho(2))^2/(R(2)));
   g =[gi1 ; gi2];

  %backward
  vk(:,k-1) = Pi*(g.*vk(:,k));
  vk(:,k-1)= vk(:,k-1)/ck(k); %normalisation avec les memes ck
 end

 
lCk = log(ck);
lvk = log(vk);
lpk = log(pk);
lvrai = cumsum(lCk,1);
 

%% DOCUMENTATION
% sauvegarde les images pour le rapport
h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['back' num2str(length(h)+1-i)], 'png');
end

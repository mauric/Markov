% test : equations de Baum forward/backward + algorithme de Viterbi
clear; close all;
T = 1000;
% T = 10;

rep = input('modele 1/2/3 : ','s');
if isempty(rep)
    rep = '1';
end
switch rep
    case '1'
        disp('trois etats');
K = 3;
A = [0.9 0.05 0.05; 0.1 0.8 0.1; 0.05 0.15 0.8];
m = [-5 1 10];
sigma2 = [1 5 10];
p = [0.9 0.1 0];
    case '2'
        disp('deux etats, facile');
K = 2;
A = [0.99 0.01; 0.02 0.98];
m = [-2 2];
sigma2 = [5 5];
p = [0.2 0.8];
    case '3'
        disp('deux etats, difficile');
K = 2;
A = [0.95 0.05; 0.05 0.95];
m = [-1 1];
sigma2 = [3 3];
p = [0.5 0.5];
end

[X,Y] = gen(A,p,m,sigma2,T);
figure(1);
for i=1:K
    axe = subplot(K+1,1,i);
    plot((Y==i),'+r');
    axis([1 T -0.5 1.5]);
    set(axe,'YTick',[0 1])
    title(['etat cache # ',num2str(i)]);
end
subplot(K+1,1,K+1), plot(X);
title('observation');
disp('pause'); pause;

[alpha,beta,dens] = ForwardBackward(X,A,p,m,sigma2,K,T);
figure(2);
for i=1:K
    axe = subplot(K,1,i);
    plot((Y==i),'+r');
    axis([1 T 0.5 1.5]);
    hold on;
    plot(alpha(:,i),'b');
    axis([1 T -0.5 1.5]);
    set(axe,'YTick',[0 1])
    title(['etat cache # ',num2str(i),' et filtre (equation forward)']);
    hold off;
end
disp('pause'); pause;

figure(3);
[u,map] = max(alpha');
miss = find(Y~=map);
lmiss = length(miss);
disp('instants de mis-classification');
fprintf('instants de mis-clasification : %i\n',length(miss));
for i=1:K
    axe = subplot(K,1,i);
    plot((Y==i),'+r');
    axis([1 T 0.5 1.5]);
    hold on;
    plot((map==i),'g');
    axis([1 T -0.5 1.5]);
    set(axe,'YTick',[0 1])
    title(['etat cache # ',num2str(i),' et MAP marginal (equation forward)']);
    plot(miss,-0.4*ones(1,lmiss),'^k','MarkerFaceColor','k');
    hold off;
end
disp('pause'); pause;

figure(4);
for i=1:K
    axe = subplot(K,1,i);
    plot((Y==i),'+r');
    axis([1 T 0.5 1.5]);
    hold on;
    plot(alpha(:,i).*(beta(i,:))','b');
    axis([1 T -0.5 1.5]);
    set(axe,'YTick',[0 1])
    title(['etat cache # ',num2str(i),' et lisseur (equations forward / backward)']);
    hold off;
end
disp('pause'); pause;

figure(5);
[u,map] = max((alpha.*(beta'))');
miss = find(Y~=map);
lmiss = length(miss);
disp('instants de mis-classification');
fprintf('instants de mis-clasification : %i\n',length(miss));
for i=1:K
    axe = subplot(K,1,i);
    plot((Y==i),'+r');
    axis([1 T 0.5 1.5]);
    hold on;
    plot((map==i),'g');
    axis([1 T -0.5 1.5]);
    set(axe,'YTick',[0 1])
    title(['etat cache # ',num2str(i),' et MAP marginal (equations forward / backward)']);
    plot(miss,-0.4*ones(1,lmiss),'^k','MarkerFaceColor','k');
    hold off;
end
disp('pause'); pause;

[valeur,ante,dens] = Viterbi(X,A,p,m,sigma2,K,T);
path = zeros(K,T);
AFFICH = 0;
if(AFFICH)
    figure(7);
end
path(:,1) = [1:K]';
for t=2:T
    path(:,1:t-1) = path(ante(t,:),1:t-1);
    path(:,t) = [1:K]';
    if(AFFICH)
        for i=1:K
            plot(path(i,1:t),'-r');
            hold on;
        end
        axis([1 T 0.5 K+0.5]);
        hold off;
        drawnow;
    end
end

[c,I] = max(valeur(T,:));
map = path(I,:);
miss = find(Y~=map);
lmiss = length(miss);
disp('instants de mis-classification');
fprintf('instants de mis-clasification : %i\n',length(miss));
figure(6);
for i=1:K
    axe = subplot(K,1,i);
    plot((Y==i),'+r');
    axis([1 T -0.5 1.5]);
    hold on;
    plot((map==i),'g');
    axis([1 T -0.5 1.5]);
    set(axe,'YTick',[0 1]);
    title(['etat cache # ',num2str(i),' et MAP (algorithme de Viterbi)']);
    plot(miss,-0.4*ones(1,lmiss),'^k','MarkerFaceColor','k');
    hold off;
end

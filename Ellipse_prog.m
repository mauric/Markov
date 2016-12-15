%% TP FILTRE DE KALMAN
clc
close all
clear all


%% Load data
load('etat_cache.mat','x');
n = size(x,2);

figure(1);
plot(x(1,:),x(2,:),'-b');
title('etat cache');
axis_store=[xlim ylim];
axis square;

%variables et data pour le 
c = .9;
alpha = -pi/8;
epsilon = .1;
H = [1 0]
load('observation_1.mat','y_1');
R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
Fk = c*R;
y=y_1;
%% Initialisation
X_ch0_pred = [0 ; 0];
Qo_v = epsilon^2;
Qo_x = eye(2);
Po_pred = Qo_x; % poly pag 63
Ho = H;
X_ch0 = X_ch0_pred+Po_pred*Ho'*(inv(Ho*Po_pred*Ho'+Qo_v))*(y(:,1)-(Ho*X_ch0_pred));
Po = Po_pred-Po_pred*Ho'*(inv(Ho*Po_pred*Ho'+Qo_v))*Ho*Po_pred;
Qk_w = (1-c^2)*eye(2);
%% prediction

 X_ch_pred = X_ch0; %initialisation
 Pk_pred = Po_pred;
 
 % correction variables  
 X_ch = zeros(2,size(x,2));
 Pk = zeros(2,2,size(x,2));


X_ch(:,1) = X_ch0;
Pk(:,:,1) = Po_pred;

 
for i=2:size(x,2)
    %calcul de modele pour X_ch_pred
    X_ch_pred = Fk*X_ch(:,i-1);
    Pk_pred = Fk*Pk(:,:,i-1)*Fk'+Qk_w;
    %correction
    Ik = y(:,i)-(H*X_ch_pred);
    X_ch(:,i) =X_ch_pred+Pk_pred*H'*(inv(H*Pk_pred*H'+Qo_v))*Ik
    Pk(:,:,i) = Pk_pred-Pk_pred*H'*(inv(H*Pk_pred*H'+Qo_v))*H*Pk_pred;
    
end

figure(3);
plot(X_ch(1,:),'-r');
hold on
plot(x(1,:),'-b');
title('etat kalman');
axis square;

figure(4);
plot(X_ch(2,:),'-r');
hold on
plot(x(2,:),'-b');
title('etat kalman');
axis square;
%%
figure(5)
for k=1:n
    plot(x(1,1:k),x(2,1:k),'-b');
    hold on;
    ellipse_conf(X_ch(:,k), Pk(:,:,k));
    plot(x(1,k),x(2,k),'.b');
    title('title');
    axis(axis_store);
    axis square
    drawnow;
    pause(0.1);
    hold off;
end




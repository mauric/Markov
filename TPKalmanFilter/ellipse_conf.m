function ellipse_conf(fmean,fcov)
%-------------------------------------------------------------------------
% calcul et représentation graphique d'une ellipse de confiance à 90% 
% pour un v.a. gaussien bi-dimensionnel
%    de moyenne m et de matrice de covariance Q
% parametres :
%    fmean : moyenne m
%    fcov : matrice de covariance Q
%-------------------------------------------------------------------------
a = 4.6052; 
% si la boite-a-outil 'statistical toolbox' est disponible
% a = chi2inv(0.9,2);
npoints = 100; 
t = linspace(0,2*pi,npoints);
U = chol(fcov,'lower');
ellipse = sqrt(a)*U*[cos(t);sin(t)];
ellipse = ellipse+repmat([fmean(1);fmean(2)],1,npoints);
plot(ellipse(1,:),ellipse(2,:),'-k');
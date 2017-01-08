function [logvrai pk vk] = markov(y ,vo,Pio,ho, Ro)

%variables
Pi = Pio
R = Ro

% alloc 
n = size(y,2);
pk = zeros(n,2);
vk = zeros(2,n);

%Initialisation de pk
gio1 = (1/sqrt(2*pi*R(1)))*exp(-1/2*abs(y(1)-ho(1))^2*inv(R(1))*(y(1)-ho(1)));
gio2 = (1/sqrt(2*pi*R(2)))*exp(-1/2*abs(y(1)-ho(2))^2*inv(R(2))*(y(1)-ho(2)));

pk(1,1) = vo(1)*gio1/(vo(1)*gio1+vo(2)*gio2);
pk(1,2) = vo(2)*gio2/(vo(1)*gio1+vo(2)*gio2);

%Initialisation de vk
vk(:,n) = [1 ; 1];


for k = 2:n
    
    %forward
    gi1 = (1/sqrt(2*pi*R(1)))*exp(-1/2*conj(y(k)-ho(1))*inv(R(1))*(y(k)-ho(1)));
    gi2 = (1/sqrt(2*pi*R(2)))*exp(-1/2*conj(y(k)-ho(2))*inv(R(2))*(y(k)-ho(2)));
    g =[gi1 ; gi2];
      
    pk(k,:) =(pk(k-1,:)*Pio)*g; %calcul de pk
    ck(k) = pk(k,:)*ones(2,1); %
    pk(k,:)= pk(k,:)/ck(k); %normalisation
end

 for k = n:2    
  %backward
  vk(:,k-1) = Pi*(g.*vk(:,k));
  vk(:,k)= vk(:,k)/ck(k); %normalisation avec les memes ck
 end
 
    
 logvrai = cumsum(log(ck));



end


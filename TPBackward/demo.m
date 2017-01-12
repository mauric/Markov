% equations de Baum forward-backward

function [alpha,beta,dens,ll] = ForwardBackward(X,A,p,m,sigma2,K,T)

% densite d'emission

dens = ones(T,K);
dens = exp(-0.5*(X'*ones(1,K)-ones(T,1)*m).^2./(ones(T,1)*sigma2))./...
    sqrt(ones(T,1)*sigma2);

% filtre avant

alpha = ones(T,K);
alpha(1,:) = p.*dens(1,:);
c(1) = sum(alpha(1,:));
alpha(1,:) = alpha(1,:)/c(1);
for t=2:T
    alpha(t,:) = alpha(t-1,:)*A;
    alpha(t,:) = alpha(t,:).*dens(t,:);
    c(t) = sum(alpha(t,:));
    alpha(t,:) = alpha(t,:)/c(t);
end
ll = cumsum(log(c));

% filtre arrière

beta = ones(K,T);
for t=T-1:-1:1
    beta(:,t) = beta(:,t+1).*(dens(t+1,:))';
    beta(:,t) = A*beta(:,t);
    beta(:,t) = beta(:,t)/(alpha(t,:)*beta(:,t));
end

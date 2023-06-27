function Ur=directsolver(A,GL,GR,s,PhiL,PhiR,ui,F,K)
         n=length(A)-1;
         sigma=zeros(n*K,1);
         Ur=zeros(2,n*K);
         P=zeros(n*K,n*K);
         for i=1:n
             sigma((i-1)*K+1:i*K)=Ch1([A(i),A(i+1)],K);
         end
         for j=1:n
             for i=1:K
                 t=(j-1)*K+i;
                 E=zeros(K,1);
                 E(i)=1;
                 for T=1:n*K
                     if T<(j-1)*K+1
                        P(T,t)=PhiR(sigma(T))*GR(sigma(t))*...
                            RIntegral1(A(j),A(j),A(j+1),E,K);
                     elseif T>j*K
                            P(T,t)=PhiL(sigma(T))*GL(sigma(t))*...
                                LIntegral1(A(j),A(j+1),A(j+1),E,K);
                     else 
                         P(T,t)=PhiR(sigma(T))*GR(sigma(t))*...
                             RIntegral1(A(j),sigma(T),A(j+1),E,K)...
                         +PhiL(sigma(T))*GL(sigma(t))*...
                         LIntegral1(A(j),sigma(T),A(j+1),E,K);
                     end
                 end
             end
         end
         P=P+eye(n*K);
         Z=sigma;
         for i=1:n*K
             Z(i)=F(sigma(i));
         end
         Sigma=P\Z;
         JL=zeros(1,n+1);
         JR=zeros(1,n+1);
         EL=zeros(n,K);
         ER=zeros(n,K);
         for j=1:n
             for i=1:K
                 EL(j,i)=GL(sigma((j-1)*K+i));
                 ER(j,i)=GR(sigma((j-1)*K+i));
             end
         end
         for j=1:n
             JL(j+1)=JL(j)+LIntegral1(A(j),A(j+1),A(j+1),EL(j,:)'.*...
                 Sigma((j-1)*K+1:j*K),K);
             JR(n+1-j)=JR(n+2-j)+RIntegral1(A(n+1-j),A(n+1-j),A(n+2-j),...
                 ER(n+1-j,:)'.*Sigma((n-j)*K+1:(n-j+1)*K),K);
         end
         for j=1:n
             for i=1:K
                 n=(j-1)*K+i;
                 Ur(2,n)=ui(sigma(n))+ER(j,i)/s*(JL(j)+ ...
                     LIntegral1(A(j),sigma(n),A(j+1),EL(j,:).*...
                     Sigma((j-1)*K+1:j*K),K))+...
                     EL(j,i)/s*(JR(j+1)+RIntegral1(A(j),sigma(n),...
                     A(j+1),ER(j,:).*Sigma((j-1)*K+1:j*K),K));
             end
         end
         Ur(1,:)=sigma';
end







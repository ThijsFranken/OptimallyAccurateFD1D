%_____________________3-PT HOMOGENEOUS_______________________
%____________________________________________________________

xm = 0;
deltax = 1;
NX=101; % the number of nodes

%xm1 = xm + deltax;
ndis = 100;                               % for discretizing steps

trialfunction=2; % 1 for linear spline, 2 for sinc
npTF=7;          % number of points in scheme (3, 5, 7)

% definition of T and H matrix

T=zeros(npTF,NX);
H1=zeros(npTF,NX);
H2=zeros(npTF,NX);


%_______________________TRIAL FUNCTIONS_______________________

if trialfunction == 1;                      %enter 1 for spline interpolation
    for i = -ndis:0;                            %for xm-1 < x < xm
        x = xm + (i/ndis)*deltax;
        phim(i+1+ndis)=(x+deltax)/deltax;
        phimderiv(i+1+ndis)=1./deltax;
    end

    for i = 0:ndis;                             %for xm < x < xm+1
        x = xm + (i/ndis)*deltax;
        phim(i+1+ndis)=(-x+deltax)/deltax;
        phimderiv(i+1+ndis)=-1./deltax;
    end

elseif trialfunction == 2;                  %enter 2 for sinc interpolation
    
    if (npTF == 3);
        
        for i = -ndis:ndis;                       
            x = xm+(i/ndis)*deltax;
            xx= (x-xm)/deltax;
            phim(i+1+ndis) = sinc(xx);
        if i == 0;
            phimderiv(i+1+ndis)=0;
        else
            %phimderiv(i+1+ndis) = (phim(i+2+ndis)-phim(i+1+ndis))/(deltax/ndis);
            phimderiv(i+1+ndis) = pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx); %Analytical solution
        end
        end
    
    elseif (npTF == 5);
        for i = -2*ndis:2*ndis;                       
            x = xm+(i/ndis)*deltax;
            xx= (x-xm)/deltax;
            phim(i+1+2*ndis) = sinc(xx);
        if i == 0;
           phimderiv(i+1+2*ndis)=0;
        else
            phimderiv(i+1+2*ndis) = pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx);
    
        end
        end
        
     elseif (npTF == 7);
        for i = -3*ndis:3*ndis;                       
            x = xm+(i/ndis)*deltax;
            xx= (x-xm)/deltax;
            phim(i+1+3*ndis) = sinc(xx);
        if i == 0;
           phimderiv(i+1+3*ndis)=0;
        else
            phimderiv(i+1+3*ndis) = pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx);
    
        end
        end
    end
    
end

%_______________________ELASTIC PROPERTIES________________________


% linearly increasing rho: rho=a_rho*X +b_rho
a_rho=0.;         % a_rho=0 for homogeneous
b_rho=1.;

% linearly increasing mu: mu=a_mu*X +b_mu
a_mu=0.;          % a_mu=0 for homogeneous
b_mu=1.;

%____________________COMPUTATION OF T AND H MATRICES______________________


% for a 3-point scheme:
if (npTF==3);                    
    
    for m=2:NX-2;               %first and last row affected by B.C.
    xm=deltax*(m-2);

    x=xm+(i/ndis)*deltax;
    rho=a_rho*x+b_rho;
    mu=a_mu*x+b_mu;

        for n=1:m+1;           % tridiagonalizing matrix
            
        % in fact for the end members m=1 and m=NX, you have to truncate
        % (ignore) the trial functions but I let you do that by yourself NF
        
        T(m,n)=0.;
        H1(m,n)=0.;
        H2(m,n)=0.;
        
            if (n == m-1)                %calculating the Tm(m-1), Hm(m-1)
        
                for i = 0:ndis-1;     
                
                  T(m,n)=T(m,n) + phim(i+1+ndis)*rho*phim(i+1)*deltax/ndis;
                  H1(m,n)=H1(m,n) + phim(i+1+ndis)*mu*phim(i+1)*deltax/ndis;
                  H2(m,n)=H2(m,n) + phimderiv(i+1+ndis)*mu*phimderiv(i+1)*deltax/ndis;
              
                end
        
             elseif(n==m)            %calculating the Tm(m), Hm(m)
        
                for i=-ndis:ndis-1;
                
                    T(m,n)=T(m,n) + phim(i+1+ndis)*rho*phim(i+1+ndis)*deltax/ndis;
                    H1(m,n)=H1(m,n) + phim(i+1+ndis)*mu*phim(i+1+ndis)*deltax/ndis;
                    H2(m,n)=H2(m,n) + phimderiv(i+1+ndis)*mu*phimderiv(i+1+ndis)*deltax/ndis;
                end
        
             elseif (n == m+1)         %calculating the Tm(m+1), Hm(m+1)
    
                for i = 0:ndis-1;
 
                T(m,n)=T(m,n) + phim(i+1+ndis)*rho*phim(i+1)*deltax/ndis;
                H1(m,n)=H1(m,n) + phim(i+1+ndis)*mu*phim(i+1)*deltax/ndis;
                H2(m,n)=H2(m,n) + phimderiv(i+1+ndis)*mu*phimderiv(i+1)*deltax/ndis;
                
                end
            end
        end
    end
 
    
% for a 5-point scheme:
    elseif (npTF==5);
    
    for m=3:NX-3;                      %first and last 2 rows affected by B.C.
    xm=deltax*(m-3);

    x=xm+(i/ndis)*deltax;
    rho=a_rho*x+b_rho;
    mu=a_mu*x+b_mu;

        for n=1:m+2;                    % tridiagonalizing matrix

        T(m,n)=0.;
        H1(m,n)=0.;
        H2(m,n)=0.;
        
            if (n == m-2)                %calculating the Tm(m-2), Hm(m-2)
        
                for i = 0:2*ndis-1;     
                
                  T(m,n)=T(m,n) + phim(i+1+2*ndis)*rho*phim(i+1)*deltax/ndis;
                  H1(m,n)=H1(m,n) + phim(i+1+2*ndis)*mu*phim(i+1)*deltax/ndis;
                  H2(m,n)=H2(m,n) + phimderiv(i+1+2*ndis)*mu*phimderiv(i+1)*deltax/ndis;
              
                end
        
             elseif(n == m-1)            %calculating the Tm(m-1), Hm(m-1)
        
                for i=-ndis:2*ndis-1;
                
                    T(m,n)=T(m,n) + phim(i+1+2*ndis)*rho*phim(i+1+ndis)*deltax/ndis;
                    H1(m,n)=H1(m,n) + phim(i+1+2*ndis)*mu*phim(i+1+ndis)*deltax/ndis;
                    H2(m,n)=H2(m,n) + phimderiv(i+1+2*ndis)*mu*phimderiv(i+1+ndis)*deltax/ndis;
                end

             elseif (n == m)         %calculating the Tm(m), Hm(m)

                for i = -2*ndis:2*ndis-1;
 
                T(m,n)=T(m,n) + phim(i+1+2*ndis)*rho*phim(i+1+2*ndis)*deltax/ndis;
                H1(m,n)=H1(m,n) + phim(i+1+2*ndis)*mu*phim(i+1+2*ndis)*deltax/ndis;
                H2(m,n)=H2(m,n) + phimderiv(i+1+2*ndis)*mu*phimderiv(i+1+2*ndis)*deltax/ndis;
                
                end
                
             elseif(n == m+1)            %calculating the Tm(m+1), Hm(m+1)
      
                for i=-ndis:2*ndis-1;
                
                    T(m,n)=T(m,n) + phim(i+1+2*ndis)*rho*phim(i+1+ndis)*deltax/ndis;
                    H1(m,n)=H1(m,n) + phim(i+1+2*ndis)*mu*phim(i+1+ndis)*deltax/ndis;
                    H2(m,n)=H2(m,n) + phimderiv(i+1+2*ndis)*mu*phimderiv(i+1+ndis)*deltax/ndis;  
                end
                
              elseif (n == m+2)                %calculating the Tm(m+2), Hm(m+2)
                  
                for i = 0:2*ndis-1;     
                
                  T(m,n)=T(m,n) + phim(i+1+2*ndis)*rho*phim(i+1)*deltax/ndis;
                  H1(m,n)=H1(m,n) + phim(i+1+2*ndis)*mu*phim(i+1)*deltax/ndis;
                  H2(m,n)=H2(m,n) + phimderiv(i+1+2*ndis)*mu*phimderiv(i+1)*deltax/ndis;
              
                end
            end
        end
    end
 
    
%for a 7-point scheme:
elseif (npTF==7);   

    for m=4:NX-4;                     %first and last 3 rows affected by B.C.
    xm=deltax*(m-4);

    x=xm+(i/ndis)*deltax;
    rho=a_rho*x+b_rho;
    mu=a_mu*x+b_mu;

        for n=1:m+3;                  % tridiagonalizing matrix

        T(m,n)=0.;
        H1(m,n)=0.;
        H2(m,n)=0.;
        
             if (n == m-3)                %calculating the Tm(m-3), Hm(m-3)
        
                for i = 0:3*ndis-1;     
                
                  T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1)*deltax/ndis;
                  H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1)*deltax/ndis;
                  H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1)*deltax/ndis;
              
                end
        
             elseif(n == m-2)            %calculating the Tm(m-2), Hm(m-2)
        
                for i=-ndis:3*ndis-1;
                
                    T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1+ndis)*deltax/ndis;
                    H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1+ndis)*deltax/ndis;
                    H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1+ndis)*deltax/ndis;
                end
        
             elseif (n == m-1)         %calculating the Tm(m-1), Hm(m-1)

                for i = -2*ndis:3*ndis-1;
 
                T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1+2*ndis)*deltax/ndis;
                H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1+2*ndis)*deltax/ndis;
                H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1+2*ndis)*deltax/ndis;
                
                end
                
             elseif(n == m)            %calculating the Tm(m), Hm(m)
        
                for i=-3*ndis:3*ndis-1;
                
                    T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1+3*ndis)*deltax/ndis;
                    H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1+3*ndis)*deltax/ndis;
                    H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1+3*ndis)*deltax/ndis;  
                end       

             elseif (n == m+1)                %calculating the Tm(m+1), Hm(m+1)
      
                for i = -2*ndis:3*ndis-1;     
                
                  T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1+2*ndis)*deltax/ndis;
                  H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1+2*ndis)*deltax/ndis;
                  H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1+2*ndis)*deltax/ndis;
              
                end              

             elseif(n == m+2)            %calculating the Tm(m+2), Hm(m+2)
 
                for i=-ndis:3*ndis-1;
                
                    T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1+ndis)*deltax/ndis;
                    H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1+ndis)*deltax/ndis;
                    H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1+ndis)*deltax/ndis;
                end
              
             elseif (n == m+3)                %calculating the Tm(m+3), Hm(m+3)

                for i = 0:3*ndis-1;     
                
                  T(m,n)=T(m,n) + phim(i+1+3*ndis)*rho*phim(i+1)*deltax/ndis;
                  H1(m,n)=H1(m,n) + phim(i+1+3*ndis)*mu*phim(i+1)*deltax/ndis;
                  H2(m,n)=H2(m,n) + phimderiv(i+1+3*ndis)*mu*phimderiv(i+1)*deltax/ndis;
              
                end
            end
        end
    end
end









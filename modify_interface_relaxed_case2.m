function [A,f]=modify_interface_relaxed_case2(A, f, u_star, meshes, flag)
%     for interface condition with relaxed parameters;
%     The case II: 1+s1*2; 2+s2*1; 

paras;

dy=meshes.dy;
m=meshes.m;
n=meshes.n;

% s1=kappa/(kappa+nu2*p1);
% s2=kappa/(kappa+nu1*p2);


for i=1:m-1
    temp=(m-1)*(n-1)+i;
    if flag==2
        A(temp, temp) = A(temp, temp)-(nu2/(dy*dy)) * (-0.5*(1-s2)*kappa+nu2/dy)/(0.5*(1-s2)*kappa+nu2/dy);

        f(temp)=f(temp)+(nu2/(dy*dy)) * ((  -(s2*nu1/dy) +1.5*(1-s2)*kappa)*u_star(i) + (  (s2*nu1/dy) - 0.5*(1-s2)*kappa )*u_star(i+m-1) )/ (0.5*(1-s2)*kappa + nu2/dy);   
            
    elseif flag==1
        A(i, i)=A(i, i) - (nu1/(dy*dy)) * (-0.5*(s1-1)*kappa - nu1/dy) / (0.5*(s1-1)*kappa - nu1/dy) ;    

        f(i)=f(i)+ (nu1/(dy*dy)) * (( 1.5*(s1-1)*kappa + s1*nu2/dy )*u_star(temp) + (-0.5*(s1-1)*kappa - s1*nu2/dy)*u_star(temp-m+1)) / (0.5*(s1-1)*kappa - nu1/dy);   
                  
    end
end


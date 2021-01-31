function [A,f]=modify_interface_relaxed(A, f, u_star, meshes, flag)
%     for interface condition with relaxed parameters;

paras;

dy=meshes.dy;
m=meshes.m;
n=meshes.n;


for i=1:m-1
    temp=(m-1)*(n-1)+i;
    if flag==1
        A(i,i) = A(i,i)-(nu1/(dy*dy)) * (-0.5*s1*kappa-nu1/dy)/(0.5*s1*kappa-nu1/dy);

        %         ---- the revised version: approximation        
        f(i)=f(i)+(nu1/(dy*dy)) * ((  ((1-s1)*nu2/dy) +1.5*s1*kappa)*u_star(temp) + ( - ((1-s1)*nu2/dy) - 0.5*s1*kappa )*u_star(temp-m+1) )/ (0.5*s1*kappa - nu1/dy); % The revised version
    elseif flag==2
        A(temp, temp)=A(temp, temp)- (nu2/(dy*dy)) * (0.5*s2*kappa + nu2/dy) / (-0.5*s2*kappa+ nu2/dy) ;

        %         ---- the original version: in draft
        f(temp)=f(temp)+ (nu2/(dy*dy)) * ((   -1.5*s2*kappa - (1-s2)*nu1/dy )*u_star(i) + (0.5*s2*kappa + (1-s2)*nu1/dy)*u_star(i+m-1)) / (-0.5*s2*kappa + nu2/dy);
        
    end
end


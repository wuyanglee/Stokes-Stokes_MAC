function [A,f]=modify_interface(A, f, u_star, nu, kappa, dx,dy, m, n, flag)   
%     for interface condition;
    
    for i=1:m-1
        temp=(m-1)*(n-1)+i;
        if flag==1
            A(i,i)=A(i,i)-(nu*(nu-0.5*dy*kappa))/(dy*dy*(nu+0.5*dy*kappa));
            f(i)=f(i)+( 2*kappa*dy*nu * (2*u_star(temp)-u_star(temp-m+1)) )/(dy*dy*(2*nu+kappa*dy));
%             f(i)=f(i)+( 2*kappa*dy*nu * u_star(temp) )/(dy*dy*(2*nu+kappa*dy));

%             2*u_star(temp)-u_star(temp-m+1)
%             u_star(temp)
        elseif flag==2
            A(temp, temp)=A(temp, temp)-(nu*(nu-0.5*dy*kappa))/(dy*dy*(nu+0.5*dy*kappa));
            f(temp)=f(temp)+(2*kappa*dy*nu* (2*u_star(i) - u_star(i+m-1)) )/(dy*dy*(2*nu+kappa*dy));
%             f(temp)=f(temp)+(2*kappa*dy*nu* u_star(i) )/(dy*dy*(2*nu+kappa*dy));

%             2*u_star(i) - u_star(i+m-1)
%             u_star(i)
    end
end


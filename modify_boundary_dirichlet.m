function [A,f]=modify_boundary_dirichlet(A,f,domain,meshes, Paras, flag)
%     Dirichlet boundary conditions: u~=0 on Gamma

nu1 = Paras.nu1; % the dynamic viscosity. for the PDEs.
nu2 = Paras.nu2;

if flag==1
    nu=nu1;
else
    nu=nu2;
end

dx=meshes.dx;
dy=meshes.dy;
m=meshes.m;
n=meshes.n;

exact_solution;

%% for left boundary condition;
x=domain.left;
if flag ==1
    y=(domain.bottom+dy/2):dy:(domain.top-dy/2);
    u1_left=solution.u11(x,y); % for u11, horizontal direction
else
    y=[(domain.bottom+dy/2):dy:domain.top-dy/2];
    u1_left=solution.u21(x,y);
end
% ------divergence free condition.
for i=1:n
    f(1+(m-1)*(i-1), 1)=f(1+(m-1)*(i-1), 1)+nu*u1_left(i)/(dx*dx);
    f(meshes.DOF_u+meshes.DOF_v+1+(i-1)*m) = f(meshes.DOF_u+meshes.DOF_v+1+(i-1)*m) + u1_left(i)/dx;
end
if flag ==1
    y=domain.bottom+dy:dy:(domain.top-dy);
    u2_left=solution.u12(x,y); % for u12, vertical direction
else
    y=domain.bottom+dy:dy:(domain.top-dy);
    u2_left=solution.u22(x,y);
end
for i=1:1:n-1
    temp=(m-1)*n+(i-1)*m+1;
    f(temp,1)=f(temp, 1)+nu*2*u2_left(i)/(dx*dx);
    A(temp, temp)= A(temp, temp)+nu/(dx*dx);
end
%%     for right boundary condition;
x=domain.right;
if flag==1
    y=(domain.bottom+dy/2):dy:(domain.top-dy/2);
    u1_right=solution.u11(x,y); % for u11, horizontal direction
else
    y=(domain.bottom+dy/2):dy:(domain.top-dy/2);
    u1_right=solution.u21(x,y);
end
%  % ------divergence free condition.
for i=1:n
    temp=i*(m-1);
    f(temp,1)=f(temp, 1)+nu*u1_right(i)/(dx*dx);
    f(meshes.DOF_u+meshes.DOF_v+i*m) = f(meshes.DOF_u+meshes.DOF_v+i*m) - u1_right(i)/dx;
end
if flag==1
    y=(domain.bottom+dy):dy:(domain.top-dy);
    u2_right=solution.u12(x,y); % for u12, vertical direction
else
    y=(domain.bottom+dy):dy:(domain.top-dy);
    u2_right=solution.u22(x,y);
end
for i=1:n-1
    temp=(m-1)*n+i*m;
    A(temp,temp)=A(temp,temp)+nu/(dx*dx);
    f(temp,1)=f(temp,1)+nu*2*u2_right(i)/(dx*dx);
end
%%  -----       for top/bottom  boundary condition;
if flag==1
    y=domain.top;
    x=(domain.left+dx):dx:(domain.right-dx);
    u1_top=solution.u11(x,y); % for u11, horizontal direction
    
    x=(domain.left+dx/2):dx:(domain.right-dx/2);
    u2_top=solution.u12(x,y); % for u12, vertical direction
    
    for i=1:m-1
        temp=(m-1)*(n-1)+i;
        f(temp,1)=f(temp,1)+2*nu*u1_top(i)/(dy*dy);
        A(temp,temp)=A(temp,temp)+nu/(dy*dy);
    end
    %      % ------divergence free condition.
    for i=1:m
        temp=(m-1)*n+m*(n-2)+i;
        f(temp)=f(temp)+nu*u2_top(i)/(dy*dy);
        f(meshes.DOF_u+meshes.DOF_v+i+(n-1)*m) = f(meshes.DOF_u+meshes.DOF_v+i+(n-1)*m) - u2_top(i)/dy;
    end
elseif flag==2
    y=domain.bottom;
    x=(domain.left+dx):dx:(domain.right-dx);
    u1_top=solution.u21(x,y); % for u21, horizontal direction
    
    x=(domain.left+dx/2):dx:(domain.right-dx/2);
    u2_top=solution.u22(x,y); % for u22, vertical direction
    
    for i=1:m-1
        temp=i;
        f(temp,1)=f(temp,1)+2*nu*u1_top(i)/(dy*dy);
        A(temp,temp)=A(temp,temp)+nu/(dy*dy);
    end
    %       % ------ divergence free condition.
    for i=1:m
        temp=(m-1)*n+i;
        f(temp)=f(temp)+nu*u2_top(i)/(dy*dy);
        f(meshes.DOF_u+meshes.DOF_v+i) = f(meshes.DOF_u+meshes.DOF_v+i) + u2_top(i)/dy;
    end
end
%% ------------- ---test---
% if flag==2
%     y=0;
%     x=dx:dx:(1-dx);
%     u1_top=solution.u21(x,y); % for u11, horizontal direction
% 
%     x=(dx/2):dx:(1-dx/2);
%     u2_top=solution.u22(x,y); % for u12, vertical direction
% 
%     for i=1:m-1
%         temp=(m-1)*(n-1)+i;
%         f(temp,1)=f(temp,1)+2*nu*u1_top(i)/(dy*dy);
%         A(temp,temp)=A(temp,temp)+nu/(dy*dy);
%     end
%     for i=1:m
%         temp=(m-1)*n+m*(n-2)+i;
%         f(temp)=f(temp)+nu*u2_top(i)/(dy*dy);
%         f(meshes.DOF_u+meshes.DOF_v+i+(n-1)*m) = f(meshes.DOF_u+meshes.DOF_v+i+(n-1)*m) - u2_top(i)/dy;  % divergence free condition.
%     end
% elseif flag==1
%     y=0;
%     x=dx:dx:(1-dx);
%     u1_top=solution.u11(x,y); % for u21, horizontal direction
% 
%     x=(dx/2):dx:(1-dx/2);
%     u2_top=solution.u12(x,y); % for u22, vertical direction
% 
%     for i=1:m-1
%         temp=i;
%         f(temp,1)=f(temp,1)+2*nu*u1_top(i)/(dy*dy);
%         A(temp,temp)=A(temp,temp)+nu/(dy*dy);
%     end
%     for i=1:m
%         temp=(m-1)*n+i;
%         f(temp)=f(temp)+nu*u2_top(i)/(dy*dy);
%         f(meshes.DOF_u+meshes.DOF_v+i) = f(meshes.DOF_u+meshes.DOF_v+i) + u2_top(i)/dy;  % divergence free condition.
%     end
% 
% end


end


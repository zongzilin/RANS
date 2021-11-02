clear;
clc;
close all

nhx = 18 + 2;
nhy = 18 + 2;

nx =  nhx - 2;
ny =  nhy - 2;

hx = 0.5;
hy = 0.006250000000000;

as = zeros(nx*ny,nx*ny);


for i = 1:nx*ny
    for j = 1:nx*ny
        ap(i,j) = 2/hx^2 + 2/hy^2 ;
    end
end


for i = nx+1:nx*ny
    for j = 1:nx*ny-nx
        as(i,j) = -1/hy^2;
    end
end

an = transpose(as);

for i = nx+2:nx*ny-1
    for j = ny+1:nx*ny - 1
        ae(i,j) = -1/hx^2;
    end
end

ae(nx*ny-nx-2:nx*ny,nx*ny-ny-1:nx*ny) = 0;
ae(nx+1:nx:nx*ny,ny:ny:nx*ny) = 0;

aw = transpose(ae);
% figure(1)
% spy(as);
% 
% figure(2)
% spy(an)

figure(3)
spy(ap)


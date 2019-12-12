
I = 13:0.5:45;
L = 1:0.5:10;


[i,l] = meshgrid(I,L);
f = ((i.^2).*l.*(0.0000013))./((0.0005.^2).*pi);

surf(i,l,f)
xlabel('Current(A)')
ylabel('Element Length(m)')
zlabel('Power Output')


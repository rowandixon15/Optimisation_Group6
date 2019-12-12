
R = 0.0005:0.00001:0.0009;
P = 0.000001:0.0000005:0.000009;


[r,p] = meshgrid(R,P);
f = ((45.^2).*10.*(p))./((r.^2).*pi);

surf(r,p,f)
xlabel('Radius(m))')
ylabel('resistivity(Ohm meter)')
zlabel('Power Output')


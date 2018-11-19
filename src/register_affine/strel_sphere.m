function se = strel_sphere(sz)
[x1,y1,z1] = ndgrid(-ceil(sz):ceil(sz));
se = (sqrt(x1.^2 + y1.^2 + z1.^2) <=sz);
end

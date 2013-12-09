nQ = [1 1 1];
v1 = [2,9,4]';
v2 = [2,9,4]';

clf;

[rot_init, diffrot_init] = diffqn(nQ);

for j=1:3
for i=1:1000

  nQ(j) = (i/100)-5;
  [rot, diffrot] = diffqn(nQ);

  var_vec(i) = (i/100)-5;
  res(i) = dot((rot*v1+rot_init*v2),diffrot(1:3,1:3)*v1) / norm(rot*v1+rot_init*v2);
endfor

%back to initial value for comparable results
nQ = [1 1 1];

plot(var_vec, res);
hold on;

endfor




close all; clear

example_simple_inverse;
close all

assert(abs(st-3.880139275426690)<1e-12);

close all; clear

example_elastic_inverse;
close all

assert(abs(errorG-0.105029741035073)<1e-12);
assert(abs(errorK-0.009002062785570)<1e-12);

disp('All tests passed!')

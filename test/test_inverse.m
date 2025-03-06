close all; clear

example_simple_inverse;
close all

assert(abs(st-3.941105683019185)<1e-12);

close all; clear

example_elastic_inverse;
close all

assert(abs(errorG-0.425355069811006)<1e-12);
assert(abs(errorK-0.016597556081363)<1e-12);

disp('All tests passed!')

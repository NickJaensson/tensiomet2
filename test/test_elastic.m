close all; clear

addpath('../app/')

gen_single_drop_elastic;

% compare to old values (gen-pendant-drop before refactoring:
eps2 = 1e-12;
assert ( abs(volume-12.7999999999999) < eps2 );
assert ( abs(area-24.3099753701003) < eps2 );
assert ( abs(p0-3.06593554364227) < eps2 );
assert ( abs(max(taus)-3.751693556096941) < eps2 );
assert ( abs(max(taup)-4.000492342729165) < eps2 );

disp('All tests passed!')
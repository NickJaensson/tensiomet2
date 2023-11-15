close all; clear

addpath('../app/')

gen_single_drop_elastic;

close all

% disp ( abs(volume-12.8000000000145) );
% disp ( abs(area-22.5156483902244) );
% disp ( abs(p0-2.02056164104927) );
% disp ( abs(max(sigmas)-3.33958761227839) );
% disp ( abs(max(sigmap)-3.86864491619739) );

% compare to old values (gen-pendant-drop before refactoring:
eps2 = 1e-10; 
assert ( abs(volume-12.8000000000001) < eps2 );
assert ( abs(area-22.5156483902096) < eps2 );
assert ( abs(itervars.p0-2.02056164104124) < eps2 );
assert ( abs(max(itervars.sigmas)-3.33958761227925) < eps2 );
assert ( abs(max(itervars.sigmap)-3.86864491619756) < eps2 );

disp('All tests passed!')
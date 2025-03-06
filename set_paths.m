% run this script from the gen-pendant-drop/ folder
addpath('app/')
addpath('app/example_parameters/')
addpath('src/')
addpath('test/')
addpath('symbolic/')

if isfolder('chebfun/')
    addpath('chebfun/')
elseif isfolder('../chebfun/')
    addpath('../chebfun/')
else
    error('Chebfun library not found')
end

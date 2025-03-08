% run this script from the gen-pendant-drop/ folder
addpath('examples/')
addpath('examples/example_parameters/')
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

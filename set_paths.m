% run this script from the tensiomet2/ folder
addpath('examples/')
addpath('examples/example_parameters/')
addpath('src/')
addpath('test/')
addpath('symbolic/')

addpath('manuscript_results/')
addpath('manuscript_results/example_parameters/')
addpath('manuscript_results/src_pp/')

if isfolder('chebfun/')
    addpath('chebfun/')
elseif isfolder('../chebfun/')
    addpath('../chebfun/')
else
    error('Chebfun library not found')
end

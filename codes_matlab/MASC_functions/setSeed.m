function setSeed(seed)
% set seed - old version 
% global RNDN_STATE  RND_STATE
% RNDN_STATE = randn('state');
% randn('state',seed);
% RND_STATE = rand('state');
% %rand('state',seed);
% rand('twister',seed);

    rng(seed);

end

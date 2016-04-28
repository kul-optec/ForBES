function cache = CacheInit(prob, x, gam)

cache.prob = prob;
cache.x = x;
cache.gam = gam;

cache.flagInit = 1;         
cache.flagEvalf = 0;
cache.flagGradStep = 0;
cache.flagProxGradStep = 0;
cache.flagFBE = 0;          
cache.flagGradFBE = 0;      
cache.flagLineSearch = 0;   

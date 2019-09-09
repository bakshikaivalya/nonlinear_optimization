function[state,options,optchanged] = gaoutputfun(options, state, flag)
 optchanged = false;

 % Note: data is stored in index Generation+1. E.g. min f for
 % generation 0, the initial population, is stored in fmin(1)
 evalin('base',['fmin(' int2str(state.Generation+1) ')=' num2str(min(state.Score)) ';'])
 evalin('base',['fmax(' int2str(state.Generation+1) ')=' num2str(max(state.Score)) ';'])
 evalin('base',['favg(' int2str(state.Generation+1) ')=' num2str(mean(state.Score)) ';'])
 D = pdist(state.Population);
 Davg = mean(D);
 evalin('base',['avgdist(' int2str(state.Generation+1) ')=' num2str(Davg) ';'])
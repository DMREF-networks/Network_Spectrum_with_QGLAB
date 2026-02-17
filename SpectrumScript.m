function []=SpectrumScript(conf, N, tess, bdry, name, equi, p1, p2)
arguments
	conf (1,1) string
	N (1,1) string
	tess (1,1) string 
	bdry (1,1) string
	name (1,1) string
	equi (1,1) string
	p1 (1,1) string = ''
	p2 (1,1) string = ''
end

proj = openProject('QGObject.prj');

% Building Quantum graphs from an edge list:

if conf == 'URL' | conf == 'A2'
        filetag = conf + '_a' + p1 + '_' + tess + '_N' + N + '_' + bdry + '_' + name;
elseif conf == 'Ti'
	filetag = conf + '_' + p1 + '_' + tess + '_N' + N + '_' + bdry + '_' + name;
else
        filetag = conf + '_' + p2 + '_' + tess + '_N' + N + '_' + bdry + '_' + name;
end


edgelist = load('../configs/edgelists/edges_'+filetag+'.txt');
s = (1+edgelist(:,1))';
t = (1+edgelist(:,2))';
[tt,It] = sort(t);
ss = s(It);
[ss1,Is] = sort(ss);
tt1 = tt(Is);
LVec = edgelist(:,3)';
LVect = LVec(It);
LVecs = LVect(Is);
nX = round(10*LVec);
nXt = nX(It);
nXs = nXt(Is);


if equi == 'False'
	Phi=quantumGraph(ss1,tt1,LVecs,'Discretization','Uniform','nxVec',nXs);  %'RobinCoeff',robinCoeff,
else
	Phi=quantumGraph(ss1,tt1,1,'Discretization','Uniform','nxVec',10);
end

nodes = load('../configs/configs/vertlocs_'+filetag+'.txt');

plotCoordFcn=@(G)plotCoordFcnFromNodes(G,nodes);
Phi.addPlotCoords(plotCoordFcn);
Phi.plot('layout')


[V,lambda]=eigs(Phi,1500);
[singles,doubles,~]=separateEigs(lambda); 

lambda

save("Esys_"+filetag+".mat","lambda", "V", "Phi")
%save("Lamb_"+filetag+".mat","lambda")

function write_graphfile(A,names,p,gfile,nfile,dag);

% Creates the graph data file (nodes, edges) and the node names 
%  file for input to graph drawing ... such as GraphExplor

% dag=0 for undirected graphs
% dag=1 for DAGs 

%write names/description file

fid = fopen(nfile,'w');
for (i=1:p)
    fprintf(fid, '%i\t', i);
    fprintf(fid, '%s\n', names(i,:));
end
fclose(fid);


%write graph or dag file


fid = fopen(gfile,'w');
for (i = 1:p)
    a=(i+1):p; 
    J=find(abs(A(i,a))>0);
    for (j=1:length(J))
        fprintf(fid, '%i\t', i);
        fprintf(fid, '%i\t', a(J(j)));
        fprintf(fid, '%6.3f\t', (1-dag)*A(i,a(J(j))));
        fprintf(fid, '%i\t', 0);
        fprintf(fid, '%i\n', dag);
    end
end
fclose(fid);
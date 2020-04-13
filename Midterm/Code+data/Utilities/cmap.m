
function cmap(map)
if (nargin>0)
   if (length(map)==2)
       if ( (((map(1)=='r' & map(2) =='g') |(map(1)=='g' & map(2) =='r'))|...
            ((map(1)=='R' & map(2) =='G') |(map(1)=='G' & map(2) =='R'))))
            map = zeros(255,3);
            map(:,1) = [zeros(1,127), (0:1:127)/127]';
            map(:,2) = [(127:-1:0)/127, zeros(1,127)]';
            %	map=zeros(64,3); map(:,2)=flipud(((0:1:63)/63)'); map(:,1)=((0:1:63)/63)'; 
       else
           map=[[1 1 1];[0 0 0]];
       end;
   end; 
   else
   map='jet';
end
colormap(map)







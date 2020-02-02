function v = portvar(a,k,q,f,K,r);
    Ka=zeros(q,q); fa=zeros(q,1);
    for i=1:k
        Ka=Ka+a(i)*K(:,:,i); fa=fa+a(i)*K(:,:,i)*f(:,i); 
    end
    Va=inv(Ka); fa=Va*fa; 
    l = ones(q,1); h = Ka*l; e = (l'*h)*(fa'*Ka*fa)-(fa'*h)^2; 
    g = (l*r-fa)/e; z = -fa'*Ka*g;  u = h'*g; 
    w = Ka*(u*fa+z*l); 
    v=w'*Va*w; 
end


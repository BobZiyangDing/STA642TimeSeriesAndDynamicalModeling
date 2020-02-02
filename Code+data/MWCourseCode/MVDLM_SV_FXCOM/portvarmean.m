function v = portvarmean(a,k,q,f,V,r);
    Va=zeros(q,q); fa=zeros(q,1);
    for i=1:k
        Va=Va+a(i)*V(:,:,i); fa=fa+a(i)*f(:,i); 
    end
    for i=1:k
        Va=Va+a(i)*(f(:,i)-fa)*(f(:,i)-fa)';
    end
    Ka=inv(Va); 
    l = ones(q,1); h = Ka*l; e = (l'*h)*(fa'*Ka*fa)-(fa'*h)^2; 
    g = (l*r-fa)/e; z = -fa'*Ka*g;  u = h'*g; 
    w = Ka*(u*fa+z*l); 
    v=w'*Va*w; 
end


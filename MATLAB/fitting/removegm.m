function objnew = removegm(obj,numremove)
% remove gaussians from the top numremove simlar pairs and return the new gmm obj

index = 1:obj.NumComponents;
simarray = zeros(obj.NumComponents * (obj.NumComponents-1 )/ 2 , 3);
count = 1;
% calculate pairwise similarity
for i = 1:obj.NumComponents-1
    for j = i+1:obj.NumComponents
        simarray(count,1) = i;
        simarray(count,2) = j;
        simarray(count,3) = gmpairsimilarity(3,obj.mu(i,:),obj.mu(j,:), obj.Sigma(:,:,i),obj.Sigma(:,:,j));
        count = count + 1;
    end
end

sortedsimarray = sortrows(simarray, -3);
removelist = [];
pos = 1;
% from the top numremove pairs, remove the one with low weight
while length(removelist) < numremove
    g1 = sortedsimarray(pos,1);
    g2 = sortedsimarray(pos,2);
    % if one gaussian in the pair has already been removed, then skip
    if ismember(g1,index) && ismember(g2,index)
        if obj.ComponentProportion(g1) <=  obj.ComponentProportion(g2)
            % remove from index list
            index = index(index~=g1);
            % write to remove list so remove from weight, mu, sigma at once
            removelist = [removelist g1];
        else
            index = index(index~=g2);
            removelist = [removelist g2];
        end
    end
    pos = pos + 1;
end

P = obj.ComponentProportion(index);
P = P / sum(P);
Mu = obj.mu(index,:);
Sigma = obj.Sigma(:,:,index);

objnew = gmdistribution(Mu,Sigma,P);

end
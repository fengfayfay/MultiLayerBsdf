function part = make_partition(part, crit, opt, W, M, R)

% function that returns a partition of the data based on a criterion
%
% part = make_partition(part, crit, opt, W, M, R)
%
% part - current partition
% crit - heuristic that decides the new partition
%        1 = 'opt' levels deeper
%        2 = expand the 'opt'*length(partition) best
%        3 = expand full partition
% opt  - (optional) options
% W, M, R - parameters of the current mixture (only needed for crit=2)
%
% part - array containing nodes that define the new partition
%
% Jan Nunnink, 2003

if crit==1
    for j=1:opt
        start_len = length(part);
        len = start_len;
        for i=1:start_len
            if part{i}.type=='node'
                len=len+1;
                part{len}=part{i}.right;
                part{i}=part{i}.left;
            end
        end
    end
end

if crit==2
    lp = length(part);
    nr = fix(lp*opt);
    for i=1:lp						% compute loglikelihood increasements
        if part{i}.type=='node'
            a = log(exp(opt_approx(part{i},M,R))'*W)*part{i}.numpoints;
            al = log(exp(opt_approx(part{i}.left,M,R))'*W)*part{i}.left.numpoints;
            ar = log(exp(opt_approx(part{i}.right,M,R))'*W)*part{i}.right.numpoints;
            ac=al+ar;
            d(i) = a/ac-1;
        end
    end
    [tmp, ind] = sort(-d);		% sort increasements
    len=lp;
    for i=1:nr					% expand the nr best increasements
        if part{ind(i)}.type=='node'
            len=len+1;
            part{len}=part{ind(i)}.right;
            part{ind(i)}=part{ind(i)}.left;
        end
    end
end

if crit==3
    f=1;
    while(f)
    	f=0;
        start_len = length(part);
        len = start_len;
        for i=1:start_len
            if part{i}.type=='node'
            	f=1;
                len=len+1;
                part{len}=part{i}.right;
                part{i}=part{i}.left;
            end
        end
    end
end





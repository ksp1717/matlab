function [Cost] = constrainedtsp(dtMatrix,Limit)
%
% [Cost] = constrainedtsp(Weights,Limit)
%
% A function for finding a feasible solution 
% of the constrained TSP, with constraint on path length.
% The solution may or may not be optimum.
%
% Cost returns the cost of the minimum tour 
% and the minimum tour except the last path.
%
% Written by: Shyamprasad Konduri

% Initialize Variables
[m,n] = size(dtMatrix);
k = 1;

for i = 2:1+Limit
    A = dtMatrix((2:m),[1:i-1,i+1:n]);
    % cost matrix [cost, from city, to city]
    C(k,1:3) = [(dtMatrix(1,i)+sum(min(A,[],2))),1,i];
    k = k + 1;
end
while(1)
    j = [];
    % find minimum of cost column and get the row
    [r,c] = find(C == min(C(:,1)));
    % find the column of the last non zero element in the row
    l = length(r);
    col = find(C(r(l),:),1,'last');
    % select current city based the current branch ending
    city = C(r(l),col-1)+1;
    % update tour matrix with min tour
    Tour=[];
    for lgt = 1:(col-1)/2
        Tour(lgt,:) = C(r(l),[2*lgt,2*lgt+1]);
        j = [j,Tour(lgt,2)];
    end
    % calculate costs from current city the bounds of for take
    % care of the length of path constraint if none of the paths
    % for the current minimum obey the constraint the row is 
    % deleted at the end of the for and the algorithm moves to
    % to the next minimum solution
    for i = (city-Limit):(city+Limit)
         % subtour elimination constraint
         go = subtour(city,i,Tour);
         costtemp = 0;
         if (i > 0 && i ~= city && isempty(find(j==i, 1)) && i<=m && go==1)
             % reduced matrix by removing selected tour rows and columns
             A = reduce(dtMatrix,Tour(:,2),m,0);
             %calculate cost upto the present branch
             for p = 1 : (city-1)
                 costtemp = costtemp + dtMatrix(p,Tour(p,2));
             end
             % cost matrix row = [cost, last branch tour, present tour]
             C(k,1:2*city+1) = [costtemp+sum(min(A,[],2)),C(r(l),2:col),...
                 city,i];
             k = k + 1;
         end
    end
    % remove the minimum cost from C
    C = reduce(C,r(l),m,1);
    k = k - 1;
    % check if all the cities are toured
    rowno = length(C(:,1));
    collgt = length(C(rowno,:));
    value = C(rowno,collgt);
    % Verify the constraint for the final two paths
    % If violated delete the correponding tours
    if (collgt == 2*n-1 && value ~= 0)
        no=1;
        while (no <= n)
            if(isempty(find(Tour(:,2)==no,1)))
                lastcity = no;
                break;
            end
            no = no + 1;
        end
        if ((value-lastcity) > Limit)
            C = reduce(C,rowno,m,1);
            k = k-1;
        else
            break;
        end
    end
end
    % return the cost function along with the tour
    Cost = C(length(C(:,1)),:);
end

function [R] = reduce(M,j,m,k)
% Sub-Function to reduce the matrices by removing the given rows
if (k == 1)
   R = M([1:j-1,j+1:length(M(:,1))],:);
else
    i = length(j);
    R = M(i+1:m,:);
    R(:,j)=[];
end
end

function [g] = subtour(cc,nc,t)
% Sub-Function to check if the current tour is a subtour
if (isempty(find(t(:,1)==nc,1)))
    g = 1;
else
    row = find(t(:,1)==nc,1);
    while(1)
        ncity = t(row,2);
        if (ncity == cc)
            g = 0;
            break;
        end
        row = find(t(:,1)==ncity,1);
        if (isempty(row))
            g = 1;
            break;
        end
    end
end
end
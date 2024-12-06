




function x=adapquad(a,b,tol)
%skilgreina hér fallið f(t) sem er verið að heilda

persistent intervals
if isempty(intervals)
    intervals = zeros(0,2); 
end

if nargin == 0
    fprintf('Total intervals used: %d\n', size(intervals, 1));
    intervals = zeros(0,2);
    x = [];
    return
end

c=(a+b)/2;
sab=simpson(a,b,1);sac=simpson(a,c,1);scb=simpson(c,b,1);

% tjékkar fyrir duplicates
if ~ismember([a, b], intervals, 'rows')
    intervals = [intervals; a, b];
end

if abs(sab-sac-scb)<10*tol 
    x=sac+scb;
else
    x=adapquad(a,c,tol/2)+adapquad(c,b,tol/2);
end
end

function s=trap(f,a,b)
s=(f(a)+f(b))*(b-a)/2;
end
function [CC,I]=differ(A,B)
[CC,I] = setdiff(A,B,'stable');
end

function result = cmult(a,b)
%     var result = new Array();
    result(1) = a(1)*b(1) - a(2)*b(2);
    result(2) = a(2)*b(1) + a(1)*b(2);
%     return result;
end

% function cmult(a,b){
%     var result = new Array();
%     result[0] = a[0]*b[0] - a[1]*b[1];
%     result[1] = a[1]*b[0] + a[0]*b[1];
%     return result;
% }
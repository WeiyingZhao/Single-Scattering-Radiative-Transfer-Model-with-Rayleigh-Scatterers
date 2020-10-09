function result = cdiv(a,b)
%     var result = new Array();

%     var c=[];
    c(1)=b(1);
    c(2)=-b(2);
    
%     var d=[];
    d=cmult(a,c);

%     var denom;
    denom = b(1)*b(1) + b(2)*b(2);

    result(1) = d(1)/denom;
    result(2) = d(2)/denom;
%     return result;
end
% 
% function cdiv(a,b){
%     var result = new Array();
% 
%     var c=[];
%     c[0]=b[0];
%     c[1]=-b[1];
%     
%     var d=[];
%     d=cmult(a,c);
% 
%     var denom;
%     denom = b[0]*b[0] + b[1]*b[1];
% 
%     result[0] = d[0]/denom;
%     result[1] = d[1]/denom;
%     return result;
% }
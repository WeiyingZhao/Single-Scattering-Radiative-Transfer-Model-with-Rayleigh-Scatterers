function result = csqrt(a)
%     var r, theta, rootr;

%     var result = new Array();
    
    r = sqrt(a(1)*a(1) + a(2)*a(2));
    theta = atan2(a(2),a(1));
    rootr = sqrt(r);

    result(1) = rootr*cos(0.5*theta);
    result(2) = rootr*sin(0.5*theta);
%     return result;

end

% function csqrt(a){
%     var r, theta, rootr;
% 
%     var result = new Array();
%     
%     r = Math.sqrt(a[0]*a[0] + a[1]*a[1]);
%     theta = Math.atan2(a[1],a[0]);
%     rootr = Math.sqrt(r);
% 
%     result[0] = rootr*Math.cos(0.5*theta);
%     result[1] = rootr*Math.sin(0.5*theta);
%     return result;
% 
% }
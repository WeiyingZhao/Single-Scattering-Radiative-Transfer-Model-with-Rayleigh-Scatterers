% //************************************************************
function result = refl_coef(the1, eps1, eps2) 
%     // calculates the v and h-polarized reflection coefficients of a plane
%     // dielectric surface
%     var result = new Array();
%     var one=Math.Complex(1.,0.);
    one = complex(1.,0.);
    n1 = sqrt(eps1);
    n2 = sqrt(eps2);

    tt = n1*sin(the1)./n2;  %tt = (n1.mul(Math.sin(the1))).div(n2);
%     //costh2 = Math.sqrt(1. - tt.mul(tt));
    costh2 = sqrt(one-tt.*tt);%(one.sub(tt.mul(tt))).sqrt();
    costh1 = cos(the1);%Math.cos(the1);

    rho_v = (n2*costh1 - n1*costh2)/(n2*costh1 + n1*costh2);%

    rho_v = rho_v*(-1.);
    rho_h = (n1*costh1 - n2*costh2)/(n1*costh1 + n2*costh2);

    result(1) = rho_v;
    result(2) = rho_h;

%     return result;
end
% 
% //************************************************************
% function refl_coef(the1, eps1, eps2) {
%     // calculates the v and h-polarized reflection coefficients of a plane
%     // dielectric surface
%     var result = new Array();
%     var one=Math.Complex(1.,0.);
% 
%     var n1 = eps1.sqrt();
%     var n2 = eps2.sqrt();
% 
%     var tt = (n1.mul(Math.sin(the1))).div(n2);
%     //costh2 = Math.sqrt(1. - tt.mul(tt));
%     var costh2 = (one.sub(tt.mul(tt))).sqrt();
%     var costh1 = Math.cos(the1);
% 
%     var rho_v =  ( (n2.mul(costh1)).sub(n1.mul(costh2))   ).div(   (n2.mul(costh1)).add(n1.mul(costh2)) );
%     rho_v = rho_v.mul(-1.);
%     var rho_h =  ( (n1.mul(costh1)).sub(n2.mul(costh2))    ).div(  (n1.mul(costh1)).add(n2.mul(costh2)) );
% 
%     result[0] = rho_v;
%     result[1] = rho_h;
% 
%     return result;
% }
% //************************************************************
function nvalue = Fresn_Refl0(eps)
%     // calculates Fresnel reflectivity at normal incidence.
%     var ss = Math.Complex(0.,0.);
%     var one = Math.Complex(1.,0.);
    one = complex(1.,0.);
    ss = sqrt(eps);
    nvalue = norm((one-ss)./(one+ss));
%     nvalue = sqrt(abs((one-ss)./(one+ss)));
%     return (one.sub(ss)).div(one.add(ss)).norm();
end

% //************************************************************
% function Fresn_Refl0(eps) {
%     // calculates Fresnel reflectivity at normal incidence.
%     var ss = Math.Complex(0.,0.);
%     var one = Math.Complex(1.,0.);
%     ss = eps.sqrt();
%     return (one.sub(ss)).div(one.add(ss)).norm();
% }
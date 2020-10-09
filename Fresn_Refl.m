
function result = Fresn_Refl(eps,theta) 
%     // calculates Fresnel reflectivities of v and h-polarizations at given set
%     // of incidence angles.
%     var result = new Array();
%     var rho = [];
%     var one=Math.Complex(1.,0.);
    one = complex(1.,0.);
    rho = refl_coef(theta, one, eps);

    result(1) = norm(rho(1));
    result(2) = norm(rho(2));

%     return result;
end

% //************************************************************
% function Fresn_Refl(eps,theta) {
%     // calculates Fresnel reflectivities of v and h-polarizations at given set
%     // of incidence angles.
%     var result = new Array();
%     var rho = [];
%     var one=Math.Complex(1.,0.);
% 
%     rho = refl_coef(theta, one, eps);
% 
%     result[0] = rho[0].norm();
%     result[1] = rho[1].norm();
% 
%     return result;
% }
% //************************************************************
function [t1, t2, gammah, gammav, t3, t4, t5, t6] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d) 
   
%     var result = new Array();

%     var theta1 = theta1d*Math.PI/180.;
    theta1 = theta1d*pi/180.;


%     var theta, h, v;
%     var eps1=[];
%     var eps2=[];
% 
%     var rhov=[];
%     var rhoh=[];
% 
%     var tauh=[];
%     var tauv=[];
% 
%     var sin_theta2=[];
%     var cos_theta2=[];
%     var one=[];
% 
    one(1)=1.0;
    one(2)=0.0;
    eps1r = real(eps1);
    eps1i = imag(eps1);
    eps2r = real(eps2);
    eps2i = imag(eps2);
    eps1(1) = eps1r;
    eps1(2) = eps1i;
    eps2(1) = eps2r;
    eps2(2) = eps2i;
    
%     var j;
 


    sin_theta2 = crmult(cdiv(csqrt(eps1), csqrt(eps2)) , sin(theta1));
    cos_theta2 = csqrt(csub(one,cmult( sin_theta2,sin_theta2)));

    rhoh = csub(crmult(csqrt(eps1),cos(theta1)), cmult(cos_theta2,csqrt(eps2)));
    rhoh = cdiv(rhoh , cadd(crmult(csqrt(eps1),cos(theta1)) , cmult(csqrt(eps2),cos_theta2)));

    rhov = csub(cmult(csqrt(eps1),cos_theta2), crmult(csqrt(eps2), cos(theta1)));
    rhov = cdiv(rhov, cadd(cmult(csqrt(eps1),cos_theta2) , crmult(csqrt(eps2), cos(theta1)) ));

    tauh = cadd(one,rhoh);
    tauv = crmult(cdiv(cadd(one , rhov),cos_theta2), cos(theta1));

    h = cabs(tauh);
    v = cabs(tauv);



    gammah = cabs_sqrd(rhoh);
    gammav = cabs_sqrd(rhov);
        
    Th = 1-gammah;
    Tv = 1-gammav;

   t1= rhoh;
   t2= rhov;
%    gammah = gammah;
%    gammav = gammav;
   t3 = tauh;
   t4 = tauv;
   t5 = Th;
   t6 = Tv;
            
%    return result;

end

% 
% //************************************************************
% function ReflTransm_PlanarBoundary(eps1r, eps1i, eps2r, eps2i, theta1d) {
%    
%     var result = new Array();
% 
%     var theta1 = theta1d*Math.PI/180.;
% 
% 
%     var theta, h, v;
%     var eps1=[];
%     var eps2=[];
% 
%     var rhov=[];
%     var rhoh=[];
% 
%     var tauh=[];
%     var tauv=[];
% 
%     var sin_theta2=[];
%     var cos_theta2=[];
%     var one=[];
% 
%     one[0]=1.0;
%     one[1]=0.0;
% 
%     eps1[0] = eps1r;
%     eps1[1] = eps1i;
%     eps2[0] = eps2r;
%     eps2[1] = eps2i;
%     
%     var j;
%  
% 
% 
%     sin_theta2 = crmult(cdiv(csqrt(eps1), csqrt(eps2)) , Math.sin(theta1));
%     cos_theta2 = csqrt(csub(one,cmult( sin_theta2,sin_theta2)));
% 
%     rhoh = csub(crmult(csqrt(eps1),Math.cos(theta1)), cmult(cos_theta2,csqrt(eps2)));
%     rhoh = cdiv(rhoh , cadd(crmult(csqrt(eps1),Math.cos(theta1)) , cmult(csqrt(eps2),cos_theta2)));
% 
%     rhov = csub(cmult(csqrt(eps1),cos_theta2), crmult(csqrt(eps2),Math.cos(theta1)));
%     rhov = cdiv(rhov, cadd(cmult(csqrt(eps1),cos_theta2) , crmult(csqrt(eps2),Math.cos(theta1)) ));
% 
%     tauh = cadd(one,rhoh);
%     tauv = crmult(cdiv(cadd(one , rhov),cos_theta2),Math.cos(theta1));
% 
%     h = cabs(tauh);
%     v = cabs(tauv);
% 
% 
% 
%     gammah = cabs_sqrd(rhoh);
%     gammav = cabs_sqrd(rhov);
%         
%     Th = 1-gammah;
%     Tv = 1-gammav;
% 
% 
%    result[0]= rhoh;
%    result[1]= rhov;
%    result[2]= gammah;
%    result[3]= gammav;
%    result[4]= tauh;
%    result[5]= tauv;
%    result[6]= Th;
%    result[7]= Tv;
%             
%    return result;
% 
% }

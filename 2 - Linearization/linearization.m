function [A,B1,B2,C,D1,D2,Ce,De1,De2]= linearization(AFRs,Ct,D,Hu,I,Pb_0,R,Tatm_0,Tman,V,Vd,a,alpha00,alpha_0,c,critpRatio,eta_b,eta_vn0,eta_vn1,eta_vn2,eta_vp1,k,mff_0,mfidot_0,n_0,patm_0,pman_0,sigma3,sigma4,sigma5,sigma6)
%% A_fcn
t2 = cos(alpha00);
t3 = cos(alpha_0);
t4 = critpRatio.*patm_0;
t5 = eta_vp1.*pman_0;
t6 = a.^2;
t7 = k+1.0;
t8 = n_0.^2;
t9 = 1.0./AFRs;
t10 = 1.0./I;
t11 = 1.0./R;
t12 = 1.0./Tman;
t13 = 1.0./V;
t14 = eta_vn1.*6.0e+1;
t15 = k-1.0;
t16 = 1.0./k;
t17 = n_0.*6.0e+1;
t18 = 1.0./n_0;
t20 = 1.0./patm_0;
t21 = sigma4-1.0;
t22 = 1.0./sigma3;
t27 = mfidot_0.*sigma5.*6.0e+1;
t34 = eta_vn2.*n_0.*7.2e+3;
t19 = 1.0./t8;
t23 = 1.0./t2;
t24 = 1.0./t3;
t25 = n_0.*t14;
t26 = -t4;
t28 = sigma5.*t17;
t29 = t16.*2.0;
t30 = pman_0.*t20;
t31 = t17.^sigma4;
t33 = 1.0./t15;
t35 = t7.*t16;
t36 = t17.^t21;
t38 = eta_vn2.*t8.*3.6e+3;
t42 = t14+t34;
t32 = pman_0+t26;
t37 = sigma6+t28-1.0;
t41 = mff_0.*t22.*t31;
t46 = mff_0.*sigma4.*t22.*t36.*6.0e+1;
t48 = eta_vn0+t5+t25+t38;
t39 = patm_0.*t32;
t40 = mfidot_0.*t37;
t45 = -t41;
t47 = -t46;
t49 = t27+t47;
t50 = t40+t45;
t52 = 1.0./t50;
t53 = t52.^2;
et1 = sqrt(2.0).*Ct.*D.^2.*R.*Tman.*t13.*k.*patm_0.*pi.*1.0./sqrt(R.*Tatm_0).*t33.*(t20.*t30.^(t29-1.0).*t29+t20.*t30.^(t35-1.0).*-1.0.*t35).*1.0./sqrt(k.*t33.*(t30.^t29-t30.^t35));
et2 = (asin(a).*6.366197723675814e-1+a.*sqrt(-t6+1.0).*6.366197723675814e-1+t23.*t3-t23.*t3.*asin(a.*t2.*t24).*6.366197723675814e-1-a.*t24.*sqrt(1.0./t24.^2-t6.*1.0./t23.^2).*6.366197723675814e-1-1.0).*(-1.0./8.0);
if ~all(cellfun(@isscalar,{t39}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (0.0 < t39)
    t0 = et1.*et2-(Vd.*n_0.*t5.*t13)./2.0-(Vd.*n_0.*t13.*t48)./2.0;
elseif (t39 < 0.0)
    t0 = Vd.*n_0.*t5.*t13.*(-1.0./2.0)-(Vd.*n_0.*t13.*t48)./2.0;
else
    t0 = NaN;
end
mt1 = [t0,0.0,0.0,Vd.*c.*n_0.*t5.*t9.*t11.*t12.*t52.*(-1.0./2.0)-(Vd.*c.*n_0.*t9.*t11.*t12.*t48.*t52)./2.0,Vd.*pman_0.*t13.*t48.*(-1.0./2.0)-(Vd.*n_0.*pman_0.*t13.*t42)./2.0];
mt2 = [Pb_0.*t10.*t19.*2.533029591058444e-2-Hu.*eta_b.*t10.*t18.*t49.*2.533029591058444e-2+Hu.*eta_b.*t10.*t19.*t50.*2.533029591058444e-2,t49,Vd.*c.*pman_0.*t9.*t11.*t12.*t48.*t52.*(-1.0./2.0)-(Vd.*c.*n_0.*pman_0.*t9.*t11.*t12.*t42.*t52)./2.0+(Vd.*c.*n_0.*pman_0.*t9.*t11.*t12.*t48.*t49.*t53)./2.0,0.0,Hu.*eta_b.*t10.*t18.*t22.*t31.*2.533029591058444e-2];
mt3 = [-t22.*t31,Vd.*c.*n_0.*pman_0.*t9.*t11.*t12.*t22.*t31.*t48.*t53.*(-1.0./2.0),0.0,0.0,0.0,-c];
A = reshape([mt1,mt2,mt3],4,4);

%% B1_fcn
t2 = cos(alpha00);
t3 = cos(alpha_0);
t4 = sin(alpha_0);
t5 = R.*Tatm_0;
t6 = critpRatio.*patm_0;
t7 = D.^2;
t8 = a.^2;
t9 = k+1.0;
t12 = 1.0./V;
t13 = k-1.0;
t14 = 1.0./k;
t15 = 1.0./patm_0;
t20 = n_0.*sigma5.*6.0e+1;
t21 = sqrt(2.0);
t10 = t2.^2;
t11 = t3.^2;
t16 = 1.0./t2;
t17 = 1.0./t3;
t19 = -t6;
t22 = 1.0./t9;
t23 = pman_0.*t15;
t25 = 1.0./t13;
t28 = 1.0./sqrt(t5);
t29 = sigma6+t20-1.0;
t18 = 1.0./t11;
t24 = pman_0+t19;
t26 = t8.*t10;
t27 = t4.*t16;
t30 = a.*t2.*t17;
t31 = patm_0.*t24;
t32 = asin(t30);
t33 = -t26;
t35 = t18.*t33;
t36 = t11+t33;
t41 = t27.*t32.*6.366197723675814e-1;
t37 = t35+1.0;
t38 = sqrt(t36);
t42 = -t41;
t39 = 1.0./t38;
t40 = 1.0./sqrt(t37);
t45 = a.*t4.*t18.*t38.*6.366197723675814e-1;
t43 = a.*t4.*t39.*6.366197723675814e-1;
t46 = a.*t4.*t17.*t40.*6.366197723675814e-1;
t44 = -t43;
if ~all(cellfun(@isscalar,{t31}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (0.0 <= t31)
    t0 = (Ct.*R.*Tman.*patm_0.*t7.*t12.*t21.*t28.*pi.*sqrt(-k.*t25.*(t23.^(t9.*t14)-t23.^(t14.*2.0))).*(t27+t42+t44+t45+t46))./4.0;
elseif (t31 < 0.0)
    t0 = (Ct.*R.*Tman.*patm_0.*t7.*t12.*t21.*t28.*pi.*sqrt(k.*t22).*(t22.*2.0).^t25.*(t27+t42+t44+t45+t46))./4.0;
else
    t0 = NaN;
end
B1 = reshape([0.0,(Hu.*eta_b.*t29.*(-2.533029591058444e-2))./(I.*n_0),sigma6+t20,(Vd.*c.*n_0.*pman_0.*t29.*1.0./(mfidot_0.*t29-(mff_0.*(n_0.*6.0e+1).^sigma4)./sigma3).^2.*(eta_vn0+eta_vn1.*n_0.*6.0e+1+eta_vp1.*pman_0+eta_vn2.*n_0.^2.*3.6e+3))./(AFRs.*R.*Tman.*2.0),t0,0.0,0.0,0.0],[4,2]);

%% B2_fcn
t2 = asin(a);
t3 = cos(alpha00);
t4 = cos(alpha_0);
t5 = R.*Tatm_0;
t6 = critpRatio.*patm_0;
t7 = D.^2;
t8 = R.^2;
t9 = a.^2;
t10 = k+1.0;
t13 = 1.0./V;
t14 = k-1.0;
t15 = 1.0./k;
t16 = 1.0./patm_0;
t21 = sqrt(2.0);
t11 = t3.^2;
t12 = t4.^2;
t17 = t16.^2;
t18 = 1.0./t3;
t19 = 1.0./t4;
t20 = -t6;
t22 = -t9;
t23 = t15.*2.0;
t24 = 1.0./t10;
t25 = pman_0.*t16;
t27 = 1.0./t14;
t30 = 1.0./sqrt(t5);
t35 = t10.*t15;
t50 = t2.*6.366197723675814e-1;
t26 = pman_0+t20;
t29 = t4.*t18;
t31 = t30.^3;
t32 = t22+1.0;
t33 = t24.*2.0;
t34 = k.*t24;
t36 = a.*t3.*t19;
t39 = t11.*t22;
t37 = patm_0.*t26;
t38 = asin(t36);
t41 = sqrt(t34);
t42 = sqrt(t32);
t44 = t12+t39;
t46 = t33.^t27;
t40 = (t37 < 0.0);
t48 = sqrt(t44);
t53 = a.*t42.*6.366197723675814e-1;
t54 = t29.*t38.*6.366197723675814e-1;
t55 = -t54;
t56 = a.*t19.*t48.*6.366197723675814e-1;
t57 = -t56;
if ~all(cellfun(@isscalar,{t37,t40}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (0.0 < t37)
    t0 = Ct.*R.*Tman.*t7.*t13.*t21.*t30.*pi.*sqrt(k.*t27.*(t25.^t23-t25.^t35)).*(t29+t50+t53+t55+t57-1.0).*(-1.0./4.0)+(Ct.*R.*Tman.*k.*patm_0.*t7.*t13.*t21.*t27.*t30.*pi.*(pman_0.*t17.*t23.*t25.^(t23-1.0)-pman_0.*t17.*t25.^(t35-1.0).*t35).*1.0./sqrt(k.*t27.*(t25.^t23-t25.^t35)).*(t29+t50+t53+t55+t57-1.0))./8.0;
elseif (t40)
    t0 = Ct.*R.*Tman.*t7.*t13.*t21.*t30.*t41.*t46.*pi.*(t29+t50+t53+t55+t57-1.0).*(-1.0./4.0);
else
    t0 = NaN;
end
mt1 = [0.0,(-2.533029591058444e-2)./(I.*n_0),0.0,0.0,t0,0.0,0.0,0.0];
if ~all(cellfun(@isscalar,{t37,t40}))
    error(message('symbolic:sym:matlabFunction:ConditionsMustBeScalar'));
end
if (0.0 <= t37)
    t0 = (Ct.*Tman.*patm_0.*t7.*t8.*t13.*t21.*t31.*pi.*sqrt(k.*t27.*(t25.^t23-t25.^t35)).*(t29+t50+t53+t55+t57-1.0))./8.0;
elseif (t40)
    t0 = (Ct.*Tman.*patm_0.*t7.*t8.*t13.*t21.*t31.*t41.*t46.*pi.*(t29+t50+t53+t55+t57-1.0))./8.0;
else
    t0 = NaN;
end
mt2 = [t0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
B2 = reshape([mt1,mt2],4,7);

%% C_fcn
C = reshape([0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0],[2,4]);

%% D1_fcn
D1 = reshape([0.0,0.0,0.0,0.0],[2,2]);

%% D2_fcn
D2 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],[2,7]);

%% Ce_fcn
Ce = reshape([0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0],[2,4]);

%% De1_fcn
De1 = reshape([0.0,0.0,0.0,0.0],[2,2]);

%% De2_fcn
De2 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,-1.0],[2,7]);

end
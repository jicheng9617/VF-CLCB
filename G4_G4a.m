function [obj,g] = G4_G4a( x,fidelity )

if ~isempty(x)
    x1=x(:,1);
    x2=x(:,2);
    x3=x(:,3);
    x4=x(:,4);
    x5=x(:,5);
    y=5.3578547*x3.^2+0.8356891*x1.*x5+37.293239*x1-40792.141;
    u1=85.334407+0.0056858*x2.*x5+0.0006262*x1.*x4-0.0022053*x3.*x5;
    u2=80.51249+0.0071317*x2.*x5+0.0029955*x1.*x2+0.002183*x3.^2;
    u3=9.300961+0.0047026*x3.*x5+0.0012547*x1.*x3+0.0019085*x3.*x4;
    g1=u1-92;
    g2=-u1;
    g3=u2-110;
    g4=90-u2;
    g5=u3-25;
    g6=20-u3;
else
    y=[];
    g1=[];
    g2=[];
    g3=[];
    g4=[];
    g5=[];
    g6=[];
end

if strcmp(fidelity,'high')
    obj = y;
    g = [g1,g2,g3,g4,g5,g6];
elseif strcmp(fidelity,'low')
%     if ~isempty(x)
%         x1=x(:,1)+0.1;
%         x2=x(:,2)-0.1;
%         x3=x(:,3)+0.1;
%         x4=x(:,4)-0.1;
%         x5=x(:,5)+0.1;
%         obj=5.3578547*x3.^2+0.8356891*x1.*x5+37.293239*x1-40792.141;
%         u1=85.334407+0.0056858*x2.*x5+0.0006262*x1.*x4-0.0022053*x3.*x5;
%         u2=80.51249+0.0071317*x2.*x5+0.0029955*x1.*x2+0.002183*x3.^2;
%         u3=9.300961+0.0047026*x3.*x5+0.0012547*x1.*x3+0.0019085*x3.*x4;
%         g1=u1-92;
%         g2=-u1;
%         g3=u2-110;
%         g4=90-u2;
%         g5=u3-25;
%         g6=20-u3;
%         g=[g1,g2,g3,g4,g5,g6];
%     else
%         obj=[];
%         g=[];
%     end
    obj = 0.9 * y + 0.5;
    gl1 = 0.9 * g1 - 0.05;
    gl2 = 0.9 * g2 - 0.05;
    gl3 = 0.9 * g3 - 0.05;
    gl4 = 0.9 * g4 - 0.05;
    gl5 = 0.9 * g5 - 0.05;
    gl6 = 0.9 * g6 - 0.05;
    g = [gl1, gl2, gl3, gl4, gl5, gl6];
else
    error('Wrong fidelity string.');
end

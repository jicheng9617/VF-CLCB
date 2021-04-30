function [best_x_hf,best_x_lf] = Criterion_VFCLCB(krigingobj_lf,krigingobj_dis,krigingcon_dis,krigingcon_lf,w,cr,num_vari,design_space,fmin)
optimizer = optimoptions('particleswarm','SwarmSize',100,'MaxIterations',100,'MaxStallIterations',100,'Display','off', 'UseVectorized', true);
num_con = length(krigingcon_dis);
best_x_hf=[]; best_x_lf=[];
c = 1.96; % the confidence level at the constraint sampling
%% generate update sample according to constraints
infill_gHF = @(x)Function_VF_GLCB(x,'high',krigingcon_dis,krigingcon_lf,cr); % infill criterion
infill_gLF = @(x)Function_VF_GLCB(x,'low',krigingcon_dis,krigingcon_lf,cr);
[gcandidate_lf, glcb_lf] = particleswarm(infill_gLF,num_vari,design_space(1,:),design_space(2,:),optimizer);
[gcandidate_hf, glcb_hf] = particleswarm(infill_gHF,num_vari,design_space(1,:),design_space(2,:),optimizer);
% determine whether to update this candidate 
if fmin == inf
    if glcb_hf <= glcb_lf 
        best_x_hf = [best_x_hf; gcandidate_hf];
    else 
        best_x_lf = [best_x_lf; gcandidate_lf];
    end
else
    if glcb_hf <= glcb_lf
        for i=1:num_con
            [pre_hf(:,i),mse_hf(:,i)] = HKpredictor(gcandidate_hf,krigingcon_dis{i},krigingcon_lf{i});
        end
        [ug_hf,index_hf]=max(pre_hf); sg_hf=sqrt(mse_hf(index_hf));
        if ug_hf.*(ug_hf+c*sg_hf)<0 | ug_hf.*(ug_hf-c*sg_hf)<0
            best_x_hf = [best_x_hf; gcandidate_hf];
        end
    else
        for i=1:num_con
            [pre_lf(:,i),mse_lf(:,i)] = HKpredictor(gcandidate_lf,krigingcon_dis{i},krigingcon_lf{i});
        end
        [ug_lf,index_lf]=max(pre_lf); sg_lf=sqrt(mse_lf(index_lf));
        if ug_lf.*(ug_lf+c*sg_lf)<0 | ug_lf.*(ug_lf-c*sg_lf)<0 
            best_x_lf = [best_x_lf; gcandidate_lf];
        end
    end
end
%% generate update sample with respect to the objective
infill_yHF = @(x)Function_VF_PLCB(x,'high',krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr);
infill_yLF = @(x)Function_VF_PLCB(x,'low',krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr);
[ycandidate_lf, plcb_lf]=particleswarm(infill_yLF,num_vari,design_space(1,:),design_space(2,:),optimizer);
[ycandidate_hf, plcb_hf]=particleswarm(infill_yHF,num_vari,design_space(1,:),design_space(2,:),optimizer);
if (plcb_hf<=plcb_lf || abs(plcb_hf - plcb_lf) < 0.0001)  && (isempty(best_x_hf) || pdist2(gcandidate_hf,ycandidate_hf)>0.001)
    best_x_hf=[best_x_hf;ycandidate_hf];
elseif plcb_hf>plcb_lf  && (isempty(best_x_lf) || pdist2(gcandidate_lf,ycandidate_lf)>0.0001)
    best_x_lf=[best_x_lf;ycandidate_lf];
end
end


function plcb_val=Function_VF_PLCB(x,fi,krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr)
w1=w(1); w2=w(2);
alpha=1e20;
num_con=length(krigingcon_dis);
nx=size(x,1);
gp_hf=zeros(nx,num_con);
gs2_hf=zeros(nx,num_con);
for j=1:num_con
    [gp_hf(:,j),gs2_hf(:,j)] = HKpredictor(x,krigingcon_dis{j},krigingcon_lf{j});
end
[yp_hf,ys2_hf,~,ys2_lf]=HKpredictor(x,krigingobj_dis,krigingobj_lf);
beta_yh = krigingobj_dis.beta;
gs_hf=sqrt(gs2_hf);
penalty=alpha*max(max(gp_hf,[],2),0);
ys_hf=sqrt(ys2_hf);
ys_lf=sqrt(ys2_lf) .* beta_yh;
if strcmp(fi,'high')
    lcb_val=w1*yp_hf-w2*ys_hf;
elseif strcmp(fi,'low')
    lcb_val=w1*yp_hf-w2*ys_lf*cr;
else
    error('Inaccurate fidelity string.');
end
plcb_val=lcb_val+penalty;
end


function glcb_val=Function_VF_GLCB(x,fi,krigingcon_dis,krigingcon_lf,cr)
a=1;
num_con=length(krigingcon_dis);
nx=size(x,1);
g_hf=zeros(nx,num_con);
gs2_hf=zeros(nx,num_con);
gs2_lf=zeros(nx,num_con);
for j=1:num_con
    [g_hf(:,j),gs2_hf(:,j),~,gs2_lf(:,j)] = HKpredictor(x,krigingcon_dis{j},krigingcon_lf{j});
end
gs_hf=sqrt(gs2_hf);
gs_lf=sqrt(gs2_lf);
glcb_val=zeros(nx,1);
for i=1:nx
    [~,index]=max(g_hf(i,:),[],2);
    if strcmp(fi,'high')
        glcb_val(i)=abs(g_hf(i,index))-a*gs_hf(i,index);
    elseif strcmp(fi,'low')
        glcb_val(i)=abs(g_hf(i,index))-a*gs_lf(i,index) * (cr);
    else
        error('Inaccurate fidelity string.');
    end
end
end

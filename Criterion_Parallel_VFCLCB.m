function [best_x_hf,best_x_lf] = Criterion_Parallel_VFCLCB(n_con,n_obj,krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr,num_vari,design_space)
options = optimoptions('particleswarm','SwarmSize',100,'MaxIterations',100,'MaxStallIterations',100,'Display','off', 'UseVectorized', true);
lb=design_space(1,:); ub=design_space(2,:);
best_x_hf=[];best_x_lf=[];
%% selection part for constraint conditions
for i=1:n_con
    [candidate_hf,glcb_hf]=particleswarm(@(x)Function_ParaVF_GLCB(x,'high',krigingcon_dis,krigingcon_lf,cr,best_x_hf,best_x_lf),num_vari,lb,ub,options);
    [candidate_lf,glcb_lf]=particleswarm(@(x)Function_ParaVF_GLCB(x,'low',krigingcon_dis,krigingcon_lf,cr,best_x_hf,best_x_lf),num_vari,lb,ub,options);
    if glcb_hf<=glcb_lf
        best_x_hf=[best_x_hf;candidate_hf];
    else
        best_x_lf=[best_x_lf;candidate_lf];
        for j=1:floor(cr-1)
            [candidate_lf,glcb_lf]=particleswarm(@(x)Function_ParaVF_GLCB(x,'low',krigingcon_dis,krigingcon_lf,cr,best_x_hf,best_x_lf),num_vari,lb,ub,options);
            best_x_lf=[best_x_lf;candidate_lf];
        end
    end
end
%% selection part for objective function
for i=1:n_obj
    [candidate_hf,plcb_hf]=particleswarm(@(x)Function_ParaVF_PLCB(x,'high',krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr,best_x_hf,best_x_lf),num_vari,lb,ub,options);
    [candidate_lf,plcb_lf]=particleswarm(@(x)Function_ParaVF_PLCB(x,'low',krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr,best_x_hf,best_x_lf),num_vari,lb,ub,options);
    if plcb_hf<=plcb_lf
        best_x_hf=[best_x_hf;candidate_hf];
    else
        best_x_lf=[best_x_lf;candidate_lf];
        for j=1:floor(cr-1)
            [candidate_lf,plcb_lf]=particleswarm(@(x)Function_ParaVF_PLCB(x,'low',krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr,best_x_hf,best_x_lf),num_vari,lb,ub,options);
            best_x_lf=[best_x_lf;candidate_lf];
        end
    end
end
end


function plcb_val=Function_ParaVF_PLCB(x,fi,krigingobj_dis,krigingobj_lf,krigingcon_dis,krigingcon_lf,w,cr,Xadd_hf,Xadd_lf)
w1=w(1);w2=w(2);
alpha=1e20;
num_con=length(krigingcon_dis);
nx=size(x,1);
gp_hf=zeros(nx,num_con);
for j=1:num_con
    gp_hf(:,j) = HKpredictorPara(x,krigingcon_dis{j},krigingcon_lf{j},Xadd_hf,Xadd_lf);
end
[yp_hf,ys2_hf,~,ys2_lf] = HKpredictorPara(x,krigingobj_dis,krigingobj_lf,Xadd_hf,Xadd_lf);
beta_yh = krigingobj_dis.beta;
penalty=alpha*max(max(gp_hf,[],2),0);
ys_hf=sqrt(ys2_hf);
ys_lf=sqrt(ys2_lf) .* beta_yh;
if strcmp(fi,'high')
    lcb_val = w1 * yp_hf - w2 * ys_hf;
elseif strcmp(fi,'low')
    lcb_val = w1 * yp_hf - w2 * ys_lf * log(cr);
else
    error('Inaccurate fidelity string.');
end
plcb_val=lcb_val+penalty;
end


function glcb_val=Function_ParaVF_GLCB(x,fi,krigingcon_dis,krigingcon_lf,cr,Xadd_hf,Xadd_lf)
a=1; % balance weight of boundary and uncertainty
num_con=length(krigingcon_dis);
nx=size(x,1);
g_hf=zeros(nx,num_con);
gs2_hf=zeros(nx,num_con);
gs2_lf=zeros(nx,num_con);
for j=1:num_con
    [g_hf(:,j),gs2_hf(:,j),~,gs2_lf(:,j)] = HKpredictorPara(x,krigingcon_dis{j},krigingcon_lf{j},Xadd_hf,Xadd_lf);
end
gs_hf=sqrt(gs2_hf);
gs_lf=sqrt(gs2_lf);
glcb_val=zeros(nx,1);
for i=1:nx
    [~,index]=max(g_hf(i,:),[],2);
    if strcmp(fi,'high')
        glcb_val(i) = abs(g_hf(i,index))-a * gs_hf(i,index);
    elseif strcmp(fi,'low')
        glcb_val(i) = abs(g_hf(i,index))-a * gs_lf(i,index) * log(cr);
    else
        error('Inaccurate fidelity string.');
    end
end
end
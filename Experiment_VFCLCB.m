clc; clearvars;close all;
addpath('dace');
%% setting problem and initializing parameters
fun_name = 'G4_G4a';
% specific information of the test problem
[num_vari,num_con,design_space,optimum,optimum_to_reach,hf_budget,cr] = TestFunctions_Multifidelity_Constrained(fun_name);
max_iteration = hf_budget * cr;
NInitHF = 3 * num_vari; % number of initial points for HF
NInitLF = 6 * num_vari;
sample_xh = lhsdesign(NInitHF, num_vari, 'criterion', 'maximin', 'iterations', 1000); % initial design points 
sample_xl = lhsdesign(NInitLF, num_vari, 'criterion', 'maximin', 'iterations', 1000);
[sample_yh,sample_gh] = feval(fun_name,sample_xh,'high'); % get the responses of initial points
[sample_yl,sample_gl] = feval(fun_name,sample_xl,'low');
iter = 0; 
HFE = 0; LFE = 0; NEFE = 0; % initialize cost records
f_min = zeros(max_iteration+1,1); % initialize optimal solution record
flag0 = 1; % initial flag params
index = sum(sample_gh <= 0, 2) == num_con; % evaluate optimal feasible solution
if sum(index) == 0
    f_min(iter+1)= inf;
    fprintf('iteration: %d, HFE: %d, LFE: %d, NEFE:%f, best solution: no feasiable solution.\n', iter, HFE, LFE, NEFE);
else
    f_min(iter+1)=min(sample_yh(index, :));
    fprintf('iteration: %d, HFE: %d, LFE: %d, NEFE:%f, best solution: %f.\n', iter, HFE, LFE, NEFE, f_min(iter+1));
end
%% Iterative procedure
while NEFE <= hf_budget && f_min(iter+1,1) >  optimum_to_reach 
    
    % build/rebuild kriging models
    krigingcon_lf=cell(1, num_con); krigingcon_hf=cell(1,num_con); 
    [krigingobj_lf, krigingobj_hf] = HKmodeling(sample_xh,sample_xl,sample_yh,sample_yl,1e0*randn(1,num_vari),1e-5*ones(1,num_vari),1e2*ones(1,num_vari));
    for ii = 1: num_con
        [krigingcon_lf{ii}, krigingcon_hf{ii}] = HKmodeling(sample_xh,sample_xl,sample_gh(:,ii),sample_gl(:,ii),1e0*randn(1,num_vari),1e-5*ones(1,num_vari),1e2*ones(1,num_vari));
    end
    w = [1, 2 + log(flag0+1)]; % weighting params in LCB function
    [best_xh,best_xl] = Criterion_VFCLCB(krigingobj_lf,krigingobj_hf,krigingcon_hf,krigingcon_lf,w,cr,num_vari,design_space,f_min(iter+1)); % infill criterion
    
    % evaluate the candidates with the real function
    best_yh=[]; best_gh=[]; best_yl=[]; best_gl=[];
    if ~isempty(best_xh)
        [best_yh,best_gh]=feval(fun_name,best_xh,'high');
        sample_xh=[sample_xh;best_xh];
        sample_yh=[sample_yh;best_yh];
        sample_gh=[sample_gh;best_gh];
    end
    if ~isempty(best_xl)
        [best_yl,best_gl]=feval(fun_name,best_xl,'low');
        sample_xl=[sample_xl;best_xl];
        sample_yl=[sample_yl;best_yl];
        sample_gl=[sample_gl;best_gl];
    end
    
    % update parameters
    iter = iter + 1;
    HFE = size(sample_xh,1) - NInitHF;
    LFE = size(sample_xl,1) - NInitLF;
    NEFE = HFE + LFE / cr;
    index = sum(sample_gh <= 0, 2) == num_con;
    if sum(index) == 0
        f_min(iter+1)= inf;
        fprintf('iteration: %d, HFE: %d, LFE: %d, NEFE:%f, best solution: no feasiable solution.\n', iter, HFE, LFE, NEFE);
    else
        f_min(iter+1)=min(sample_yh(index, :));
        fprintf('iteration: %d, HFE: %d, LFE: %d, NEFE:%f, best solution: %f.\n', iter, HFE, LFE, NEFE, f_min(iter+1));
    end
    relative_error = abs((f_min(iter+1,1)-optimum)/optimum)*100;
    absolute_error = abs(f_min(iter+1,1)-optimum);
    if f_min(iter+1) == f_min(iter)
        flag0 = flag0 + 1;
    else 
        flag0 = 1;
    end

end
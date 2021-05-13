function [q_con,q_obj] = Calculation_NumPara(q_con,q_obj,best_xh,best_xl,best_gh,best_gl,krigingcon_dis,krigingcon_lf,beta,f_min,iteration)
num_con=length(krigingcon_dis);
n_q=q_con+q_obj;
nh=size(best_xh,1);
nl=size(best_xl,1);

conInd_hf = 0; conInd_lf = 0;
if nh > 0
    pg_xh = zeros(nh,num_con);
    for i=1:num_con
        pg_xh(:,i)=HKpredictor(best_xh,krigingcon_dis{i},krigingcon_lf{i});
    end
    conInd_hf=sum(best_gh.*pg_xh<=0,2)>0;
end
if nl > 0
    pg_xl = zeros(nl,num_con);
    for i=1:num_con
        [~,~,pg_xl(:,i)]=HKpredictor(best_xl,krigingcon_dis{i},krigingcon_lf{i});
    end
    conInd_lf=sum(best_gl.*pg_xl<=0,2)>0;
end 
conInd=sum(conInd_hf)+sum(conInd_lf);

% if conInd >= beta * (nh + nl) && (f_min(iteration+1) == f_min(iteration))
if conInd >= beta * (nh + nl)
    q_con = min(n_q,q_con+1);
    q_obj = n_q-q_con;
elseif f_min(iteration+1) == f_min(iteration)
    q_obj = min(n_q,q_obj+1);
    q_con = n_q-q_obj;
end

if q_obj == 0 
    q_obj = q_obj + 1;
    q_con = q_con - 1; 
end
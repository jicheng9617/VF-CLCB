function [num_vari,num_con,design_space,optimum,optimum_to_reach,hf_budget,cr,max_iter] = TestFunctions_Multifidelity_Constrained( fun_name )
hf_budget=[];
cr = 4;
switch fun_name
    
    % numerical cases 
    case 'Gano2_Gano2a'; num_vari=2;num_con=1;design_space=[0.1,0.1;2,2];optimum=5.668365;optimum_to_reach=5.670; hf_budget=50;max_iter = 20;
    case 'G5MOD_bifidelity'; num_vari=4; num_con=5; design_space=[0,0,-0.55,-0.55;1200,1200,0.55,0.55]; optimum=5126.50; optimum_to_reach=5130; hf_budget=50;max_iter = 20;
    case 'G4_G4a'; num_vari=5; num_con=6;design_space=[78,33,27,27,27;102,45,45,45,45];optimum=-30665.539;optimum_to_reach=-30665.000; hf_budget=50;max_iter = 20;
    case 'G24_bifidelity'; num_vari=2; num_con=2; design_space=[0,0;3,4]; optimum=-5.5080; optimum_to_reach=-5.5070; hf_budget=50;max_iter = 20;
    case 'Hesse_bifidelity'; num_vari=6; num_con=6; design_space=[0,0,1,0,1,0;5,4,5,6,5,10]; optimum=-310.00; optimum_to_reach=-309; hf_budget=50;max_iter = 20;
    case 'Reducer_bifidelity'; num_vari=7; num_con=11; design_space=[2.6,0.7,17,7.3,7.3,2.9,5;3.6,0.8,28,8.3,8.3,3.9,5.5]; optimum=2994.42; optimum_to_reach=2995; hf_budget=50;max_iter = 20;
    case 'Pressure_Vessel_bifidelity'; num_vari=4; num_con=3; design_space=[1,0.625,25,25;1.375,1,150,240]; optimum=7.0068e+03; optimum_to_reach=7007; hf_budget=50;max_iter = 20;
    case 'G1'; num_vari=13; design_space=[zeros(1,13);1,1,1,1,1,1,1,1,1,5,5,5,1]; optimum=-15; num_con=9; optimum_to_reach=-14.8; hf_budget=100; max_iter = 30;
    case 'G6'; num_vari=2; design_space=[13,0;100,100]; optimum=-6961.81388; num_con=2; optimum_to_reach=-6960; hf_budget=50; max_iter = 20;
    case 'G7'; num_vari=10; design_space=repmat([-10;10],1,10); optimum=24.3062; num_con=8; optimum_to_reach=28; hf_budget=100; max_iter = 30;
    case 'G8'; num_vari=2; design_space=[0.0001,0.0001;10,10]; optimum=-0.095825; num_con=2; optimum_to_reach=-0.0957; hf_budget=100; max_iter = 30;
    case 'G9'; num_vari=7; design_space=repmat([-10;10],1,num_vari); optimum=680.6301; num_con=4; optimum_to_reach=1000; hf_budget=200; max_iter = 50;
        
    % engineering case
    case 'stiffened_cylindrical_shell'; num_vari=9; num_con=9; design_space=[18 14 80 12 250 22 130 20 450;28 26 120 22 300 35 170 30 500]; optimum = 0; optimum_to_reach = -1e10; hf_budget = 50; cr = 3; max_iter = 20;
        
end

end


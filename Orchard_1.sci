// Initialization
  
  clear
  mode(0);
  ieee(1);
  clearglobal;
  
//    T_wd = "E:\Laurent-AA02870\boulot\prod_sci\En_cours\Article_BEST_2K\scilab\7_BEST-DP_versus_SP";
    T_wd = "C:\Users\Utilisateur\Documents\travail_maison\z_September_2018\Article_BEST_2K\scilab\7_BEST-DP_versus_SP";
    T_functions = T_wd + "\BEST_functions-2";
    T_data = T_wd + "\Data_files";  
    
    i_fig = 1;
    
// Sample and conditions

//    names = ["Forest_1";"Forest_2";"Pasture_1";"Pasture_2";"Orchard_1";"Orchard_2"];

   i0 = 5;
   
   file_name = "O_run_5";
   
    gammaa = 0.75;
    betaa = 0.6;
        
/// functions for modelling hydraulic HCWRF with vGBBC model
    
    chdir(T_functions);
    
    exec('0_general_tool_functions.sci',-1);
    exec('2_BEST_functions.sci',-1);
    exec('3_It_functions_2K.sci',-1);
    exec('4_HCWRF_comp.sci',-1);  
    

// nomemclature & cleanup

   names = ["Forest_1";"Forest_2";"Pasture_1";"Pasture_2";"Orchard_1";"Orchard_2"];
   
   ntot = size(names,1);
   
   name = names(i0);
   
   close_windows(20);
   
// data collection for cumulative water infiltration
   
   // Upload of hydraulic properties
   
   chdir(T_data);
   
   best_dp= read('BEST-DP.txt',-1,ntot);        // units = mm & mm/s
   
   ph_m = best_dp(1:6,i0);
   wf = best_dp(7,i0);
   ph_f = best_dp(8:13,i0);
    
   best_sp = read('BEST-SP.txt',-1,ntot);        // units = mm & mm/s
   
   ph_sp = best_sp(1:6,i0);
   
   // Upload of cumulative infiltration
   
   CI30_exp_table = read("I30_exp.txt",-1,2*ntot);                // units = mm & mm/s
   CI0_exp_table = read("I0_exp.txt",-1,2*ntot);
   
   a_30 = CI30_exp_table(2:size(CI30_exp_table,1),(2*i0-1));
   [a b] = min(a_30);
   t30_exp = CI30_exp_table(1:b,(2*i0-1));
   I30_exp = CI30_exp_table(1:b,2*i0);
    

   a_0 = CI0_exp_table(2:size(CI0_exp_table,1),(2*i0-1));        // units = mm & mm/s
   [a b] = min(a_0);
   t0_exp = CI0_exp_table(1:b,(2*i0-1));
   I0_exp = CI0_exp_table(1:b,2*i0);
   
   
   // BEST entry
    
    chdir(T_data+"\BEST_data");
    
    best_entry = read(file_name+'.txt',-1,2);           //--> data file
    nbl = size(best_entry,1);                           // number of line of data file

    theta0s = best_entry(2,1);                          // initial water content  for Beerkan run [%]
    theta030 = best_entry(2,2);                         // initial water content  for Tension infil. [%]
    theta30 = best_entry(3,1);                          // initial water content [%]
    rds = best_entry(4,1);                              // ring radius for Beekan exp [mm]
    rd30 = best_entry(4,2);                             // ring radius for TI 30 [mm]
    
// Computation of hydraulic conductivity and water retention curves

    // water retention and hydraulic conductivity curves
    
    xh = (-logspace(-1,4,200))';
    
    theta_m = theta_vGB(xh,ph_m(3),ph_m(1),ph_m(2),ph_m(5));
    theta_f = theta_vGB(xh,ph_f(3),ph_f(1),ph_f(2),ph_f(5));
    theta_dp = (1-wf)*theta_m+wf*theta_f;
    theta_mat = (1-wf)*theta_m;
    
    theta_sp = theta_vGB(xh,ph_sp(3),ph_sp(1),ph_sp(2),ph_sp(5));

    K_m = K_vGB(ph_m(4),theta_m,ph_m(1),ph_m(2),ph_m(6));
    K_f = K_vGB(ph_f(4),theta_f,ph_f(1),ph_f(2),ph_f(6));
    K_dp = (1-wf)*K_m+wf*K_f;
    K_mat = (1-wf)*K_m;
    
    K_sp = K_vGB(ph_sp(4),theta_sp,ph_sp(1),ph_sp(2),ph_sp(6));
    
    // determination of initial water pressure heads

    h0_sp = wcr_vGB(theta030,ph_sp(1),ph_sp(2),ph_sp(3),ph_sp(5));
    
    [a b] = min(abs(theta_dp-theta030));
    h0_dp = xh(b);
    
// computation of original and final states

    // initial water contents and for B

    theta0_TI_sp = theta_vGB(h0_sp,ph_sp(3),ph_sp(1),ph_sp(2),ph_sp(5));
    
    theta0_TI_m = theta_vGB(h0_dp,ph_m(3),ph_m(1),ph_m(2),ph_m(5));
    theta0_TI_f = theta_vGB(h0_dp,ph_f(3),ph_f(1),ph_f(2),ph_f(5));
    theta0_TI_dp = wf*theta0_TI_f+(1-wf)*theta0_TI_m;               // prediction of water content at -30
    
    R_theta0_TI = theta0_TI_dp/theta030;                            // precision of the estimate
    
    theta0_B_sp = theta0_TI_sp;
    theta0_B_m = theta0_TI_sp;
    theta0_B_f = theta0_TI_sp;
    
    // final water contents for TI 
    
    thetaf_TI_sp = theta_vGB(-30,ph_sp(3),ph_sp(1),ph_sp(2),ph_sp(5));
    
    thetaf_TI_m = theta_vGB(-30,ph_m(3),ph_m(1),ph_m(2),ph_m(5));
    thetaf_TI_f = theta_vGB(-30,ph_f(3),ph_f(1),ph_f(2),ph_f(5));
    thetaf_TI_dp = wf*thetaf_TI_f+(1-wf)*thetaf_TI_m;               // prediction of water content at -30
    
    R_thetaf_TI = thetaf_TI_dp/theta30;                            // precision of the estimate
    
    // final water contents for B
    
    thetaf_B_sp = theta_vGB(0,ph_sp(3),ph_sp(1),ph_sp(2),ph_sp(5));
    
    thetaf_B_m = theta_vGB(0,ph_m(3),ph_m(1),ph_m(2),ph_m(5));
    thetaf_B_f = theta_vGB(0,ph_f(3),ph_f(1),ph_f(2),ph_f(5));
    thetaf_B_dp = wf*thetaf_B_f+(1-wf)*thetaf_B_m;       // prediction of water content at -30
    
    R_thetaf_B = thetaf_B_dp/thetaf_B_sp;                // precision of the estimate thetas_B_sp référence car défini par 1-rho_d/rho_s
    
    // initial hydraulic conductivities for boths TI and B, similar initial water contents)

    K0_TI_sp = K_vGB(ph_sp(4),theta0_TI_sp,ph_sp(1),ph_sp(2),ph_sp(6));
    
    K0_TI_m = K_vGB(ph_m(4),theta0_TI_m,ph_m(1),ph_m(2),ph_m(6));
    K0_TI_f = K_vGB(ph_f(4),theta0_TI_f,ph_f(1),ph_f(2),ph_f(6));
    
    K0_B_sp = K0_TI_sp;
    K0_B_m = K0_TI_m;
    K0_B_f = K0_TI_f;
    
    // final hydraulic conductivities for TI

    Kf_TI_sp = K_vGB(ph_sp(4),thetaf_TI_sp,ph_sp(1),ph_sp(2),ph_sp(6));
    
    Kf_TI_m = K_vGB(ph_m(4),thetaf_TI_m,ph_m(1),ph_m(2),ph_m(6));
    Kf_TI_f = K_vGB(ph_f(4),thetaf_TI_f,ph_f(1),ph_f(2),ph_f(6));
    
    Kf_TI_dp = (1-wf)*Kf_TI_m+wf*Kf_TI_f;
    
    // final hydraulic conductivities for B

    Kf_B_sp = K_vGB(ph_sp(4),thetaf_B_sp,ph_sp(1),ph_sp(2),ph_sp(6));
    
    Kf_B_m = K_vGB(ph_m(4),thetaf_B_m,ph_m(1),ph_m(2),ph_m(6));
    Kf_B_f = K_vGB(ph_f(4),thetaf_B_f,ph_f(1),ph_f(2),ph_f(6));
    
    // sorptivities

    S_TI_sp = (S_2_vGB(h0_sp,-30,ph_sp(1),ph_sp(2),ph_sp(3),ph_sp(4),ph_sp(5),ph_sp(6)))^(1/2);
    S_B_sp = (S_2_vGB(h0_sp,0,ph_sp(1),ph_sp(2),ph_sp(3),ph_sp(4),ph_sp(5),ph_sp(6)))^(1/2);
    
    S_TI_m = (S_2_vGB(h0_dp,-30,ph_m(1),ph_m(2),ph_m(3),ph_m(4),ph_m(5),ph_m(6)))^(1/2);
    S_TI_f = (S_2_vGB(h0_dp,-30,ph_f(1),ph_f(2),ph_f(3),ph_f(4),ph_f(5),ph_f(6)))^(1/2);
    
    S_B_m = (S_2_vGB(h0_dp,0,ph_m(1),ph_m(2),ph_m(3),ph_m(4),ph_m(5),ph_m(6)))^(1/2);
    S_B_f = (S_2_vGB(h0_dp,0,ph_f(1),ph_f(2),ph_f(3),ph_f(4),ph_f(5),ph_f(6)))^(1/2);
    
    
//  Computation of cumulative infiltrations

    // Cumulative infiltration into the SP soils
    
    I_TI_sp = I_3D(t30_exp,Kf_TI_sp,K0_TI_sp,S_TI_sp,rd30,thetaf_TI_sp,theta0_TI_sp,betaa,gammaa);
    I_B_sp = I_3D(t0_exp,Kf_B_sp,K0_B_sp,S_B_sp,rds,thetaf_B_sp,theta0_B_sp,betaa,gammaa);
    
    // Cumulative infiltration into the SP soils
    
    I_TI_m = I_3D(t30_exp,Kf_TI_m,K0_TI_m,S_TI_m,rd30,thetaf_TI_m,theta0_TI_m,betaa,gammaa);
    I_TI_f = I_3D(t30_exp,Kf_TI_f,K0_TI_f,S_TI_f,rd30,thetaf_TI_f,theta0_TI_f,betaa,gammaa);
    I_TI_dp = wf*I_TI_f+(1-wf)*I_TI_m;
    
    I_B_m = I_3D(t0_exp,Kf_B_m,K0_B_m,S_B_m,rds,thetaf_B_m,theta0_B_m,betaa,gammaa);
    I_B_f = I_3D(t0_exp,Kf_B_f,K0_B_f,S_B_f,rds,thetaf_B_f,theta0_B_f,betaa,gammaa);
    I_B_dp = wf*I_B_f+(1-wf)*I_B_m;
    
//////////////////////// Choix des couleurs et eppaisseur //////////////
//Rq : la definition des colormap doit etre avant les fonctions scf et clf !

     markers = ['+';'+';'<';'<';'o';'o'];
     color_ = ["green";"green";"blue";"blue";"red";"red"];
     
//////// Figure 1 //////////////////////////////////////////////////

     font_title = 3;
     font_xylabel = 3;
     font_axs = 2;
     font_leg = 2;
     
     clf(i_fig);
     scf(i_fig);
     
     i_fig = i_fig+1;
     
     subplot(221)
     
      plot(abs(xh),theta_sp,"color","blue");
      plot(abs(xh),theta_dp,"color","red");
      plot(30,thetaf_TI_sp,"bo",30,thetaf_TI_dp,"ro");
          
      title("BEST-DP "+name,"fontsize",font_title, "color", "black");
      xlabel("$h\ [mm]$","fontsize",font_xylabel, "color", "black");
      ylabel("$\theta\ [-]$","fontsize",font_xylabel, "color", "black");
      
      g = gca();
//      g.data_bounds = [1,0;10^4,0.8];
      g.log_flags = "lnn";
      g.font_size = font_axs;

     subplot(222)
     
      plot(abs(xh),K_sp,"color","blue");
      plot(abs(xh),K_dp,"color","red");
      plot(30,Kf_TI_sp,"bo",30,Kf_TI_dp,"ro");
      
      title("BEST-DP","fontsize",font_title, "color", "black");
      xlabel("$h\ [mm]$","fontsize",font_xylabel, "color", "black");
      ylabel("$K\ [mm/min]$","fontsize",font_xylabel, "color", "black");
      
      g = gca();
      g.log_flags = "lln";
//      g.data_bounds = [1,10^-6;10^4,10];
      g.font_size = font_axs;
      
      legend("BEST-SP","BEST-DP",3);
      
     subplot(223)
     
      plot(t30_exp,I30_exp,"color","black",'LineSt',"none","marker",markers(i0),"marksize",8);
      plot(t30_exp,I_TI_sp,"color","blue");
      plot(t30_exp,I_TI_dp,"color","red");

      xlabel("$t_{TI}\ [min]$","fontsize",font_xylabel);
      ylabel("$I_{TI}\ [mm]$","fontsize",font_xylabel);
      title('Cumulative infiltrations I30',"fontsize",font_title);
      
      g = gca();
      g.log_flags = 'nnn';
//      g.y_location = 'right';
//      g.data_bounds = [0,0;125,200];

     subplot(224)
     
      plot(t0_exp,I0_exp,"color","black",'LineSt',"none","marker",markers(i0),"marksize",8);
      plot(t0_exp,I_B_sp,"color","blue");
      plot(t0_exp,I_B_dp,"color","red");

      xlabel("$t_B\ [min]$","fontsize",font_xylabel);
      ylabel("$I_{B}\ [mm]$","fontsize",font_xylabel);
      title('Cumulative infiltrations I0',"fontsize",font_title);
      
      g = gca();
      g.log_flags = 'nnn';
//      g.y_location = 'right';
//      g.data_bounds = [0,0;125,200];
      legend("exp. data","BEST-SP","BEST-DP",2);


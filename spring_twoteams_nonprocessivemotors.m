clear; clc; % nonprocessive kinesin - 1
% run with cal_v_lim_Bi.m
% sometime if one motor bound and perform a step before new binding  then
% the simulation can stop
tic
threshold = 10^-5;  % consider 10^-4 as absolute 0 for V of cargo;
nrun=100;


Num_lane=[1,20,40,60,80,100];  % no. of left lane
[w,e]=size(Num_lane);
Nr =10; % no. of right lane
Nl =10; % no. of left lane
M=20;
nstep =10^6;
k = 0.52; % spring constant nN/microm
step_size = 8; %nm
%stepping rates
p=100;    % forward rate
q=10;  % backward rate
delta=0.5;

%binding rates
kon=5;  %constant binding rate 
koff=1;  % constant unbinding rate


for d= 1:e

    Nlane_l = Num_lane(d);
    Nlane_r = Num_lane(d);
    M_r=20*Nlane_r;   % the total number of binding sites on each side
    M_l=20*Nlane_l; 
    Vr0= cal_v_lim_Bi(Nr,M,0,kon,koff,p,q,delta)*step_size; % velocity on each lane
    Vl0= cal_v_lim_Bi(Nl,M,0,kon,koff,p,q,delta)*step_size;
    %Vr0 = cal_v_Bi(Nr,0,kon,koff,p,q,delta)*step_size;
    %Vl0 = cal_v_Bi(Nl,0,kon,koff,p,q,delta)*step_size;
   % Vr0=cal_v_pro(Nr,0,p,q,delta)*step_size;
   % Vl0=cal_v_pro(Nl,0,p,q,delta)*step_size;


    for round =1:nrun

        t=0;
        nr=zeros(Nlane_r,1);
        nl=zeros(Nlane_l,1);
        V_r=zeros(Nlane_r,1);
        V_l=zeros(Nlane_l,1);
        rec_Vr=[]; rec_Vl=[];
        rec_nr=[]; rec_nl=[];
        rec_aa=[]; rec_dx_move=[]; T=[];
        rec_dx_stretch=[]; rec_fT=[]; rec_fm=[]; rec_state=[];
        rec_xr=[]; rec_xl=[]; rec_bound_r=[]; rec_bound_l=[];
        x_bound_r=[]; % position of only bound motors on the right
        x_bound_l=[]; % position of only bound motors on the left
       % M_site=[];

        dxr=0;
        dxl=0;


        M=20; %cargo width and the number of binding sites
        x_r = zeros(Nlane_r,Nr); % position of all motors
        x_l = zeros(Nlane_l,Nl); % position of all motors
        bound_r = zeros(Nlane_r,Nr);% binding (1) and unbinding (0)
        bound_l = zeros(Nlane_l,Nl);% binding (1) and unbinding (0)
        rec_dxl = [];
        rec_dxr = [];
        rec_M_r = [];
        rec_M_l = [];
        state=0;
        tau=0;
        center_r=1000;
        center_l=-1000;
        

        for j=1:nstep

            if ~isempty(x_bound_r)
                dxr=(max(x_bound_r)-xr0)*step_size; %nm
            else
                dxr=0;
            end

            if ~isempty(x_bound_l)
                dxl=(xl0-min(x_bound_l))*step_size; %nm
            else
                dxl=0;
            end
            
            f_spring = k*(dxr+dxl); % pN
            f_motor = f_spring*1.87; %dimensionless unit


            for ir = 1:Nlane_r
                if nr(ir) > 0
                    fr_lane = f_motor/Nlane_r;
                    p1_r=p*exp(-fr_lane*delta);
                    q1_r=q*exp(fr_lane*(1-delta));
                    V_r(ir) = cal_v_lim_Bi(nr(ir),M_r_lane(ir),fr_lane,kon,koff,p,q,delta)*step_size;
                    %V_r(ir) = cal_v_Bi(nr(ir),fr_lane,kon,koff,p,q,delta)*step_size;
                   % V_r(ir) = cal_v_pro(nr(ir),fr_lane,p,q,delta)*step_size;
                else
                    fr_lane = f_motor/Nlane_r;
                    p1_r=p*exp(-fr_lane*delta);
                    q1_r=q*exp(fr_lane*(1-delta));
                    V_r(ir) = Vr0;
                end
            end


            for il =1:Nlane_l
                if nl(il) > 0
                    fl_lane = f_motor/Nlane_l;
                    p1_l=p*exp(-fl_lane*delta);
                    q1_l=q*exp(fl_lane*(1-delta));
                    V_l(il) = cal_v_lim_Bi(nl(il),M_l_lane(il),fl_lane,kon,koff,p,q,delta)*step_size;
                    %V_l(il) = cal_v_Bi(nl(il),fl_lane,kon,koff,p,q,delta)*step_size;
                    %V_l(il) = cal_v_pro(nl(il),fl_lane,p,q,delta)*step_size;
                else
                    fl_lane = f_motor/Nlane_l;
                    p1_l=p*exp(-fl_lane*delta);
                    q1_l=q*exp(fl_lane*(1-delta));
                    V_l(il) = Vl0;
                end
            end

            frac_Vr = V_r./Vr0;
            frac_Vl = V_l./Vl0;
            

            if any(frac_Vl < threshold) && any(frac_Vr < threshold) 
                fs_system = f_motor; %dimensionless unit
            end

              
                
            if all(frac_Vr > threshold) || all(frac_Vl > threshold)


                rec_Vl(:,j) = frac_Vl;
                rec_Vr(:,j) = frac_Vr;
                
            
               
                %check possibility of stepping
                %right team
                if all(x_r==0)
                    xr_all = nonzeros(x_r);
                    xr_lead_forward = nonzeros(x_r(:,1))+1;
                    xr_lead_backward = nonzeros(x_r(:,1))-1;
                    xr_fol_forward = nonzeros(x_r(:,2:Nr))+1;
                    xr_fol_backward = nonzeros(x_r(:,2:Nr))-1;
                else
                    new_xr= x_r;
                    new_xr(new_xr== 0) = NaN;
                    xr_all = nonzeros(x_r);
                    xr_lead = max(new_xr, [], 2);
                    xr_lead=(xr_lead(~isnan(xr_lead)));
                    xr_lead_forward = xr_lead+1;
                    xr_lead_backward = xr_lead-1;
                    xr_fol = setdiff(xr_all,xr_lead,'stable');
                    xr_fol_forward = xr_fol+1;
                    xr_fol_backward = xr_fol-1;
                end

                if any(frac_Vr < threshold)  
                xr_lead_forward =  xr_all;  
                xr_lead_backward =  xr_all;  
                xr_fol_forward = xr_all;
                xr_fol_backward =xr_all;
                end
                
               
                [xr_lead_p,kr_lead_fore] = setdiff(xr_lead_forward,xr_all,'stable');
                [xr_lead_b,kr_lead_back] = setdiff(xr_lead_backward,xr_all,'stable');
                [xr_fol_p,kr_fol_fore] = setdiff(xr_fol_forward,xr_all,'stable');
                [xr_fol_b,kr_fol_back] = setdiff(xr_fol_backward,xr_all,'stable');

                %left team
                if all(x_l==0)
                    xl_all = nonzeros(x_l);
                    xl_lead_forward = nonzeros(x_l(:,1))-1;
                    xl_lead_backward = nonzeros(x_l(:,1))+1;
                    xl_fol_forward = nonzeros(x_l(:,2:Nl))-1;
                    xl_fol_backward = nonzeros(x_l(:,2:Nl))+1;
                else
                    new_xl= x_l;
                    new_xl(new_xl== 0) = NaN;
                    xl_all = nonzeros(x_l);
                    xl_lead = min(new_xl, [], 2);
                    xl_lead=(xl_lead(~isnan(xl_lead)));
                    xl_lead_forward = xl_lead-1;
                    xl_lead_backward = xl_lead+1;
                    xl_fol = setdiff(xl_all,xl_lead,'stable');
                    xl_fol_forward = xl_fol-1;
                    xl_fol_backward = xl_fol+1;
                end

                if any(frac_Vl < threshold)
                xl_lead_forward =  xl_all;  
                xl_lead_backward =  xl_all;  
                xl_fol_forward = xl_all;
                xl_fol_backward =xl_all;
                end


                [xl_lead_p,kl_lead_fore] = setdiff(xl_lead_forward,xl_all,'stable');
                [xl_lead_b,kl_lead_back] = setdiff(xl_lead_backward,xl_all,'stable');
                [xl_fol_p,kl_fol_fore] = setdiff(xl_fol_forward,xl_all,'stable');
                [xl_fol_b,kl_fol_back] = setdiff(xl_fol_backward,xl_all,'stable');

                if isempty(kr_lead_fore)
                    aa(1)=0;
                else
                    aa(1)=size(kr_lead_fore,1);
                end

                if isempty(kr_lead_back)
                    aa(2)=0;
                else
                    aa(2)=size(kr_lead_back,1);
                end

                if isempty(kr_fol_fore)
                    aa(3)=0;
                else
                    aa(3)=size(kr_fol_fore,1);
                end

                if isempty(kr_fol_back)
                    aa(4)=0;
                else
                    aa(4)=size(kr_fol_back,1);
                end

                if isempty(kl_lead_fore)
                    aa(5)=0;
                else
                    aa(5)=size(kl_lead_fore,1);
                end

                if isempty(kl_lead_back)
                    aa(6)=0;
                else
                    aa(6)=size(kl_lead_back,1);
                end

                if isempty(kl_fol_fore)
                    aa(7)=0;
                else
                    aa(7)=size(kl_fol_fore,1);
                end

                if isempty(kl_fol_back)
                    aa(8)=0;
                else
                    aa(8)=size(kl_fol_back,1);
                end

                rec_aa(:,j) = aa;

                %Right team
                %Binding
                
                alpha(1,1) = (Nr*Nlane_r-sum(nr))*((M_r-sum(nr))/(M_r))*kon;  % N-n = number of unbound motors 
                alpha(1,2) = sum(nr)*koff;   % n = number of bound motors
                alpha(1,3) = aa(1)*p1_r;    % leading motor's forward rate
                alpha(1,4) = aa(2)*q1_r;    % leading motor's backward rate
                alpha(1,5) = aa(3)*p;       % following motors' forward rate
                alpha(1,6) = aa(4)*q;       % following motors' backward rate
                %Left team
                %Binding
               
                alpha(1,7) = (Nl*Nlane_l-sum(nl))*((M_l-sum(nl))/(M_l))*kon;  % N-n = number of unbound motors  
                alpha(1,8) = sum(nl)*koff;   % n = number of bound motors
                alpha(1,9) = aa(5)*p1_l;    % leading motor's forward rate
                alpha(1,10) = aa(6)*q1_l;    % leading motor's backward rate
                alpha(1,11) = aa(7)*p;       % following motors' forward rate
                alpha(1,12) = aa(8)*q;       % following motors' backward rate
                alpha0 = sum(alpha);  % total rate

                r1= rand(1); r2 = rand(1);   % generate two random number
                % r1 is for time step, r2 is for binding

                tau = (1/alpha0)*(log(1/r1)); % time step

                % twelve possibility happening
                % state = 1 = binding (right),
                % state = 2 = unbindinh (right)
                % state = 3 = leading forward(right),
                % state = 4 = leading backward(right)
                % state = 5 = following forward(right),
                % state = 6 = following backward(right)
                % state = 7 = binding (left),
                % state = 8 = unbindinh (left)
                % state = 9 = leading forward(left),
                % state = 10 = leading backward(left)
                % state = 11 = following forward(left),
                % state = 12 = following backward(left)

                if r2 <= alpha(1,1)/alpha0

                    state=1; %binding of right team
                    [row, column] = find(bound_r == 0); % find the unbound motor

                    if all(bound_r(row(1),:)==0)

                        a=center_r-1;
                        b=center_r+1;
                        
                        if row(1)>1
                            a=max(x_bound_r)-3;
                            b=max(x_bound_r)-1;  
                        end
                        xtobind =  ceil(a + (b-a).*rand(1));

                    else
                        xr_lane=nonzeros(x_r(row(1),:)); % the bound is started from index 1
                        a=min(xr_lane)-3;
                        b=max(xr_lane)-1;

                        if any(bound_r(row(1),2:Nr) == 1) && column(1) == 1 % for the index 1 is rebound again 
                                 a=max(xr_lane)+1;
                                 b=max(xr_lane)+1;            
                        end

                        xtobind =  ceil(a + (b-a).*rand(1));

                        while any(xtobind ==  xr_lane)   % not binding on the same site
                            xtobind =  ceil(a + (b-a).*rand(1));
                            if xtobind ~=  xr_lane
                                break
                            end
                        end


                    end
                    

                     if all(bound_r(:)==0)
                          xr0=xtobind;
                     end

                     bound_r(row(1),column(1)) =1;
                     x_r(row(1),column(1)) = xtobind;
                     nr(row(1)) = nr(row(1)) + 1;

                         

                elseif r2 <= sum(alpha(1:2))/alpha0

                    state =2; % state of unbinding
                    [row, column] = find(bound_r == 1); % find the bound motor

                    if ~isempty(row)
                        bound_r(row(1),column(1)) = 0;
                        x_r(row(1),column(1)) = 0;
                        nr(row(1)) = nr(row(1)) - 1;  % unbinding of right team
                    end

                elseif r2 <= sum(alpha(1:3))/alpha0

                    state = 3; % the leading motor moves forward
                    [row, column] = find(xr_lead_p(1)-1==x_r);
                    x_r(row(1),column(1)) = xr_lead_p(1);


                elseif r2 <= sum(alpha(1:4))/alpha0

                    state = 4; % the leading motor moves backward
                    [row, column] = find(xr_lead_b(1)+1==x_r);
                    x_r(row(1),column(1)) = xr_lead_b(1);


                elseif r2 <= sum(alpha(1:5))/alpha0

                    state = 5; % the following motor moves forward
                    [row, column] = find(xr_fol_p(1)-1==x_r);
                    x_r(row(1),column(1)) = xr_fol_p(1);


                elseif r2 <= sum(alpha(1:6))/alpha0

                    state = 6; % the following motor moves backward
                    [row, column] = find(xr_fol_b(1)+1==x_r);
                    x_r(row(1),column(1)) = xr_fol_b(1);


                elseif r2 <= sum(alpha(1:7))/alpha0

                    state=7; %binding of left team
                    [row, column] = find(bound_l == 0); % find the unbound motor

                    if all(bound_l(row(1),:)==0)
                        a=center_l+1;
                        b=center_l-1;

                        if row(1)>1
                            a=min(x_bound_l)+3;
                            b=min(x_bound_l)+1;  
                        end

                        xtobind =  ceil(a + (b-a).*rand(1));
                    else

                        xl_lane=nonzeros(x_l(row(1),:));

                        a=max(xl_lane)+3; % the bound motor is started from index 1
                        b=min(xl_lane)+1;

                        if any(bound_l(row(1),2:Nl) == 1) && column(1) == 1 % the index 1 is rebound again
                                 a=min(xl_lane)-1;
                                 b=min(xl_lane)-1;            
                        end


                        xtobind =  ceil(a + (b-a).*rand(1));

                        while any(xtobind ==  xl_lane)   % not binding on the same site
                            xtobind =  ceil(a + (b-a).*rand(1));
                            if xtobind ~=  xl_lane
                                break
                            end
                        end
  

                    end

                     if all(bound_l(:)==0)
                          xl0=xtobind;
                     end


                     bound_l(row(1),column(1)) =1;
                     x_l(row(1),column(1)) = xtobind;
                     nl(row(1)) = nl(row(1)) + 1;

                                 

                elseif r2 <= sum(alpha(1:8))/alpha0

                    state =8; % state of unbinding
                    [row, column] = find(bound_l == 1); % find the bound motor

                    if ~isempty(row)
                       % x_del=x_l(row(1),column(1));
                        bound_l(row(1),column(1)) = 0;
                        x_l(row(1),column(1)) = 0;
                        nl(row(1)) = nl(row(1)) - 1;
                    end

                elseif r2 <= sum(alpha(1:9))/alpha0

                    state = 9; % the leading motor moves forward
                    [row, column] = find(xl_lead_p(1)+1==x_l);
                    x_l(row(1),column(1)) = xl_lead_p(1);

                elseif r2 <= sum(alpha(1:10))/alpha0

                    state = 10; % the leading motor moves backward
                    [row, column] = find(xl_lead_b(1)-1==x_l);
                    x_l(row(1),column(1)) = xl_lead_b(1);

                elseif r2 <= sum(alpha(1:11))/alpha0

                    state = 11; % the following motor move forward
                    [row, column] = find(xl_fol_p(1)+1==x_l);
                    x_l(row(1),column(1)) = xl_fol_p(1);

                else

                    state = 12; % the following motor move backward
                    [row, column] = find(xl_fol_b(1)-1==x_l);
                    x_l(row(1),column(1)) = xl_fol_b(1);

                end


                x_bound_r = nonzeros(x_r);
                x_bound_l = nonzeros(x_l);

                if ~isempty(x_bound_r)
                    center_r = (max(x_bound_r)+min(x_bound_r))/2;
                end

                if ~isempty(x_bound_l)
                    center_l = (max(x_bound_l)+min(x_bound_l))/2;
                end


               
                rec_xr(:,:,j)= x_r;
                rec_bound_r(:,:,j)= bound_r;
                rec_bound_l(:,:,j)= bound_l;
                rec_dxr(:,j) = dxr;
                rec_xl(:,:,j)= x_l;
                rec_dxl(:,j) = dxl;
                rec_state(:,j) = state;


                t = t+tau;
                T(1,j) = t;
                rec_dx_move(j) = dxr-dxl; % plus = the center moves to the right
                rec_dx_stretch(j) = dxr+dxl; % minus = the center moves to the left
                rec_fT(j) = f_spring;
                rec_fm(j) = f_motor;
                rec_nr(:,j) = nr;
                rec_nl(:,j) = nl;

                % Find available binding sites for right and left team 
                %right team 
                XR = x_r;
                XR(XR == 0) = NaN;
                b_r = max(XR,[],2);
                a_r = min(XR,[],2);
                M_r_lane = (b_r-a_r)+2;
                M_r_lane(isnan(M_r_lane)) = 20;
                M_r = sum(M_r_lane);

                %left team 
                XL = x_l;
                XL(XL == 0) = NaN;
                b_l = min(XL,[],2);
                a_l = max(XL,[],2);
                M_l_lane = (a_l-b_l)+2;
                M_l_lane(isnan(M_l_lane)) = 20;
                M_l = sum(M_l_lane);
                
                rec_M_r(:,j) = M_r_lane; 
                rec_M_l(:,j) = M_l_lane;
           
            else
                break;
            end

 
        end
        fname = sprintf('laner%d_lanel%d_Nr%dNl%d_%.1f_%d.mat', Nlane_r, Nlane_l, Nr, Nl, k, round);
        save(fname)

    end
end
toc

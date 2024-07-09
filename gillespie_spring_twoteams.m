clear; clc;
% Gillespie Algorithm applying with motors' stepping 
% To compare with the SEP of processive motors

K=0.1;  % spring constant 
[w,e]=size(K);
rate_break = 0.1;  % consider 10^-4 as absolute 0 for V of cargo; 

for d= 1:e
    
    nrun=1;
    
    for round =1:nrun
%% parameters

        nstep =10^6;
        k = K(d);  % spring constant 

        M=10;
        Nr=2;   % total number of motors on right team
        Nl=2;   % total number of motors on left team 
        t=0;   % initial time
        
        center = M/2; % initial center of the cargo
        
        % For following motor
        p=1;    % forward rate
        q=0.1;  % backward rate
        delta=0.5;
        
        % For the leading  motor (the 1st motor)

        x_r = sort(randperm(M, Nr),'descend'); % the array for right team's stepping position at time t
        x_l = -sort(randperm(M, Nl),'descend');  % the array for left team's stepping position at time t
        x_r0=x_r;
        x_l0=x_l;
        x0 = max(x_r)-min(x_l);
        x_bound_r = x_r;
        x_bound_l = x_l;
        %x_bound = cat(2,x_r,x_l); % the array for bound position at time t (not moved yet)
        
        for j=1:nstep
            
            dx=max(x_r)-min(x_l)-x0;
            dxr=max(x_r)-max(x_r0);
            dxl=min(x_l0)-min(x_l);
            ddx=dxr+dxl;
            
            f=k*(ddx);
            
            % Leading motors' rates are keep changing in every step
            p1=p*exp(-f*delta);
            q1=q*exp(f*(1-delta));
            
            % right team 
            
            top = (1-exp(f)*(q/p)^Nr)*(1-q/p);
            bottom = exp(f*delta)*(1-q/p)+exp(f)*(q/p-(q/p)^Nr);
           % format longG
            V_r = p*top/bottom;
           
            % left team 
            
            top = (1-exp(f)*(q/p)^Nl)*(1-q/p);
            bottom = exp(f*delta)*(1-q/p)+exp(f)*(q/p-(q/p)^Nl);
            %format longG
            V_l = p*top/bottom;
            
            
            
            %lim_Vr= V_cargo
            
            if V_r >= rate_break || V_l >= rate_break
                
                rec_fs(j) = f;  % equivalent to stall force
                rec_p1(j) = p1;
                rec_q1(j) = q1;
                rec_Vr(1,j) = V_r;
                rec_Vl(1,j) = V_l;
               
               
                r1= rand(1); r2 = rand(1);   % generate two random number
                % r1 is for time step, r2 is for choosing event happening
                
                % right team 
                if  V_r >= rate_break
                x_fore_r = x_bound_r+1;  % check forward destination
                x_back_r = x_bound_r-1;  % check backward destination
                [x_fore_r,ia]=setdiff(x_fore_r,x_bound_r,'stable');  % possible forward move
                % ia is the index of the value in x_plus before deleting
                % stable =   max ----->   min
                [x_back_r,ib]=setdiff(x_back_r,x_bound_r,'stable'); % possible backward move
                else 
                    ia=0;
                    ib=0;
                end
                
                % left team
                if V_l >= rate_break
                x_fore_l = x_bound_l-1; % forward 
                x_back_l = x_bound_l+1; % backward
                [x_fore_l,ic]=setdiff(x_fore_l,x_bound_l,'stable');  % possible forward move
                [x_back_l,id]=setdiff(x_back_l,x_bound_l,'stable'); % possible backward move
                else
                    ic=0;
                    id=0;
                end
                
                aa = [any(ia==1) any(ib==1) sum(ia>1) sum(ib>1) ...
                    any(ic==1) any(id==1) sum(ic>1) sum(id>1)];
                
                rec_aa(:,j) = aa;
                    
                %Right Team
                alpha(1,1) = aa(1)*p1;    % leading motor's forward rate
                alpha(1,2) = aa(2)*q1;    % leading motor's backward rate
                alpha(1,3) = aa(3)*p;       % following motors' forward rate
                alpha(1,4) = aa(4)*q;       % following motors' backward rate
                %Left Team
                alpha(1,5) = aa(5)*p1;    % leading motor's forward rate
                alpha(1,6) = aa(6)*q1;    % leading motor's backward rate
                alpha(1,7) = aa(7)*p;       % following motors' forward rate
                alpha(1,8) = aa(8)*q;       % following motors' backward rate
                
                
                
                
                
                alpha0 = sum(alpha);  % total rate
                
                tau = (1/alpha0)*(log(1/r1)); % time step
                
                % eight possibility happening
                
                % state = 1 = leading forward(right),
                % state = 2 = leading backward(right)
                % state = 3 = following forward(right),
                % state = 4 = following backward(right)
                % state = 5 = leading forward(left),
                % state = 6 = leading backward(left)
                % state = 7 = following forward(left),
                % state = 8 = following backward(left)
                
                if r2 <= alpha(1,1)/alpha0
                    % the leading motor moves forward
                    state =1;
                    
                    x_r(1) = x_r(1)+1;
                    
                elseif r2 > alpha(1,1)/alpha0  && r2 <= sum(alpha(1:2))/alpha0
                    % the leading motor is moving backward
                    state=2;
                    x_r(1) = x_r(1)-1;
                    
                    
                    
                elseif r2 > sum(alpha(1:2))/alpha0  && r2 <= sum(alpha(1:3))/alpha0
                    %following motor is moving forward
                    state=3;
                    ii=randperm(Nr,1); % ii is the index in the array of x_bound
                    
                    % to choose the index of the following motor
                    while ii > 0
                        ii=randperm(Nr,1);
                        if ii >1 &&  any(ii == ia)  %(not the leading motor)
                            break
                        end
                    end
                    
                    inx3 = find(x_r==x_bound_r(ii));
                    % x_bound = the matrix at the previous step
                    % x_step = the matrix at current step
                    
                    x_r(inx3) = x_r(inx3)+1; % which is equal a value in x_plus
                    
                elseif r2 > sum(alpha(1:3))/alpha0  && r2 <= sum(alpha(1:4))/alpha0
                    %following motor is moving backward
                    state=4;
                    
                    ii=randperm(Nr,1); % ii is the index in the array of x_occupied
                    
                    % to choose the index of the following motor
                    while ii > 0
                        ii=randperm(Nr,1);
                        if ii >1  && any(ii== ib)  %(not the leading motor)
                            break
                        end
                    end
                    
                    inx4= find(x_r==x_bound_r(ii));
                    x_r(inx4) = x_r(inx4)-1; % which is equal a value in x_plus
                    
                elseif r2 > sum(alpha(1:4))/alpha0  && r2 <= sum(alpha(1:5))/alpha0
                    % the leading motor moves forward
                    state =5;
                    
                    x_l(1) = x_l(1)-1;
                elseif r2 > sum(alpha(1:5))/alpha0  && r2 <= sum(alpha(1:6))/alpha0
                    state =6;
                    
                    x_l(1) = x_l(1)+1;
                elseif r2 > sum(alpha(1:6))/alpha0  && r2 <= sum(alpha(1:7))/alpha0
                    
                    state = 7;
                    ii=randperm(Nl,1); % ii is the index in the array of x_bound
                    
                    % to choose the index of the following motor
                    while ii > 0
                        ii=randperm(Nl,1);
                        if ii >1 &&  any(ii == ic)  %(not the leading motor)
                            break
                        end
                    end
                    
                    inx7 = find(x_l==x_bound_l(ii));
                    
                    x_l(inx7) = x_l(inx7)-1; % which is equal a value in x_plus
                else
                     state = 8;
                     ii=randperm(Nl,1); % ii is the index in the array of x_occupied
                    
                    % to choose the index of the following motor
                    while ii > 0
                        ii=randperm(Nl,1);
                        if ii >1  && any(ii== id)  %(not the leading motor)
                            break
                        end
                    end
                    
                    inx8= find(x_l==x_bound_l(ii));
                    x_l(inx8) = x_l(inx8)+1; % which is equal a value in x_plus  
                end
                
                rec_xr(:,j)= x_r;
                rec_xl(:,j)= x_l;
                t = t+tau;
                T(1,j) = t;
                rec_xbound_r(:,j) = x_bound_r;
                rec_xbound_l(:,j) = x_bound_l;
                x_bound_r = x_r;
                x_bound_l = x_l;
                rec_state(:,j) = state;
                rec_dx(:,j) = dx;
                rec_ddx(:,j) = ddx;
                rec_dxr(:,j) = dxr;
                rec_dxl(:,j) = dxl;
            else
                break;
            end
            
            
            
        end
        
        fname = sprintf('sping_cargo_twoteams_%.1f_%d.mat', k, round);
        save(fname)
        
    end
end
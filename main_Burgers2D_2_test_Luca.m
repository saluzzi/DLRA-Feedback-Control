clear
close all

addpath '../Cecilia''s code'

folder = "DLRA/";

strategy = "10_";


Nh1 = 50;
Nh = Nh1^2;
isFOM = 1;
isDiagonal = 1;
isControlled = 1;
isIcare = 0;
isLieGroup = 1;
isSave = 1;
T = 2e-3;
c1 = 20;
c2 = 20;
c = 20;
a = -10;
b = 10;
x = linspace(a,b,Nh1);
dx = x(2)-x(1);
D = sparse((-eye(Nh1)+diag(ones(Nh1-1,1),-1))/dx);
DD = sparse((eye(Nh1)-diag(ones(Nh1-1,1),1))/dx);
D2basic = kron(D,eye(Nh1))+kron(eye(Nh1),D);
[X1,X2] = meshgrid(x,x);
xx = X1(:);
yy = X2(:);
if isDiagonal
    D2 = c2*kron(D,eye(Nh1))+c1*kron(eye(Nh1),D);
else
    D2 = max(-yy,0).*kron(D,eye(Nh1))+max(xx,00).*kron(eye(Nh1),D);
    D2 = D2 + min(-yy,0).*kron(DD,eye(Nh1))+min(xx,00).*kron(eye(Nh1),DD);
    D2 = c*D2;
end
dt = 1.e-3;
nt = floor(T/dt);
alpha = 0.1;
coeff_burger = 1;

np_k = 4;
np_kk = 4;

isAdaptive = 1;
isCheckC_adaptive = 0;
isK = 0;
isDP = 0;
isP0 = 1;

tol_ratio = 0.9999;

%%

tol_trunc = 1e-9;
linesearch = 1;
it_linesearch = 1;
outer_it = 30;

a_B = [-5 -5];
b_B = [0 0];
a_Q = [0 0];
b_Q = [5 5];

vec_B = (xx<=b_B(1)).*(xx>=a_B(1)).*(yy<=b_B(2)).*(yy>=a_B(2));
vec_Q = (xx<=b_Q(1)).*(xx>=a_Q(1)).*(yy<=b_Q(2)).*(yy>=a_Q(2));

m_B = sum(vec_B>0);
m_Q = sum(vec_Q>0);
B = sparse(Nh,m_B);
B(vec_B>0,:) = speye(m_B);
Q = sparse(Nh,Nh);
Q(vec_Q>0,vec_Q>0) =  speye(m_Q)*dx;
R = alpha*speye(m_B)*dx;
invRB = R\B';
By = B*invRB;
A = sparse(D2);
snap_save = 1;

%% tensore %%

TT = sparse(Nh,Nh^2);
TT(1,1) = -2;
for k = 2:Nh1
    TT(k,(k-1)*Nh+k) = -2;
    TT(k,(k-1)*Nh+k-1) = 1;
end

for k = Nh1+1:Nh
    TT(k,(k-1)*Nh+k) = -2;
    if mod(k-1,Nh1) ~= 0
        TT(k,(k-1)*Nh+k-1) = 1;
    end
    TT(k,(k-1)*Nh+k-Nh1) = 1;
end
TT = (coeff_burger/dx)*TT;

%%

%y0fun = @(x,y,k, kk) k*exp(-kk*((x).^2+(y).^2));
y0fun = @(x,y,k, kk) k*exp(-kk*((x+2).^2+(y+2).^2));
k_vec = linspace(0.01,0.05,np_k);
kk_vec = linspace(0.1,0.3,np_kk);
y0 = [];
for i = 1:np_k
    for j = 1:np_kk
        y0 = [y0 y0fun(xx,yy,k_vec(i),kk_vec(j))];
    end
end
np = np_k*np_kk;
s_count = 1;
%% FOM
nrmQ = norm(Q,'fro');
res_A = @(Ay,X) norm(Ay'*X+X*Ay-X*By*X+Q,'fro');
control_FOM = @(Y,Guess) control_full(A,B,Q,R,By,D2basic,Y,Guess,tol_trunc,linesearch,it_linesearch,outer_it,coeff_burger,invRB,isIcare,res_A);
    
if isFOM

    Ny_full = @(Y) coeff_burger*(Y.*(D2basic*Y));
    Fy = @(Y,control) A*Y + Ny_full(Y)+B*control;

    Yj = y0(:,1);
    Y = Yj;

    Y_full = cell(nt,1);
    control = cell(nt-1,1);
    Y_full{1} = y0;
    Y_to_pass = y0;
    Bfull = full(B);
    Qfull = full(Q);
    Rfull = full(R);
    total_cost_FOM = zeros(np,1);
    Pguess = speye(Nh);
    Guess = Pguess;
    P_full = zeros(Nh,Nh,np);
    GuessP_mu0 = Pguess;

    t0 = tic;

    cost_fom = zeros(nt,np);

    for j = 1:np

        cost_fom(1,j) = y0(:,j)'*Q*y0(:,j);

    end

    it_nk_FOM = zeros(nt-1,1);

    for time = 1:nt-1
        time

    %    if mod(time, snap_save) == 0

            Y_full{s_count+1} = zeros(Nh,np);
            control{s_count} = zeros(m_B,np);
        
 %       end

        it_nk_FOM_mu = zeros(np,1);

        for j = 1

            Yj = Y_to_pass(:,j);
            Ay = A+coeff_burger*diag(D2basic*Yj);

            if isControlled

                if mod(j-1,np_kk)==0

                    Guess = Pguess;

                else

                    Guess = P;

                end

                if j == 1

                    Guess = GuessP_mu0;

                end


                [control_j, P, it_nk_FOM_mu(j)] = control_FOM(Yj,Guess);
                res_A(Ay,P)
                
                P_full(:,:,j) = P;

                if mod(j-1,np_kk) == 0

                    Pguess = P;

                end

                if j == 1

                    GuessP_mu0 = P;

                end

            else

                control_j = zeros(m_B,1);

            end
            Y12 = Yj+dt*0.5*Fy(Yj,control_j);
            Ynew_j = Yj + dt*Fy(Y12,control_j);
            cost_fom(time + 1, j) = Yj'*Q*Yj + control_j'*R*control_j;
            Y_to_pass(:,j) = Ynew_j;
        %    if mod(time,snap_save) == 0
            
                control{s_count}(:,j) = control_j;
                Y_full{s_count+1}(:,j) = Ynew_j;
            
            
         %   end
            
        end

        if isControlled
          %  if mod(time,snap_save) == 0
                P_in_time{s_count} = P_full;
                s_count = s_count + 1
           % end
            it_nk_FOM(time) = mean(it_nk_FOM_mu);

        end

    end

    total_cost_FOM = dt*(sum(cost_fom(1:end-1,:),1)+sum(cost_fom(2:end,:),1))/2;
       
        if isSave

            save(folder + strategy + "FOM_sol.mat", "Y_full");
            save(folder + strategy + "FOM_cost.mat", "total_cost_FOM")

            if isControlled

                save(folder + strategy + "FOM_control.mat", "control")
                save(folder + strategy + "FOM_iterations.mat", "it_nk_FOM")

            end

        end
    
    
    full_time_T = toc(t0)
else
    
    if isDiagonal 
		
		load("old_test/" + "DIAGFOM_sol.mat")
		load("old_test/" + "DIAGFOM_control.mat")
	    
	else
		
		load("old_test/" + "ROUNDFOM_sol.mat")
		load("old_test/" + "ROUNDFOM_control.mat")
	end

end


%% DLRA




% NyZ = @(U,Z) coeff_burger*(U.*(D2basic*U))*((Z'.*Z')*Z);
% NyU = @(U,Z) coeff_burger*(Z.*Z)*((U'.*(D2basic*U)')*U);
NyZ = @(U,Z) (TT*(kron(U,U)))*(compute_Zkron(Z)*Z);
NyU = @(U,Z) (U'*TT*(kron(U,U))*compute_Zkron(Z))';



%l = 3;
[U_tot,S_tot,V_tot] = svd(y0);
diagStot = diag(S_tot).^2;
sum_diag = sum(diagStot);
ratio = diagStot(1)/sum_diag;
l_y0 = 1;
while ratio < tol_ratio
    l_y0 = l_y0+1;
    ratio = sum(diagStot(1:l_y0))/sum_diag;
end
l_y0
U = U_tot(:,1:l_y0);


if isAdaptive
    P0 = [];
    for kk = [1]
        [~, P01] = control_FOM(y0(:,kk),zeros(Nh));
        P0 = [P0 P01];
    end
    if isP0
        MM = full(P0);
    elseif isDP
        MM = full(dP);
    else
        K = (-invRB*P0)';
        MM = full(K);
    end
    [Udp2,Sdp2,~] = svd((eye(Nh)-U*U')*MM);
    diagSp2 = diag(Sdp2).^2;
    sum_diag = sum(diagSp2);
    ratio = tol_ratio-1;
    l_dP = 0;
    Z = S_tot*V_tot';
    Z = Z(1:l_y0,:)';
    C = Z'*Z;
    l = l_y0;


    while ratio < tol_ratio && (rank(C) == l || ~isCheckC_adaptive)
        l_dP = l_dP+1;
        U = [U Udp2(:,l_dP)];
        Z = (U'*y0)';
        C = Z'*Z;
        l = l_y0+l_dP;
        ratio = sum(diagSp2(1:l_dP))/sum_diag;
    end

       if rank(C) ~= l && isCheckC_adaptive
        l_dP = l_dP-1;
        l = l-1;
        U = U(:,1:end-1);
        Z = (U'*y0)';
       end

       else
    l = l_y0;
    Z = S_tot*V_tot';
    Z = Z(1:l_y0,:)';
end

U0 = U;
Z0 = Z;

s_count_dlra = 1;

%% RUNGE-KUTTA 2TH ORDER

    U = U0;
    Z = Z0;
    
    F2 = @(U,Z,C,control) (eye(Nh)-U*U')*(A*U+(NyZ(U,Z)+control)*pinv(C));
    Fz2 = @(U,Z,control) Z*(U'*A'*U)+NyU(U,Z)+control;
    
    if isSave
        
        Ymid{1} = U*Z';
    
    end
    
    cost_dlra = zeros(nt,np);
    it_nk_dlra = zeros(nt-1,1);
    
    for mu=1:np
        
        Y_to_use = Ymid{1}(:,mu);
        cost_dlra(1,mu) = Y_to_use'*Q*Y_to_use;
         
    end
    
    Guess_d = eye(l);
    umat = @(U,Z,Guess_d,l) control_proj_T(A,B,Q,R,By,U,Z',tol_trunc,linesearch,it_linesearch,outer_it,np,np_kk,m_B,invRB,isIcare,TT,l, Guess_d);
    t0 = tic;
    rank_dlra = zeros(nt-1,1);
    control_dlra = cell(nt-1,1);

    if ~isControlled

        control_and_P{1} = 0;
        control_and_P{2} = 0;
    
    end

    for time = 1:nt-1
        time

           lvec(time) = l;

        C = Z'*Z;
        rank_dlra(time) = rank(C);
         if isControlled

            control_and_P = umat(U,Z,Guess_d,l);
            % keyboard
            % if mod(time, snap_save)
                 control_dlra{s_count_dlra} = control_and_P{5};
            % end

        else
            % if mod(time, snap_save)
            %     control_dlra{s_count_dlra} = zeros(m_B,np);
            % end

        end

        % dx*norm(control_and_P{5}(:,1)-control{time}(:,1))
        % dx*norm(control_and_P{5}-control{time},'fro')/np
        if isLieGroup
            FF = @(U,Z)  F2(U,Z,C,control_and_P{2});
            U1 = second_uflow(U,Z,dt,FF);
            Z12 = Z+dt*0.5*Fz2(U,Z,control_and_P{1});
            U12 = second_uflow(U,Z,dt/2,FF);
            Z = Z+dt*Fz2(U12,Z12,control_and_P{1});
            U = U1;
        else
            FF = @(U,Z)  F2(U,Z,C,control_and_P{2});
            U12 = U+dt*0.5*FF(U,Z);
            Z12 = Z+dt*0.5*Fz2(U,Z,control_and_P{1});
            Z = Z+dt*Fz2(U12,Z12,control_and_P{1});
            U = U+dt*FF(U12,Z12);
            [U,Rz] = qr(U,'econ');
            Z = Z*Rz';
        end



        if isSave
        
            Y1 = U*Z';
            
         %   if mod(time,snap_save) == 0

                Ymid{s_count_dlra+1} = Y1;
        
            
                if isControlled
            
                    Pdlra{s_count_dlra} = control_and_P{3};
            
                end
                
                s_count_dlra = s_count_dlra + 1
            
          %  end
        
        end
        
        %%%%% compute cost
        
        for mu = 1:np
            
            Y_to_use = Y1(:,mu);
            
        
        
        if isControlled
            control_to_use = control_and_P{5}(:,mu);
            cost_mu = Y_to_use'*Q*Y_to_use + control_to_use'*R*control_to_use;
            cost_dlra(time+1,mu) = cost_mu;
        
        
            Guess_d = control_and_P{4};
            it_nk_dlra(time) = control_and_P{6};
            cost_dlra_time(time) = mean(cost_dlra(time+1,:));
        
        end

        end
       
    end

    total_cost_dlra = dt*(sum(cost_dlra(1:end-1,:),1)+sum(cost_dlra(2:end,:),1))/2;

    dlra_time_T = toc(t0);
    
        if isSave
        
        save(folder + num2str(l) + "_" + strategy + "_T_DLRA_sol.mat", "Ymid")
        save(folder + num2str(l) + "_" + strategy + "_T_DLRA_cost.mat", "total_cost_dlra")
        save(folder + num2str(l) + "_" + strategy + "_T_DLRA_cost_time.mat", "cost_dlra_time")
        save(folder + num2str(l) + "_" + strategy + "_T_DLRA_rank_C.mat", "rank_dlra")
    
        if isControlled
        
            save(folder + num2str(l) + "_" + strategy + "_T_DLRA_control.mat", "control_dlra")
            save(folder + num2str(l) + "_" + strategy + "_T_DLRA_iterations.mat", "it_nk_dlra")
        
        end
         
        err_fro = zeros(s_count_dlra,1);
        err_control = zeros(s_count_dlra,1);
        err_param = zeros(np,1);
        

        for time = 1:s_count_dlra-1
            time
    
            err_fro(time) = sqrt(dx)*norm(Y_full{time} - Ymid{time},'fro');
            err_control(time) = sqrt(dx)*norm(control{time} - control_dlra{time},'fro');
            

            for j=1:np

                err_param(j) = err_param(j) + norm(Y_full{time}(:,j) - Ymid{time}(:,j));
         
            end    
        
        end

        err_param = err_param/nt;
        err_fro = err_fro/np;
        err_control = err_control/np;

        save(folder + num2str(l) + "_" + strategy + "_DLRA_err_fro.mat", "err_fro")
        save(folder + num2str(l) + "_" + strategy + "_DLRA_err_param.mat", "err_param")
        save(folder + num2str(l) + "_" + strategy + "_DLRA_err_control.mat", "err_control")
            

        if isControlled
        
            err_proj_P = zeros(nt-1,1);

            for time = 1:s_count_dlra - 1
            
                dlra_proj_mean = mean(Pdlra{time},3);
                FOM_proj_mean = mean(P_in_time{time},3);
            
                err_proj_P(time) = norm(dlra_proj_mean - FOM_proj_mean,'fro');
            
            end

            save(folder + num2str(l) + "_" + strategy + "_DLRA_err_proj.mat", "err_proj_P")
        
        end
            
    end
        

%%
Nh_plot = 200;
xplot = linspace(a,b,Nh_plot);
[Xq,Yq] = meshgrid(xplot,xplot);
Ymidnt = mean(Ymid{s_count_dlra},2);
Yfullnt = mean(Y_full{s_count_dlra},2);
Ymid_int = interp2(X1,X2,reshape(Ymidnt,Nh1,Nh1),Xq,Yq);
Yfull_int = interp2(X1,X2,reshape(Yfullnt,Nh1,Nh1),Xq,Yq);
close all
max_err_fro = max(err_fro)
max_err_control = max(err_control)
%norm_control
err_cost = mean(abs(total_cost_dlra-total_cost_FOM))
mesh(Xq,Yq,reshape(Ymid_int,Nh_plot,Nh_plot))
figure
mesh(Xq,Yq,reshape(Yfull_int,Nh_plot,Nh_plot))
figure
semilogy(err_fro)

%%
Ymid_int= [];
Yfull_int = [];
for i = 1:nt
    Ymidnt = mean(Ymid{i},2);
    Yfullnt = mean(Y_full{i},2);
    Ymid_int{i} = reshape(interp2(X1,X2,reshape(Ymidnt,Nh1,Nh1),Xq,Yq),Nh_plot,Nh_plot);
    Yfull_int{i} = reshape(interp2(X1,X2,reshape(Yfullnt,Nh1,Nh1),Xq,Yq),Nh_plot,Nh_plot);
end


%%
figure
for i = 1:nt
    subplot(1,2,1)
    contourf(Xq,Yq,Yfull_int{i})
    title('FOM')
    colorbar
    subplot(1,2,2)
    contourf(Xq,Yq,Ymid_int{i})
    title('DLRA')
    colorbar
    drawnow
end
%%
figure
for i = 1:floor(nt/20):nt
    contourf(X1,X2,reshape(mean(abs(Y_full{i}-Ymid{i}),2),Nh1,Nh1))
    title('Error')
    colorbar
    drawnow
end


function Zkron = compute_Zkron(Z)
    
    np = length(Z(:,1));
    r = length(Z(1,:));
    Ztilde = Z';
    Zkron = zeros(r*r,np);
    
    for mu=1:np
        Zkron(:,mu) = kron(Ztilde(:,mu),Ztilde(:,mu));     
    end

end


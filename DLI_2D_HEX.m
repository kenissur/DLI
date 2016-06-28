classdef DLI_2D_HEX < DLI
    %DLI_2D_HEX Solves the DLI model in a 2D hexagonal lattice geometry
    %   The lattice is in the XY plane (?=?/2)
    
    properties
        neighb_ind          % a [n_cells,6] matrix with the indices of the lattice neighbors, sorted clockwise starting from the top neighbor
        n_harms             % total number of spherical harmonics
        h_cells             % number of cells in the horizontal dimension
        v_cells             % number of cells in the vertical dimension
        l_mat               % matrix of l values
        m_mat               % matrix of m values
        motl                % matrix of (-1)^l
        delta_l             % matrix of delta_l0
        sqrt2lplus1         % matrix of sqrt(2l+1)
        Clm                 % spherical harmonics coefficients
    end
    
    methods
        
        % constructor
        function obj=DLI_2D_HEX(varargin)
            p = inputParser;
            boundary_conditions =   {'periodic','sharp'};
            initial_conditions =    {'uniform', 'zero' ,'biased', 'one_higher'};
            
            %set default values
            addParameter(p,'solver',         @ode45,    @(x)isa(x,'function_handle'));% handle to ODE solver to use
            addParameter(p,'h_cells',        8,         @isnumeric);% number of cells in the horizontal dimension
            addParameter(p,'v_cells',        8,         @isnumeric);% number of cells in the vertical dimension
            addParameter(p,'l_max',          10,        @isnumeric);% number of highest spherical harmonic l index representing each distribution
            addParameter(p,'D0',             100,       @isnumeric);% normalized Delta steady-state levels
            addParameter(p,'N0',             100,       @isnumeric);% normalized Notch steady-state levels
            addParameter(p,'DiffD',          0.1,       @isnumeric);% normalized Delta diffusion rate
            addParameter(p,'DiffN',          0.1,       @isnumeric);% normalized Notch diffusion rate
            addParameter(p,'k',              0.001,     @isnumeric);% normalized interaction strength
            addParameter(p,'alpha',          1,         @isnumeric);% ?_D / ?_N
            addParameter(p,'t_max',          100,       @isnumeric);% simulation time
            % addParameter(p,'cells_to_plot',  3,         @isnumeric);
            addParameter(p,'N_theta',        100,       @isnumeric);
            addParameter(p,'N_phi',          200,       @isnumeric);
            addParameter(p,'BC',             'periodic',@(x) any(validatestring(x,boundary_conditions)));  %boundary conditions
            addParameter(p,'IC',             'uniform', @(x) any(validatestring(x,initial_conditions)));   %initial conditions
            parse(p,varargin{:});
            params = p.Results;
            
            %work out additional parameters from the given input
            params.tspan = [0,params.t_max];                            %simulation time span
            params.n_cells = params.h_cells * params.v_cells;           %total number of cells
            params.n_harms = (1+(params.l_max*2+1))*(params.l_max+1)/2; %total number of spherical harmonics
            
            %matricize kinetic parameters
            params.D0 =     DLI_2D_HEX.matricize(params.D0,    params.n_cells, params.n_harms);
            params.N0 =     DLI_2D_HEX.matricize(params.N0,    params.n_cells, params.n_harms);
            params.DiffD =  DLI_2D_HEX.matricize(params.DiffD, params.n_cells, params.n_harms);
            params.DiffN =  DLI_2D_HEX.matricize(params.DiffN, params.n_cells, params.n_harms);
            
            %sph harm coeffs
            % sph_coeff_fun = @(l,m)sqrt((2*l+1)/(4*pi));                           %this is the spherical harmonic normalization coefficient when m=0
            % params.sph_coeff = sph_coeff_fun(params.l_mat);                     %spherical harmonic norm coeffs for m=0
            l_vec=[];m_vec=[];
            for l=0:params.l_max
                l_vec=[l_vec,l*ones(1,2*l+1)];                                  %#ok<AGROW> %generate vector of l values
                m_vec=[m_vec,-l:l];                                             %#ok<AGROW> %generate vector of l values
            end
            params.m_mat = repmat(m_vec,params.n_cells,1);                      %matrix of m values
            params.l_mat = repmat(l_vec,params.n_cells,1);                      %matrix of l values
            params.motl = (-1).^params.l_mat;                                   %matrix of (-1)^l
            params.delta_l = params.l_mat==0;                                   %matrix of ?_l?_m
            params.sqrt2lplus1 = sqrt(2*params.l_mat+1);                        %matrix of sqrt(2l+1)
            params.Clm = sqrt(2*params.l_mat+1.*factorial(params.l_mat-...      %Clm coefficient of the spherical harmonic normalization
                params.m_mat)./factorial(params.l_mat+params.m_mat));
            
            %boundary conditions
            if(strcmp(params.BC,'periodic'))
                %hexagonal lattice - the matrix neighb_ind will contain the indices the
                %the neighbors of the ith cell in the ith row, ordered clockwise as:
                %up,up-right,up-left,down,down-left,up-left
                params.neighb_ind = zeros(params.n_cells,6);
                for j = 1:params.n_cells
                    params.neighb_ind(j,:) = DLI_2D_HEX.findneighborhex(j,params.v_cells,params.h_cells);
                end
            end
            
            % copy from params to obj
            p = properties(obj);
            fn = fieldnames(params);
            for i = 1:length(p)                 % iterate properties
                if(ismember(p{i},fn))           % if a property is a field in params
                    obj.(p{i}) = params.(p{i}); % copy from params to obj property
                end
            end
            
            % set results to empty values
            obj.t = [];
            obj.y_de = [];
            obj.y_ss = [];
        end
        
        
        
        % get time derivatives vector
        function dy = dy(obj,t,y)
            %get N and D in matrix format
            [N,D] = obj.y2nd(y);
            % get interaction terms
            [Lambda_N_ijlm, Lambda_D_ijlm] = obj.get_lambda(y);
            % calculate differentials
            dN =  obj.delta_l - N./obj.N0 - obj.DiffN .* obj.l_mat .* (obj.l_mat+1) .* N - obj.k .* Lambda_N_ijlm;
            dD = (obj.delta_l - D./obj.D0 - obj.DiffD .* obj.l_mat .* (obj.l_mat+1) .* D - obj.k .* Lambda_D_ijlm) ./ obj.alpha;
            % discars imaginary part (it is always small from what I tested)
            dN = real(dN);
            dD = real(dD);
            % return equations in vector form
            dy = obj.nd2y(dN,dD);
        end
        
        % get Lambda interaction terms
        function [Lambda_N_ijlm, Lambda_D_ijlm] = get_lambda(obj,y)
            %Now we calculate the concentrations of N and D on the interfaces
            theta = pi/2;                                                   %the interface is in the ?=?/2
            phi = pi/3*(0:5);                                               %hexagonal lattice
            [N_raw,D_raw,ylm_raw] = obj.ND_distribution(...
                y,'theta',theta,'phi',phi);                          %N,D concentrations on the interfaces, [n_cells,1,6]
            ylm_tmp = repmat(permute(ylm_raw,[4,1,2,3]),1,obj.n_cells,1);       % take CC and permute dimensions to the shape [6,n_cells,n_harms]
            N_dist = permute(N_raw,[3,1,2]);                                    % permute dimensions to the shape [6,n_cells]
            D_dist = permute(D_raw,[3,1,2]);                                    % permute dimensions to the shape [6,n_cells]
            % in this part we will calculate the indices to work out the
            % neighbor multiplications of the form Nij*Dpq
            ij_index = repmat((1:obj.n_cells),6,1);                             % the index of the "(i,j)" cell, the one we are "looking at"
            pq_index = obj.neighb_ind';                                         % the index of the neighbor cells, sorted clockwise
            int_index = repmat(1:6,1,obj.n_cells)';                             % the order of the interfaces - [u,ur,dr,d,dl,ul] (u=up,d=down,r=right,l=left)
            shifted_int_index = repmat(circshift(1:6,3,2),1,obj.n_cells)';      % the order of the interfaces shifted by 180deg - [d,dl,ul,u,ur,dr] (u=up,d=down,r=right,l=left)
            my_indices = sub2ind(size(N_dist),int_index(:),ij_index(:));        % the indices of the ij cell sorted clockwise ("my index")
            neighb_indices = sub2ind(size(N_dist),shifted_int_index(:),pq_index(:)); % the indices of the pq cell (neighbors to ij) sorted clockwise but shifted by 180deg (the circshift of the 1:6 array does that)
            NijDpq = reshape( N_dist(my_indices) .* D_dist(neighb_indices),6,obj.n_cells); % the Nij(phi_0)Dpq(phi_0+pi) terms in the sum, shaped to [6,n_cells]
            DijNpq = reshape( D_dist(my_indices) .* N_dist(neighb_indices),6,obj.n_cells); % the Dij(phi_0)Npq(phi_0+pi) terms in the sum, shaped to [6,n_cells]
            
            % now we calculate the lambda interaction factors of the form "sigma-Nij*Dpq*Ylm "
            % LaTeX: \Lambda^N_{ijlm}=\sum_{(p,q,\phi_0)}N_{ij}(\frac{\pi}{2},\phi_0)D_{pq}(\frac{\pi}{2},\phi_0+\pi)Y^*_{lm}(\frac{\pi}{2},\phi_0)
            Lambda_N_ijlm = squeeze(sum(...                                      % we will sum over the first dimension of 6 neighbors
                repmat(NijDpq,1,1,obj.n_harms) .* ...                           % we extend the dimensions of this expression to n_harms, it does not depend on lm
                ylm_tmp,1));                                                    % multiply by Y*lm and we are finished
            Lambda_D_ijlm = squeeze(sum(...                                      % we will sum over the first dimension of 6 neighbors
                repmat(DijNpq,1,1,obj.n_harms) .* ...                           % we extend the dimensions of this expression to n_harms, it does not depend on lm
                ylm_tmp,1));                                                    % multiply by Y*lm and we are finished
%             % this part was used until 28.6.16, I think it has a mistkae
%             N_dist = repmat(permute(N_raw,[1,4,2,3]),1,obj.n_harms,1,1);        % extend dims to [n_cells,n_harms,1,6]
%             D_dist = repmat(permute(D_raw,[1,4,2,3]),1,obj.n_harms,1,1);        % extend dims to [n_cells,n_harms,1,6]
%             ylm_cc = repmat(conj(ylm_raw),obj.n_cells * 6,1,1,1);               % take CC and extend dims to [n_cells*6,n_harms,1,6]
%             %both the pq and ij indices are shaped so that when we take (:) we get each 6 consecutive indices matching the same cell
%             %Now we calculate the Lambda factors of the interaction
%             %What is done in those lines is that we spread the indices along the first
%             %dim (n_cells) so each 6 consecutive indices are ij indices of the same
%             %cells and the different pq indices of the neighbours. We take the product
%             %of N*D*Y_lm and then reshape to another dimension, so we get
%             %[6,n_cells,n_harms,1,6]. Then we sum over the first and last dimension and
%             %squeeze - that is the (p,q,phi_0) sum from the derivation (15c,15d)
%             Lambda_N_ijlm =squeeze(sum(sum(reshape(...
%                 N_dist(ij_index(:),:,:,:) .* ...        %this part is N_ij(theta=pi/2,phi_0)
%                 circshift(D_dist(pq_index(:),:,:,:),3,4) .* ...   %this part is D_pq(theta=pi/2,phi_0 + pi) - the circular shift does the same as adding pi to the phi values. The arguments are: 3 - by how much the matrix is shifted, 4 - the dimension of the shift
%                 ylm_cc...                               %times Y*_lm. The summation is over the pq indices
%                 ,[6,obj.n_cells,obj.n_harms,1,6]),5),1));       %those are the reshape and sum arguments
%             
%             Lambda_D_ijlm =squeeze(sum(sum(reshape(...
%                 D_dist(ij_index(:),:,:,:) .* ...        %this part is D_ij(theta=pi/2,phi_0)
%                 circshift(N_dist(pq_index(:),:,:,:),3,4) .* ...   %this part is N_pq(theta=pi/2,phi_0 + pi) - the circular shift does the same as adding pi to the phi values. The arguments are: 3 - by how much the matrix is shifted, 4 - the dimension of the shift
%                 ylm_cc...                               %times Y*_lm. The summation is over the pq indices
%                 ,[6,obj.n_cells,obj.n_harms,1,6]),5),1));       %those are the reshape and sum arguments
%             
        end
        
        % get initial conditions vector
        function y_ic = get_ic(obj)
            %return y vector of initial conditions.
            %options - 'uniform' or 'zero' or 'biased' or 'one_higher'
            N = zeros(obj.n_cells, obj.n_harms);
            D = zeros(obj.n_cells, obj.n_harms);
            switch(obj.IC)
                case 'zero'
                case 'uniform'%setting the first harmonic in each cell to the steady state value (when k=D=0)
                    N = obj.N0 .* obj.delta_l / 3;  % division by 3 is an approximation to the steady state concentrations in the presence of interactions, otherwise it takes a long time for the solution to reach steady state
                    D = obj.D0 .* obj.delta_l / 3;
                case 'one_higher'%the first cell has 10% more N and D than the steady state, the rest are uniform
                    N = obj.N0 .* obj.delta_l / 2;
                    D = obj.D0 .* obj.delta_l / 2;
                    N(1,1) = N(1,1) * 1.1;
                    N(1,1) = N(1,1) * 1.1;
                case 'biased'%setting the l=1,m=-1 harmonic to be nonzero so we get bias
                    N = obj.N0 .* obj.delta_l / 3;
                    N(2,:) = N(1,:) * 0.5;
                    D = obj.D0 .* obj.delta_l / 3;
                    D(2,:) = D(1,:) * (-0.5);
                otherwise
                    Exception('invalid mode')
            end
            y_ic = obj.nd2y(N,D);
        end
        
        % convert N,D matrices to vector
        function [N,D] = y2nd(obj,y)
            %y2nd  converts the long y vector of differential equations to Notch and Delta
            %matrices, so that N(i,l) is the (l-1)-th spherical harmonic in the i-th cell
            N_all = y(1:obj.n_cells * obj.n_harms);
            D_all = y(obj.n_cells * obj.n_harms + 1 : 2* obj.n_cells * obj.n_harms);
            N = reshape(N_all,obj.n_harms,obj.n_cells)';
            D = reshape(D_all,obj.n_harms,obj.n_cells)';
        end
        
        % convert vector to N,D matrices
        function y = nd2y(~,N,D)
            %nd2y  converts N,D matrices to long vector format
            y = [reshape(N',1,[]),reshape(D',1,[])]';
        end
        
        % get angular distributions from coefficients
        function [N_dist,D_dist,ylm] = ND_distribution(obj,y,varargin)
            %ND_distribution  converts the N and D coefficients to functions of theta
            %ylm, the spherical harmonic itself on the corresponding theta,phi values, is also returned
            %varargin can contain a 'theta',theta_vector and/or 'phi',phi_vector
            p = inputParser;
            addParameter(p,'theta', linspace(0,pi,obj.N_theta), @isnumeric);% theta values
            addParameter(p,'phi',   linspace(0,2*pi,obj.N_phi), @isnumeric);% phi values
            parse(p,varargin{:});
            vparams = p.Results;
            
            [N,D] = obj.y2nd(y);
            theta =   vparams.theta;
            phi =     vparams.phi;
            N_theta = length(vparams.theta);
            N_phi =   length(vparams.phi);
            ylm = zeros(1,obj.n_harms,N_theta,N_phi);                      %[cells,harms,theta,phi]
            counter=1;                                                 %counter is the index of the next value to be filled
            for l = 0:obj.l_max
                ylm(1,counter:counter+2*l,:,:) = DLI.Y_lm(l,theta,phi,[]); %[m,theta,phi]
                counter = counter + 2*l+1;
            end
            ylm_mat = repmat(ylm,obj.n_cells,1,1,1);
            N_mat = repmat(N,1,1,N_theta,N_phi);                            %strech the dimensions to cover theta
            D_mat = repmat(D,1,1,N_theta,N_phi);                            %strech the dimensions to cover theta
            %after the two following lines we are left with dimensions [n_cells,n_theta,n_phi]
            N_dist = sum(N_mat .* ylm_mat,2);                               %calculate the spherical harmonics sum
            D_dist = sum(D_mat .* ylm_mat,2);                               %calculate the spherical harmonics sum
            sz = size(N_dist);
            N_dist = reshape(N_dist,[sz(1),sz(3),sz(4)]);                   %squeeze out the second dimension only (n_harms)
            D_dist = reshape(D_dist,[sz(1),sz(3),sz(4)]);                   %squeeze out the second dimension only (n_harms)
        end
        
        % plot angular distributions
        function ax = plot_distributions(obj)
            % plot_distributions  TODO DOC
            ss_exists = ~isempty(obj.y_ss); %steady state solution exists
            de_exists = ~isempty(obj.y_de); %ODE solution exists
            if(~ss_exists && ~de_exists)
                throw(MException('DLI_2D_HEX:plot_distributions','Attempted to plot before getting results'));
            end
            DLI.fullscreen_figure();
            
            plot_times = ceil([1,find(obj.t>0.1*obj.t(end),1),find(obj.t>0.5*obj.t(end),1),length(obj.t)]);%times where we will plot the distributions
            i_mid = round(obj.N_theta/2);
            theta_vec = linspace(0,pi,obj.N_theta);
            phi_vec = linspace(0,2*pi,obj.N_phi);
            
            %get steady state
            if(ss_exists)
                [N_dist,D_dist] = obj.ND_distribution(obj.y_ss,'theta',theta_vec,'phi',phi_vec);
                Nss = squeeze(real(N_dist(1,:,:)));
                Dss = squeeze(real(D_dist(1,:,:)));
            end
            
            %plot distributions
            counter=0;
            ax = zeros(1,length(plot_times));
            for it = plot_times
                counter = counter+1;
                ax(counter) = subplot(2,2,counter);
                box on,hold on
                if(ss_exists)
                    hN = polar(phi_vec,Nss(i_mid,:),  '--b'); %  hN.LineWidth=1;hN.Color(4)=0.2;
                    hD = polar(phi_vec,Dss(i_mid,:),  '--r'); %  hD.LineWidth=1;hD.Color(4)=0.2;
                end
                if(de_exists)
                    [N_dist,D_dist] = obj.ND_distribution(obj.y_de(it,:));
                    N_vec = squeeze(real(N_dist(1,:,:)));
                    D_vec = squeeze(real(D_dist(1,:,:)));
                    hN = polar(phi_vec,N_vec(i_mid,:),'b');   %  hN.LineWidth=1;hN.Color(4)=0.6;
                    hD = polar(phi_vec,D_vec(i_mid,:),'r');   %  hD.LineWidth=1;hD.Color(4)=0.6;
                end
                title(sprintf('t=%.1f',obj.t(it)),'FontSize',18)
                set(gca,'FontSize',14)
                axis equal
            end
            drawnow
            
            % 3D plots of the steady state solution
            hfig = DLI.fullscreen_figure();
            if(ss_exists)
                N_plot = real(Nss);
                D_plot = real(Dss);
                set(hfig,'name','Steady State')
            else
                N_plot = real(squeeze(N_dist(1,:,:)));
                D_plot = real(squeeze(D_dist(1,:,:)));
                set(hfig,'name',sprintf('end of ODE solution (t=%d)',obj.t(end)))
            end
            ax2(1) = subplot(2,2,1);hold on
            hs=surf(sin(theta_vec')*sin(phi_vec),sin(theta_vec')*cos(phi_vec),cos(theta_vec')*ones(1,obj.N_phi),N_plot);
            hs.LineStyle='none';
            axis equal,box on;colorbar,view(15,30),title('Notch')
            ax2(3) = subplot(2,2,3);hold on
            contour(sin(theta_vec')*sin(phi_vec),sin(theta_vec')*cos(phi_vec),N_plot,15)
            axis equal,box on;colorbar,title('Notch')
            ax2(2) = subplot(2,2,2); hold on
            hs = surf(sin(theta_vec')*sin(phi_vec),sin(theta_vec')*cos(phi_vec),cos(theta_vec')*ones(1,obj.N_phi),D_plot);
            hs.LineStyle='none';
            axis equal,box on;colorbar,view(15,30),title('Delta')
            ax2(4) = subplot(2,2,4);hold on
            contour(sin(theta_vec')*sin(phi_vec),sin(theta_vec')*cos(phi_vec),D_plot,15)
            axis equal,box on;colorbar,title('Delta')
            ax = [ax,ax2]; %append axes for the return value
            
        end
        
        % plot total number of molecules
        function ax = plot_total(obj)
            ss_exists = ~isempty(obj.y_ss); %steady state solution exists
            de_exists = ~isempty(obj.y_de); %ODE solution exists
            if(~ss_exists && ~de_exists)
                throw(MException('DLI_2D_HEX:plot_distributions','Attempted to plot before getting results'));
            end
            
            if(de_exists)
                t_ind = ceil([1,find(obj.t>0.1*obj.t(end),1),find(obj.t>0.5*obj.t(end),1),length(obj.t)]);
                %? times where we will plot the distributions
                t_plot = obj.t(t_ind);
                Ntot = zeros(length(t_plot),obj.n_cells);             %allocate arrays
                Dtot = zeros(length(t_plot),obj.n_cells);
                for j = 1:length(t_plot)                                 %iterate over t samples and perform integrations to obtain total amounts
                    [N_dist,D_dist] = obj.ND_distribution(obj.y_de(t_ind(j),:)) ;
                    Ntot(j,:) = DLI_2D_HEX.integral_func(real(N_dist));
                    Dtot(j,:) = DLI_2D_HEX.integral_func(real(D_dist));
                end
                
                max_Ntot = max(Ntot(:));
                min_Ntot = min(Ntot(:));
                
                if(any(Ntot(:)<0) || any(Dtot(:)<0))                %skip this plot if there are negative values in N, which is an indication of ODE solver problems
                    disp('Negative values in N distribution - skipping total_plot')
                    ax = [];
                    return
                end
            end
            
            % get steady state total numbers
            if(ss_exists)
                [N_dist,D_dist] = obj.ND_distribution(obj.y_ss);
                Nss = DLI_2D_HEX.integral_func(real(N_dist));                        %steady state
                Dss = DLI_2D_HEX.integral_func(real(D_dist));                        %steady state
            end
            
            % create figure for lattice plot
            hFig=figure;drawnow;warning off;jFig = get(handle(hFig), 'JavaFrame');jFig.setMaximized(true);warning on; %#ok<WNON,WNOFF>
            ax = zeros(1,5);
            ax(1)=subplot(3,2,1);       set(gca,'FontSize',14);hold on;box on
            ax(2)=subplot(3,2,2);       set(gca,'FontSize',14);hold on;box on
            ax(3)=subplot(3,2,3);       set(gca,'FontSize',14);hold on;box on
            ax(4)=subplot(3,2,4);       set(gca,'FontSize',14);hold on;box on
            ax(5)=subplot(3,2,[5,6]);   set(gca,'FontSize',14);hold on;box on
            for i = 1:4
                set(ax(i),'XTick',[],'YTick',[],'FontSize',14);
                title(ax(i),sprintf('t=%.1f',t_plot(i)))
            end
            title(ax(5),'SS')
            set(ax(5),'XTick',[],'YTick',[],'FontSize',14)
            
            % do_plot is a function that plots to ax(k) circle #j
            do_plot = @(k,j)DLI_2D_HEX.plotCircle(...
                ax(k),j,(Ntot(k,j)-min_Ntot+0.05)/max_Ntot,obj.v_cells);
            for j = 1:obj.n_cells   %time
                if(de_exists)
                    do_plot(1,j);
                    do_plot(2,j);
                    do_plot(3,j);
                    do_plot(4,j);
                end
                if(ss_exists)
                    DLI_2D_HEX.plotCircle(ax(5),j,(Nss(j)-min(Nss)+0.05)/max(Nss(:)),obj.v_cells);
                end
            end
            
            % create figure for concentration vs. time plot
            if(de_exists)
                hFig=figure;drawnow;warning off;jFig = get(handle(hFig), 'JavaFrame');jFig.setMaximized(true);warning on; %#ok<WNON,WNOFF>
                t_ind = round(linspace(1,length(obj.t),100));%we only work with 100 time samples in this function
                t_plot = obj.t(t_ind);
                for j = 1:length(t_plot)                                 %iterate over t samples and perform integrations to obtain total amounts
                    [N_dist,D_dist] = obj.ND_distribution(obj.y_de(t_ind(j),:)) ;
                    Ntot(j,:) = DLI_2D_HEX.integral_func(real(N_dist));
                    Dtot(j,:) = DLI_2D_HEX.integral_func(real(D_dist));
                end
                ax1=subplot(2,1,1);hold on,box on,set(gca,'FontSize',14)
                if(ss_exists)
                    plot([t_plot(1) t_plot(end)],mean(Nss) * [1 1],'--b','LineWidth',2)
                end
                plot(t_plot,Ntot,'-.')
                ylabel('Total Notch in cell')
                xlabel('\tau')
                ax2=subplot(2,1,2);hold on,box on,set(gca,'FontSize',14)
                if(ss_exists)
                    plot([t_plot(1) t_plot(end)],mean(Dss) * [1 1],'--r','LineWidth',2)
                end
                plot(t_plot,Dtot,'-.')
                ylabel('Total Delta in cell')
                xlabel('\tau')
                linkaxes([ax1 ax2],'x')
                ax = [ax,[ax1,ax2]]; %append those axes for the return value
            end
            
        end
       
        % plot a display of the steady state Nijlm coefficients
        function ax = plot_coeffs(obj,varargin)
            ss_exists = ~isempty(obj.y_ss); %steady state solution exists
            de_exists = ~isempty(obj.y_de); %ODE solution exists
            if(~ss_exists && ~de_exists)
                throw(MException('DLI_2D_HEX:plot_coeffs','Attempted to plot before getting results'));
            end
            if(ss_exists)
                [N,~] = obj.y2nd(obj.y_ss);
            else
                [N,~] = obj.y2nd(squeeze(obj.y_de(end,:)));
            end
            if(~isempty(varargin) && ishandle(varargin{1}))
                ax = varargin{1};
            else
                DLI.fullscreen_figure();
                ax = gca;
            end
            scatter3(ax,obj.l_mat(1,:),obj.m_mat(1,:),real(N(1,:)),36,real(N(1,:)),'filled','marker','s');
            box on
            view(0,90)
            colorbar
            axis equal
            xlabel('l')
            ylabel('m')
            if(ss_exists)
                title('N_1_1_l_m (steady state)')
            else
                title('N_1_1_l_m (end of ODE solution)')
            end
        end
        
    end
    
    methods(Static)
        % convert input to matrix form
        function y = matricize(x,c,h)
            %matricize  make sure the output is a a (n_cells,n_harms) sized matrix, whether the
            %input is a matrix or a scalar
            %%% converts x to 2D matrix, whether x is scalar or matrix of the right size
            if(length(x)==1)                    %scalar
                y = x * ones(c,h);
            elseif(size(x,1)==c && size(x,2)==h)%mat
                y = x;
            elseif(size(x,1)==h && size(x,2)==c)%transposed mat
                y = x';
            elseif(size(x,1)==1 && size(x,2)==c)%transposed vector
                y = repmat(x',1,h);
            elseif(size(x,1)==c && size(x,2)==1)%vector
                y = repmat(x,1,h);
            else
                throw(MException('BadInputDimensions'))
            end
        end
        
        % perform an integral on an angular distribution
        function integral_res = integral_func(X_dist)
            %integral_func   performs a theta and phi intergartion on an input of the
            %form [N_cells,N_theta,N_phi], assuming it is uniformly sampled in
            %theta and phi
            [N_cells,N_theta,N_phi] = size(X_dist);
            theta = linspace(0,pi,N_theta);
            % phi = linspace(0,2*pi,N_phi);
            sintheta = repmat(sin(theta),N_cells,1,N_phi);
            dtheta = pi/N_theta;
            dphi = 2*pi/N_phi;
            integral_res = reshape( trapz(trapz(X_dist.*sintheta,3),2)*dtheta*dphi , 1,N_cells);
        end
        
        % hexagon indexing function adjusted from David Sprinzak's code
        function out = findneighborhex(ind,P,Q)
        % findneighborhex  returns the indices of the 6 neighbors in a
        % hexagonal lattice, sorted clockwise.
            % This function finds the 6 neighbors of cell ind
            [p,q] = DLI_2D_HEX.ind2pq(ind,P);
            
            % above and below:
            out(1) = DLI_2D_HEX.pq2ind(mod(p,P)+1,q,P);
            out(4) = DLI_2D_HEX.pq2ind(mod(p-2,P)+1,q,P);
            
            % left and right sides:
            qleft = mod(q-2,Q)+1;
            qright = mod(q,Q)+1;
            
            if q/2~=round(q/2),
                pup = p;
                pdown = mod(p-2,P)+1;
            else
                pup = mod(p,P)+1;
                pdown = p;
            end;
            out(6) = DLI_2D_HEX.pq2ind(pup,qleft,P);
            out(5) = DLI_2D_HEX.pq2ind(pdown,qleft,P);
            out(2) = DLI_2D_HEX.pq2ind(pup,qright,P);
            out(3) = DLI_2D_HEX.pq2ind(pdown,qright,P);
        end
        
        % convert single index to (p,q) indices
        function ind=pq2ind(p,q, P)
            ind = p + (q-1)*P;
        end
        
        % convert (p,q) indices to single index
        function [p,q]=ind2pq(ind, P)
            q = 1+floor((ind-1)/P);
            p = ind - (q-1)*P;
        end
        
        % get the geometric center of the (p,q) hex
        function [pc,qc]=get_center(ind,P)
        % get_center  returns the geometric center of the (p,q) hex
            
            %I copied this function from David's code without really understanding the
            %computation. There is a smart way to determine the center but it's
            %unimportant
            [p0,q0]=DLI_2D_HEX.ind2pq(ind, P);
            s32 = sqrt(3)/4;
            q = q0*3/4;
            p = p0*2*s32;
            if q0/2 == round(q0/2),
                p = p+s32;
            end;
            
            x(1)=q-.5; x(2)=q-.25; x(3)=q+.25;
            x(4)=q+.5; x(5)=q+.25; x(6)=q-.25;
            
            y(1)=p ; y(2)=p+s32; y(3)=p+s32;
            y(4)=p; y(5)=p-s32; y(6)=p-s32;
            pc=mean(x);
            qc=mean(y);
        end
        
        % plot a circle representing the ind cell
        function plotCircle(ax,ind,c,v_cells)
        % plotCircle  plots to ax a circle representing the ind cell
        % c is a scalar between 0 and 1
            [p0,q0]=DLI_2D_HEX.get_center(ind, v_cells);               %get center coordinates    
            % if(length(c)==1),c=c*ones(1,3);end              %if c is a scalar treat it as greyscale
            axes(ax);                                       %switch axes
            r = sqrt(3)/2;                                  %cell radius in the display
            rectangle('Position',[p0,q0,r,r],'Curvature',...%plot the circle
                [1 1],'FaceColor',[1-c,1-c,1],'EdgeColor','k','LineWidth',1)
            axis equal
        end
        
        
    end
    
end
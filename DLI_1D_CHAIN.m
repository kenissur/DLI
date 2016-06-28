classdef DLI_1D_CHAIN < DLI
    %DLI_1D Solves the DLI model in a 1D chain geometry
    %   Using only m=0 harmonics, neighboring cells touch only in ?=0 or ?=?
    
    properties
        cells_to_plot       % determines how many cells are plotted in plot_distributions
        i_up                % vector of the index of the up neighbor to the ith cell
        i_down              % vector of the index of the down neighbor to the ith cell
        l_mat               % matrix of l values
        motl                % matrix of (-1)^l
        delta_l             % matrix of ?_l 
        sqrt2lplus1         % matrix of sqrt(2l+1)
        sph_coeff           % matrix of the spherical harmonic coefficient
    end
    
    methods
        % constructor
        function obj = DLI_1D_CHAIN(varargin)
            % DLI_1D_CHAIN is a constructor for the class
            %   obj = DLI_1D_CHAIN() creates an object with the default values
            %
            %   obj = DLI_1D_CHAIN('param',value) creates an object with the
            %   default values, except the field 'param' is replaced with value
            %
            %   Input parameters:
            p = inputParser;
            boundary_conditions =   {'periodic','sharp'};
            initial_conditions =    {'uniform', 'zero' ,'biased', 'one_higher'};
            
            %set default values
            addParameter(p,'n_cells',        10,         @isnumeric);% number of cells in the chain
            addParameter(p,'l_max',          30,         @isnumeric);% number of spherical harmonics representing each distribution
            addParameter(p,'D0',             100,        @isnumeric);% normalized Delta steady-state levels
            addParameter(p,'N0',             100,        @isnumeric);% normalized Notch steady-state levels
            addParameter(p,'alpha',          1,          @isnumeric);% the ration between notch and delta production, ?_D/?_N
            addParameter(p,'DiffD',          0.1,        @isnumeric);% normalized Delta diffusion rate
            addParameter(p,'DiffN',          0.1,        @isnumeric);% normalized Notch diffusion rate
            addParameter(p,'k',              0.001,      @isnumeric);% normalized interaction strength
            addParameter(p,'t_max',          100,        @isnumeric);% simulation time
            addParameter(p,'cells_to_plot',  3,          @isnumeric);
            addParameter(p,'N_theta',        100,        @isnumeric);
            addParameter(p,'BC',             'periodic', @(x) any(validatestring(x,boundary_conditions)));  %boundary conditions
            addParameter(p,'IC',             'uniform',  @(x) any(validatestring(x,initial_conditions)));   %initial conditions
            addParameter(p,'solver',         @ode45,     @(x)isa(x,'function_handle'));% handle to ODE solver to use
            parse(p,varargin{:});
            params = p.Results;
            %vectorize kinetic parameters
            params.D0 =     DLI_1D_CHAIN.vectorize(params.D0,    params.l_max, params.n_cells);
            params.N0 =     DLI_1D_CHAIN.vectorize(params.N0,    params.l_max, params.n_cells);
            params.DiffD =  DLI_1D_CHAIN.vectorize(params.DiffD, params.l_max, params.n_cells);
            params.DiffN =  DLI_1D_CHAIN.vectorize(params.DiffN, params.l_max, params.n_cells);
            params.alpha =  DLI_1D_CHAIN.vectorize(params.alpha, params.l_max, params.n_cells);
            
            %work out additional parameters from the given input
            %cell to plot
            params.tspan = [0,params.t_max];
            params.cells_to_plot = min(3,params.n_cells);
            %sph harm coeffs and other matrices
            sph_coeff_fun = @(l)sqrt((2*l+1)/(4*pi));                           %this is the spherical harmonic normalization coefficient when m=0
            params.l_mat = repmat((0:params.l_max-1)',1,params.n_cells);        %matrix of l values
            params.motl = (-1).^params.l_mat;                                   %matrix of (-1)^l
            params.delta_l = params.l_mat==0;                                   %matrix of ?_l (he integral over one spherical harmonic is ?_l0*?_m0     source: http://functions.wolfram.com/Polynomials/SphericalHarmonicY/21/02/02/0001/
            params.sqrt2lplus1 = sqrt(2*params.l_mat+1);                        %matrix of sqrt(2l+1)
            params.sph_coeff = sph_coeff_fun(params.l_mat);                     %spherical harmonic norm coeffs for m=0
            %boundary conditions
            if(strcmp(params.BC,'periodic'))
                %note - I swapped i_up and i_down so that now i=1 is the topmost cell
                params.i_up = mod((1:params.n_cells)-2,params.n_cells)+1;     %down with periodic BC (try this: mod((1:10)-2,10)+1  )
                params.i_down = mod(1:params.n_cells,params.n_cells)+1;       %up with periodic BC (try this: mod(1:10,10)+1  )
            elseif(strcmp(params.BC,'sharp'))                                 %here the index 0 will stand for "no neighbor"
                params.i_up = (1:params.n_cells)-1;                           %
                params.i_down = mod((1:params.n_cells)+1,params.n_cells+1);   %up with sharp BC (try this: mod((1:10)+1,10+1)  )
            end
            params.i_up = repmat(params.i_up,params.l_max,1);                 %extend to the full matrix form
            params.i_down = repmat(params.i_down,params.l_max,1);
            
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
        
%         % solve steady state
%         function [obj,msg] = solve_ss(obj)          
%         % solve_ss  solve the SS system, store results in obj
%             [obj.y_ss,~,exitflag,output] = fsolve(@(y)obj.dy(0,y),obj.get_ic());
%             msg = sprintf('Exit Flag = %g\n\nMessage:\n%s',exitflag,output.message);
%         end
%         
%         % solve ODE
%         function [obj,msg] = solve_ode(obj)           
%         % solve_ode  solve the ODE, store results in obj
%             [obj.t,obj.y_de] = ode45(@(t,y)obj.dy(t,y),obj.tspan,obj.get_ic());  %solve differential equations
%             msg = 'ODE solution succesful\n';
%             [msgstr, msgid] = lastwarn;
%             lastwarn([msgstr,'_'],[msgid,'_']);         % update lastwarn so we cannot read it twice
%             if strncmpi(msgid, 'MATLAB:ode',10)
%                 msg = sprintf('%s:%s',msgid,msgstr);
%             end
%         end       
        
        % get time derivatives vector
        function dy = dy(obj,~,y)
            %get N and D in matrix format
            [N,D] = obj.y2nd(y);
            
            %Now we calculate the concentrations of N and D on top (?=0) and bottom (?=?)
            %We utilise the face that when m=0 legendre polynomials for x=1 are always 1 and for
            %x=-1 are (-1)^l. Can be verified with the following line:
            %for l=0:30,a=legendre(l,1);b=legendre(l,-1);fprintf('l=%d: P(x=%d) = %d,P(x=%d) = %d\n',l,1,a(1),-1,b(1)),end
            N_top = sum(N .* obj.sph_coeff,1);                           %top is theta=0
            N_bottom = sum(N .* obj.sph_coeff .* (-1).^obj.l_mat,1);  %bottom is theta=pi
            D_top = sum(D .* obj.sph_coeff,1);
            D_bottom = sum(D .* obj.sph_coeff .* (-1).^obj.l_mat,1);
            
            
            
            %calculate boundary conditions
            % BC_N = N_top(i) * D_bottom(i_up) * (i_up~=0) + (-1)^l*N_bottom(i) * D_top(obj.i_down) * (obj.i_down~=0)
            BC_N =  repmat(N_top,obj.l_max,1).*     repmat(circshift(D_bottom,[0 1]),obj.l_max,1) .*    (obj.i_up~=0)  + ...
                repmat(N_bottom,obj.l_max,1).*  repmat(circshift(D_top,[0 -1]),  obj.l_max,1) .*    (obj.i_down~=0) .* obj.motl;
            % BC_D = D_top(i) * N_bottom(obj.i_up) * (obj.i_up~=0) + (-1)^l*D_bottom(i) * N_top(obj.i_down) * (obj.i_down~=0)
            BC_D =  repmat(D_top,obj.l_max,1).*     repmat(circshift(N_bottom,[0 1]),obj.l_max,1) .*    (obj.i_up~=0)  + ...
                repmat(D_bottom,obj.l_max,1).*  repmat(circshift(N_top,[0 -1]),obj.l_max,1) .*      (obj.i_down~=0) .* obj.motl;
            % calculate differentials
            dN =  obj.delta_l - N./obj.N0 - obj.DiffN .* obj.l_mat .* (obj.l_mat+1) .* N - obj.k .* obj.sqrt2lplus1 .* BC_N;
            dD = (obj.delta_l - D./obj.D0 - obj.DiffD .* obj.l_mat .* (obj.l_mat+1) .* D - obj.k .* obj.sqrt2lplus1 .* BC_D) ./ obj.alpha;
            %return equations in vector form
            dy = obj.nd2y(dN,dD);
        end
        
        % get initial conditions vector
        function y_ic = get_ic(obj)
            % get_ic  return y vector of initial conditions.
            % options - 'uniform' or 'zero' or 'biased'
            N = zeros(obj.l_max,obj.n_cells);
            D = zeros(obj.l_max,obj.n_cells);
            switch(obj.IC)
                case 'zero'
                case 'uniform'%setting the first harmonic in each cell to the steady state value (when k=D=0)
                    N = obj.N0 .* obj.delta_l;
                    D = obj.D0 .* obj.delta_l;
                case 'one_higher'%the first cell has 10% more N and D than the steady state, the rest are uniform
                    N = obj.N0 .* obj.delta_l;
                    D = obj.D0 .* obj.delta_l;
                    N(1,1) = N(1,1) * 1.1;
                    N(1,1) = N(1,1) * 1.1;
                case 'biased'%setting the l=1 harmonic to be nonzero so we get bias
                    N = obj.N0 .* obj.delta_l;
                    N(2,:) = N(1,:) * 0.5;
                    D = obj.D0 .* obj.delta_l;
                    D(2,:) = D(1,:) * (-0.5);
                otherwise
                    Exception('invalid mode')
            end
            y_ic = obj.nd2y(N,D);
        end
        
        % convert vector to N,D matrices
        function [N,D] = y2nd(obj,y)
            %converts the long y vector of differential equations to Nothc and Delta
            %matrices, so that N(l,i) is the (l-1)-th spherical harmonic in the i-th cell
            
            %make sure the dimensions are correct - two equations for each l in each
            %cell - one for Notch and one for Delta
            % assert(length(y) == n_cells * l_max * 2)
            N_all = y(1:obj.n_cells * obj.l_max);
            D_all = y(obj.n_cells * obj.l_max + 1 : 2* obj.n_cells * obj.l_max);
            N = reshape(N_all,obj.l_max,obj.n_cells);
            D = reshape(D_all,obj.l_max,obj.n_cells);
        end
        
        % convert N,D matrices to vector
        function y = nd2y(~,N,D)
            % nd2y  convert N,D matrices to long vector format
            y = [N(:);D(:)];
        end
        
        % plot angular distributions
        function ax = plot_distributions(obj)
        % plot_distributions  TODO DOC
            ss_exists = ~isempty(obj.y_ss); %steady state solution exists
            de_exists = ~isempty(obj.y_de); %ODE solution exists
            if(~ss_exists && ~de_exists)
               throw(MException('DLI_1D_CHAIN:plot_distributions','Attempt to plot before getting results')); 
            end
            DLI.fullscreen_figure();
            plot_times = ceil([1,find(obj.t>0.1*obj.t(end),1),find(obj.t>0.5*obj.t(end),1),length(obj.t)]);%times where we will plot the distributions
       
            %get steady state
            if(ss_exists)
                [N_dist,D_dist] = obj.ND_distribution(obj.y_ss);
                Nss = N_dist(1:obj.cells_to_plot,:)';Nss=Nss(:);
                Dss = D_dist(1:obj.cells_to_plot,:)';Dss=Dss(:);
            end
            
            %plot distributions
            counter=0;
            N_times = length(plot_times);
            theta_vec = repmat(linspace(0,pi,obj.N_theta),1,obj.cells_to_plot) + ...%vector of theta values
                sort(repmat(0:obj.cells_to_plot-1,1,obj.N_theta))*pi;%this part offsets the angles of the cells
            theta_vec = theta_vec/pi*180;%rad to deg
            ax = zeros(1,length(plot_times));
            for it = plot_times
                counter = counter+1;
                ax(counter) = subplot(N_times,1,counter);
                box on,hold on
                if(ss_exists)
                    hN = plot(theta_vec,Nss,    '--b',  'LineWidth',2,'DisplayName','Notch steady state');hN.Color(4)=0.2;
                    hN = plot(theta_vec,Dss,    '--r',  'LineWidth',2,'DisplayName','Delta steady state');hN.Color(4)=0.2; 
                end
                if(de_exists)
                    [N_dist,D_dist] = obj.ND_distribution(obj.y_de(it,:));
                    N_vec = N_dist(1:obj.cells_to_plot,:)';
                    D_vec = D_dist(1:obj.cells_to_plot,:)';
                    hN = plot(theta_vec,N_vec(:),'b',   'LineWidth',2,'DisplayName','Notch'             );hN.Color(4)=0.6;
                    hD = plot(theta_vec,D_vec(:),'r',   'LineWidth',2,'DisplayName','Delta'             );hD.Color(4)=0.6;
                end
                xlim(ax(counter),[0,obj.cells_to_plot*180])
                set(ax(counter),'XTick',0:45:theta_vec(end))
                set(ax(counter),'XTickLabel',{''})
                ylabel(sprintf('Level [a.u]\nt=%.1f',obj.t(it)),'FontSize',24)
                set(gca,'FontSize',14)
            end
            x_tick_labels = repmat({'180|0','45','90','135'},1,obj.cells_to_plot);x_tick_labels{1}='0';x_tick_labels{end}='180';
            set(ax(end),'XTickLabel',x_tick_labels)
            xlabel(ax(end),'\theta[deg]','FontSize',24)
            linkaxes(ax)
            
            %plot cell boundaries
            y_lims = get(gca,'YLim');
            for ax0=ax
                for j = 1:obj.cells_to_plot
                    plot(ax0,[1,1]*180*j,y_lims,'--k')
                end
            end
            
            %adjust subplot spacing
            p1=get(ax(1),'Pos');p2=get(ax(2),'Pos');p3=get(ax(end),'Pos');
            dp = abs(p2(2)-p1(2));      %get vertical space between plots
            dh = p3(2)-0.08;
            for ax0=ax %set vertical spacing to that value with 0.02 spacing
                p0=get(ax0,'Pos');
                p0(4)=dp-0.02;
                p0(2)=p0(2)-dh;
                set(ax0,'Pos',p0);
            end
            drawnow
            
        end
        
        % plot total number of molecules
        function ax = plot_total(obj)
        % plot_total TODO DOC
            ss_exists = ~isempty(obj.y_ss); %steady state solution exists
            de_exists = ~isempty(obj.y_de); %ODE solution exists
            if(~ss_exists && ~de_exists)
                throw(MException('DLI_1D_CHAIN:plot_total','Attempt to plot before getting results'));
            end
            % get totals for ODE solution
            if(de_exists)
                t_ind = round(linspace(1,length(obj.t),100));%we only work with 100 time samples in this function
                t_plot = obj.t(t_ind);
                Ntot = zeros(length(t_plot),obj.n_cells);
                Dtot = zeros(length(t_plot),obj.n_cells);
                theta_vec = repmat(linspace(0,pi,obj.N_theta),obj.n_cells,1);
                dtheta = pi/obj.N_theta;
                integral_func=@(X_dist)2*pi * sum(X_dist.*sin(theta_vec),2)*dtheta;
                for j = 1:length(t_plot)
                    [N_dist,D_dist] = obj.ND_distribution(obj.y_de(t_ind(j),:)) ;
                    Ntot(j,:) = integral_func(N_dist);%2*pi * sum(N_dist(1:cells_to_plot,:).*repmat(sin(theta_vec) * dtheta,obj.cells_to_plot,1),2);
                    Dtot(j,:) = integral_func(D_dist);%2*pi * sum(D_dist(1:cells_to_plot,:).*repmat(sin(theta_vec) * dtheta,obj.cells_to_plot,1),2);
                end
            end
            % get totals for SS solution
            if(ss_exists)
                [N_dist,D_dist] = obj.ND_distribution(obj.y_ss);
                Nss = integral_func(N_dist);                        %steady state
                Dss = integral_func(D_dist);                        %steady state
                med_N = mean(Nss);                                  %mean level in the steady state
                med_D = mean(Dss);                                  %mean level in the steady state
            end
            % create plot
            hFig=figure;drawnow;warning off;jFig = get(handle(hFig), 'JavaFrame');jFig.setMaximized(true);warning on; %#ok<WNON,WNOFF>
            ax = zeros(1,4);
            ax(1)=subplot(3,2,1);       set(gca,'FontSize',14);hold on;box on
            ax(2)=subplot(3,2,2);       set(gca,'FontSize',14);hold on;box on
            ax(3)=subplot(3,2,[3 4]);   set(gca,'FontSize',14);hold on;box on
            ax(4)=subplot(3,2,[5 6]);   set(gca,'FontSize',14);hold on;box on
            % steady state diff plot
            if(ss_exists)
                width1 = 0.8;
                bar(ax(1),1:obj.n_cells,abs((Nss-med_N)/med_N),width1,'FaceColor','b')
                bar(ax(2),1:obj.n_cells,abs((Dss-med_D)/med_D),width1,'FaceColor','r')
                % axes(ax(2));[~,h1,h2]=plotyy(1:obj.n_cells,(Nss-med_N)/med_N,1:obj.n_cells,(Dss-med_D)/med_D);set(h1,'Marker','x','Color','b');set(h2,'Marker','x','Color','r')
                set(ax(1),'XTick',1:obj.n_cells,'XLim',[0,obj.n_cells+1],'YScale','log')
                set(ax(2),'XTick',1:obj.n_cells,'XLim',[0,obj.n_cells+1],'YScale','log')
                grid(ax(1),'minor')
                xlabel(ax(1),'Cell #','FontSize',24);
                xlabel(ax(2),'Cell #','FontSize',24);
                ylabel(ax(1),'\(log_{10}(|\frac{N_{ss}-N_{med}}{N_{med}}|)\)','FontSize',24,'interpreter','latex')
                ylabel(ax(2),'\(log_{10}(|\frac{D_{ss}-D_{med}}{D_{med}}|)\)','FontSize',24,'interpreter','latex')
            end
            % time plot
            markers = ['x','o','^','+','*','s','v','<','>','p','h'];n_markers=length(markers);
            % hN = plot(ax(2),[t(1),t(end)],Nss,'--b','LineWidth',2,'DisplayName','Notch_S_S');hN.Color(4)=0.2;   %steady state
            % hD = plot(ax(3),t,Dss,'--r','LineWidth',2,'DisplayName','Delta_S_S');hD.Color(4)=0.2;   %steady state
            if(de_exists)
                for j = 1:obj.n_cells   %time
                    hN = plot(ax(3),t_plot,squeeze(Ntot(:,j)),['-b',markers(mod(j-1,n_markers)+1)],'LineWidth',0.2,'DisplayName',['Notch',num2str(j)]);hN.Color(4)=0.6;
                    hD = plot(ax(4),t_plot,squeeze(Dtot(:,j)),['-r',markers(mod(j-1,n_markers)+1)],'LineWidth',0.2,'DisplayName',['Delta',num2str(j)]);hD.Color(4)=0.6;
                end
                %titles
                % xlabel(ax(3),'\tau [a.u.]','FontSize',24)
                xlabel(ax(4),'\tau [a.u.]','FontSize',24)
                ylabel(ax(3),'Notch total [a.u.]','FontSize',24)
                ylabel(ax(4),'Delta total [a.u.]','FontSize',24)
                linkaxes(ax(3:4),'x')
            end
        end
        
        % get angular distributions from coefficients
        function [N_dist,D_dist] = ND_distribution(obj,y)
        % ND_distribution  converts the N and D coefficients to functions of theta
        % [N_dist,D_dist] = obj.ND_distribution(y) return vectors for N
        %  and for D with dimensions corresponding to obj.N_theta.
            [N,D] = obj.y2nd(y);
            ylm_vec = zeros(obj.l_max,1,obj.N_theta);
            for l = 1:obj.l_max
                ylm_vec(l,1,:) = squeeze(DLI.Y_lm(l-1,linspace(0,pi,obj.N_theta),0));
            end
            ylm_mat = repmat(ylm_vec,1,obj.n_cells,1);
            N_mat = repmat(N,1,1,obj.N_theta);%strech the dimensions to cover theta
            D_mat = repmat(D,1,1,obj.N_theta);%strech the dimensions to cover theta
            %after the following line we are left with dimensions [n_cells,N_theta]
            N_dist = squeeze(sum(N_mat .* ylm_mat,1));
            D_dist = squeeze(sum(D_mat .* ylm_mat,1));
        end
        
    end
    
    methods(Static)
        % convert input to vector form
        function y = vectorize(x,l,n)
            % vectorize receives input that is either in scalar, vector or
            % matrix format and return it in matrix format, with dimensions
            % [l,n]
            if(length(x)==1)
                y=x*ones(l,n);
            elseif(size(x,1)==1 && size(x,2)==n)
                y = repmat(x,l,1);
            elseif(size(x,1)==n && size(x,2)==1)
                y = repmat(x',l,1);
            else
                throw(MException('BadInputDimensions'))
            end
            
        end
    end
    
end


classdef (Abstract) DLI
    %DLI Abstract class for a Diffusion Lateral Inhibition model solver
    %   This class contains all the relevant variables and methods to solve
    %   the Diffusion Lateral Inhibition model dynamically and in steady
    %   state, and to plot and save the results. Classes solving for
    %   specific geometries will inherit from this class.
    
    properties
        n_cells     % number of cells in the tissue
        l_max       % l order of highest spherical harmonic
        D0          % normalized Delta steady-state levels
        N0          % normalized Notch steady-state levels
        alpha       % the ration between notch and delta production, ?_D/?_N
        DiffD       % normalized Delta diffusion rate
        DiffN       % normalized Notch diffusion rate
        k           % normalized interaction strength
        t_max       % simulation time
        tspan       % [0 t_max]
        N_theta     % number of point in the theta axis when displaying results
        N_phi       % number of point in the phi axis when displaying results
        BC          % a string representing boundary conditions, can be 'periodic','sharp', or additional definitions in the inheriting class
        IC          % a string representing initial conditions, can be 'uniform', 'zero' ,'biased', 'one_higher' or additional definitions in the inheriting class
        t           % a vector of time samples of the ODE solution
        y_de        % a vector of the results of the ODE solution
        y_ss        % a vector of the results of the steady state solution
        solver      % function handle to ode solver (default is ODE45)
    end
    
    methods (Abstract)
        %         [obj,msg] =         solve_ode(obj)          % solve the ODE system, store results in obj
        %         [obj,msg] =         solve_ss(obj)           % solve the steady state, store results in obj
        dy =                dy(obj,t,y)             % get a vector dy of time derivatives of the model
        y_ic =              get_ic(obj)             % get vector of initial conditions
        ax =                plot_distributions(obj) % plot molecule distributions on the cells
        ax =                plot_total(obj)         % plot total number of molecules in the cells
        %         axarr =             plot_all(obj)           % run all plotting functions
    end
    
    methods
        
        % solve both ODE and SS, and calls all the plotting functions
        function [obj,msg,axarr] = run_and_plot(obj,varargin)
            % run_and_plot  solves both ODE and SS, and calls all the plotting functions
            %   [msg,axarr] = run_and_plot(obj) runs the solution and
            %    returns the output message in msg and an array of axes of
            %    the plots in axarr
            % varargin can contain the string "skip_ss" to avoid steady state solution
            % varargin can contain the string "skip_ode" to avoid ODE solution
            % varargin can contain the string "skip_plot" to avoid plotting the results
            % varargin can contain the string "skip_print" to avoid printing output message
            
            t0 = tic;                       % start timer
            if(~ismember('skip_ode',varargin))
                [obj,msg1] = obj.solve_ode; % solve ODE
            else
                msg1 = 'Running without ODE solution';
            end
            t_ode = toc(t0);                % stop ode solution timer
            t1 = tic;                       % start timer
            if(~ismember('skip_ss',varargin))
                [obj,msg2] = obj.solve_ss;  % solve SS
            else
                msg2 = 'Running without Steady State solution';
            end
            t_ss = toc(t1);                 % stop steady state solution timer
            if(~ismember('skip_plot',varargin))
                axarr = obj.plot_all();     % run all plotting functions
            else
                axarr = [];
            end
            msg = sprintf([msg1,'\n',msg2,'\n\nRun finished succesully\n\n',...
                'Elpased time [sec]:\nSteady state solution:\t%.2f\nODE solution:\t\t\t%.2f\nTotal time elapsed:\t\t%.2f\n'],t_ss,t_ode,toc(t0));
            if(~ismember('skip_print',varargin))
                fprintf(msg)             % print output message
            end
        end
        
        % solve steady state
        function [obj,msg] = solve_ss(obj)
            % solve_ss  solve the SS system, store results in obj
            [obj.y_ss,~,exitflag,output] = fsolve(@(y)obj.dy(0,y),obj.get_ic());
            msg = sprintf('Exit Flag = %g\n\nMessage:\n%s',exitflag,output.message);
        end
        
        % solve ODE
        function [obj,msg] = solve_ode(obj)
        % solve_ode  solve the ODE, store results in obj
            opts = odeset('RelTol',1e-6);
            [obj.t,obj.y_de] = obj.solver(@(t,y)obj.dy(t,y),obj.tspan,obj.get_ic(),opts);  %solve differential equations
            msg = 'ODE solution succesful\n';
            [msgstr, msgid] = lastwarn;
            if strncmpi(msgid, 'MATLAB:ode',10)
                lastwarn([msgstr,'_'],[msgid,'_']);         % update lastwarn so we cannot read it twice
                msg = sprintf('%s:%s',msgid,msgstr);
            end
        end
        
        % run all plotting functions
        function ax = plot_all(obj)
            % plot_all runs all plotting functions
            ax1 = obj.plot_distributions();
            ax2 = obj.plot_total();
            ax = [ax1,ax2];
        end
        
    end
    
    methods(Static)
        
        function ylm = Y_lm(l,theta,phi,varargin)
            % Spherical harmonic function - returns answer as a 3d mat with dimensions
            % [m,theta,phi]. If varargin is not empty, return results for all m. Otherwise,
            % only for m=0.
            %   ylm = Y_lm(l,theta,phi) returns a 2D matrix of values
            %    corresponding to theta and phi vector of the spherical harmonic l,0
            %
            %   ylm = Y_lm(l,theta,phi,[]) returns a 3D matrix of values
            %    corresponding to theta and phi vector of the spherical harmonic
            %    l,m, and the third dimension corresponds to m values (-l:l)
            
            %make sure we have theat as a column vector and phi a row vector so the
            %outer product works
            if(~iscolumn(theta)),   theta = theta'; end
            if(iscolumn(phi)),      phi = phi';     end
            
            norm = sqrt((2*l+1)/(4*pi));
            p = legendre(l,cos(theta));                                 %returns answers for m=0,1,..,l
            if(isempty(varargin))                                       %m = 0
                p = p(1,:)';
                ylm = norm * p * ones(size(phi));                       %return answer for m=0 only
            else                                                        %all m values
                m = (0:l)';
                factorial_fun=@(mm)sqrt(1./prod(l-abs(mm)+1:l+abs(mm)));%avoid calculating the parts of the factorial that cancel out
                factorial_norm = arrayfun(factorial_fun,m);             %this is the sqrt coefficient
                m_theta_part =  bsxfun(@times,factorial_norm,p);        %this is the part that depends on m and theta
                m_phi_part =    exp(1i*m*phi);                          %this is the part that depends on m and phi
                ylm = norm * bsxfun(@times,m_theta_part,...             %in this line we shift the dimension and multiply
                    permute(m_phi_part,[1,3,2]));                       %   so that the result is [m,theta,phi]
                ylm = bsxfun(@times,[conj(ylm(end:-1:2,:,:));ylm],...   %in this line we expand the results to negative m values
                    (sign(-l:l).^abs(-l:l))');                          %   we for negative m we take CC and multiply by (-1)^|m|
            end
        end
        
        % return a full screen figure
        function hFig = fullscreen_figure()
            hFig=figure;drawnow;warning off;jFig = get(handle(hFig), 'JavaFrame');jFig.setMaximized(true);warning on;%#ok<WNON,WNOFF>
        end
    end
    
end


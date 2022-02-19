function [f,sigma_f,order,condition,BIC] = trend_estimation(u,y,width,order,Nlevels,verbose)

% [f,sigma_f,order,condition,BIC] 
%           = solar_trend_estimation(u,y,width,order,Nlevels,order_max,verbose)
%
% Trend correction by comparing pairs of events with the same value in time series u
% and using these to detect and correct trends in time series y. Both time
% series are expected to have daily values. The corrected time series 
% is  y_corr = y./f;
%
% u        Reference signal (array with daily values)
% y        Signal to be corrected (array with daily values)
% width    Width of gaussian lowpass filter in days, required to remove solar
%          rotational variability. Default (0) is width = 81 days
% order    Order of the model to be adjusted to the data. 
%          If order<-1 the model scans orders 0 to order_max and
%                    selects the best one automatically
%          If order=-1 (default) the model scans orders 0 to order_max and 
%                    calls for a manual selection of order, showing the condition number
%          If order>=0 the model uses the specified order.
% Nlevels  number of amplitude levels used, typically 100-200. 
%          By default Nlevels = 150
% verbose  Display no results (verbose=0), main results only (verbose=1), and
%          all details (verbose=2). Default is verbose=1
%
% f        Response function of y (vertical array)
% sigma_f  One standard deviation of f, obtained by standard bootstrapping
%          (vertical array)
% order    Order used by the model
% condition  Condition number for the chosen order
% BIC      Array of values of the BIC (provided only if order<0)
%
% version 1.0   T. Dudok de Wit  Feb. 2022     ddwit@cnrs-orleans.fr
% for supporting information, see T. Dudok de Wit, "Detecting undocumented trends in 
% solar irradiance observations", J Space Weather Space Climate, 2022



% default values
if nargin<6,    verbose = 1;    end
if nargin<5,    Nlevels = 150;  end
if nargin<4,    order = -1;     end
if nargin<3,    width = 81;     end

if Nlevels<=0,  Nlevels = 150;  end
if width<=0,    width   = 81;   end


shortest_separation = 365*2;      % intervals<shorter_separation (in days) are ignored
Nboot = 100;                % number of evaluations for bootstrapping 
width_median = 27;          % width [days] of median filter for removing solar rotation
condition_max = 1000;       % maximum allowed condition number
order_max = 50;             % largest order when running the model automatically 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sanity checks

if size(y,1)>1 && size(y,2)>1, error('** y must be an array **'); end
if size(u,1)>1 && size(u,2)>1, error('** u must be an array **'); end
y = y(:);
u = u(:);
Nsamples = length(u);                      % number of samples in each record
if Nsamples~=length(y), error('** u and y must have same length **'); end

if Nlevels<=20, error('** ERROR : Number of levels should be > 20 **'); end 
if Nlevels>600, disp('** BEWARE : Number of levels is likely to be too large **'); end 
if width<30, disp('** BEWARE : use larger width to exclude rotational variability **'); end 

% time is the time array in days, with time(1) = 0
time_days = (0:Nsamples-1)';                    % absolute time in days
time_reduced = linspace(0,1,Nsamples)';         % 0 < time_reduced < 1

% keep copy of original signal
u_original = u;
y_original = y;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove solar rotational variability. First apply running median filter, followed
% by lowpass filter
u = movmedian(u,width_median,'endpoints','shrink');  % running median filter
y = movmedian(y,width_median,'endpoints','shrink');
u = smoothdata(u,'gaussian',width*2);                % lowpass filter
y = smoothdata(y,'gaussian',width*2);


% bin the reference signal into Nlevel equidistant levels
u_level = linspace(min(u), max(u), Nlevels+1)';
delta_u = u_level(2)-u_level(1);           % vertical bin size for u
delta_y = (max(y)-min(y))/Nlevels;         % vertical bin size for y

w = 1:Nsamples-1;
index = [];         % contains the first and last index of each crossing
u_value = [];       % contains the corresponding values of u
y_value = [];       % contains the corresponding values of y
slope   = [];
Tsamp = 1;          % sampling period = 1 day


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% determine and store the number of crossings

% plot progress with pairs detected for each level if verbose > 1
if verbose>1
    fprintf('\n')
    disp('Estimating pairs of values with the same level')
    clf        
    subplot(211)
    h1 = plot(time_days,u,time_days(1),u(1),'*');
    ylabel('reference signal r(t)')
    xlabel('time [days]')

    subplot(212)
    h2 = plot(time_days,y,time_days(1),y(1),'*');
    ylabel('test signal s(t)')
    xlabel('time [days]')
    drawnow
end

for i1=1:Nlevels
    a = u_level(i1);
    
    % find indices {k} for which u(k) crosses amplitude level a
    k = find((u(w)>=a & u(w+1)<a) | (u(w)<=a & u(w+1)>a));
    
    % use linear interpolation to refine the location of the crossing
    k_exact = k + (a-u(k))./(u(k+1)-u(k)); 
    
    % find corresponding value in y by linear interpolation 
    y_slope = y(k+1)-y(k);
    y_exact = y(k) + y_slope./(u(k+1)-u(k)).*(a-u(k));    
    
    Nk_exact = length(k_exact);       % number of level crossings (need at least two)
    
    % collect all the results in tables
    if Nk_exact>1
        for j1=1:Nk_exact
            for j2=j1+1:Nk_exact
                if k_exact(j2)-k_exact(j1) > shortest_separation;
                    index   = [index; k_exact(j1) k_exact(j2)];     % timing
                    u_value = [u_value; a a];                       % crossing of u
                    y_value = [y_value; y_exact(j1) y_exact(j2)];   % crossing of y
                    slope   = [slope; y_slope(j1) y_slope(j2)];
                end
            end
        end
    end
    

    % show detailed progress if necessary, level by level
    if verbose>1
        if i1==1, 
            clf
            subplot(211)
            h1 = plot(time_days,u,k_exact,a*ones(Nk_exact,1),'*-');
            ylabel('reference signal r(t)')
            xlabel('time [days]')

            subplot(212)
            h2 = plot(time_days,y,k_exact,y_exact,'*',k_exact,y_exact,'-');
            ylabel('test signal s(t)')
            xlabel('time [days]')
            drawnow
        else
            h1(2).XData =  k_exact;
            h1(2).YData =  a*ones(Nk_exact,1);
            h2(2).XData =  k_exact;
            h2(2).YData =  y_exact;
            h2(3).XData =  k_exact;
            h2(3).YData =  y_exact;
            drawnow
        end
    end
end

% check whether there are enough pairs
Npairs = size(index,1);         % number of detected pairs
if Npairs<10, warning('** warning **    number of crossings is too small '); end


% the effective nr of pairs is actually smaller because of serial correlation 
s = sort(index(:,1));
s(diff(s)<0.01) = [];
step = median(diff(s));
correction = width/step;
if correction<1, correction = 1; end

% estimate the effective nr of pairs, compounded by serial correlation
Npairs_eff = Npairs/correction;     


if verbose>1
    disp(['  Number of detected pairs  : ',int2str(Npairs)])
    disp(['  Effective number of pairs : ',int2str(Npairs_eff)])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the regression matrix, the weight is inversely proportional to the
% local slope ==> flat portions of data put a stronger constraint

offset = delta_y*delta_y / 10000;
weight = 1.0./((slope(:,1).*slope(:,1) + slope(:,2).*slope(:,2))*(Tsamp*Tsamp) + offset);

ampl = mean(y);
y1 = y_value(:,1)/ampl;
y2 = y_value(:,2)/ampl;

t1 = (index(:,1)-1)/(Nsamples-1);   % first crossing in reduced time
t2 = (index(:,2)-1)/(Nsamples-1);   % second crossing in reduced time



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% if order<0 determine order manually by scanning orders from 0 to order_max

Npairs_valid = ceil(0.9*Npairs);      % cross-validation: use 90% of values for testing 

if order<0
    if verbose>1
        fprintf('\n')
        disp('Estimating and testing model for different orders')
    end
    order_test = (0:order_max)';    % array with the numbers of orders
    Norder = length(order_test);
    
    m = y1-y2;
    condition = zeros(Norder,1);    % condition number
    BIC = zeros(Norder,1);          % Bayesian Information Criterion
    logNpairs = log(Npairs_eff);
    m_variance = var(m);
    
    i1 = 1;
    while i1<=Norder
        n = order_test(i1)+1;
        
        % build regression matix
        M = zeros(Npairs,n);
        M(:,1) = y2.*t1 - y1.*t2;
        for i=1:order_test(i1)
            M(:,i+1) = y2.*sin(i*pi*t1) - y1.*sin(i*pi*t2);
        end
        M_weight = (weight*ones(1,n)).*M;
        m_weight = weight.*m;
        condition(i1) = cond(M_weight);         % condition of linear system
        
        coeff = M_weight\m_weight;
        resid = M*coeff - m;
        RSS = sum(resid'*resid);                % residual sum of squares      
        BIC(i1) = logNpairs*(i1+1) + Npairs_eff*log(RSS/Npairs);
        
        if condition(i1)>condition_max
            Norder = i1;
            BIC = BIC(1:Norder);
            condition = condition(1:Norder);
        end
        i1 = i1+1;
    end
    
    [BIC_min,k_min] = min(BIC);
    
    % display the results
    clf
    subplot(211)
    semilogy(1:Norder,condition,'o-',[1 Norder],[condition_max condition_max],'k--');
    grid on
    xlabel('order')
    ylabel('condition nr')
    
    subplot(212)
    plot(1:Norder,BIC,'.-',k_min,BIC_min,'*')
    grid on
    xlabel('order')
    ylabel('BIC')
    drawnow


    if order<-1
        order = k_min;
        if verbose>0, disp(['  Order = ',int2str(k_min)]); end
    else
        fprintf('\n')
        disp( 'Manually select the order of the trend approximation')
        disp(['  Minimum error occurs for order ',int2str(k_min)])
        order = input('  Select order : ');
        if isempty(order), order = 0; end
    end
else
    BIC = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Solve system to find model coefficients

% build the regression matrix
M_weight = zeros(Npairs,order+1);
M_weight(:,1) = weight.*(y2.*t1-y1.*t2);
for i=1:order
    M_weight(:,i+1) = weight.*(y2.*sin(i*pi*t1)-y1.*sin(i*pi*t2));
end

condition = cond(M_weight);
condition_min = condition;
if verbose
    disp(['  Order : ',int2str(order)]); 
    disp(['  Condition number : ',num2str(condition_min)]); 
end

m_weight = weight.*(y1-y2);
coeff = M_weight\m_weight;    % solve linear system

if verbose>1
    clf
    fit = (M_weight*coeff)./weight;
    R = corrcoef(y1-y2,fit);
    plot(y1-y2,fit,'+')
    title(['Pearson correlation  R = ',num2str(R(1,2),4)])
    xlabel('y_1(t)-y_2(t) observed')
    ylabel('y_1(t)-y_2(t) modelled')
    grid on
    fprintf('\n')
    disp('-- Hit any key to continue -- ')
    pause
end
    
% regress contains the different regressors: linear trend + sine functions
regress = zeros(Nsamples,order+1);
regress(:,1) = time_reduced;
for i=1:order
    regress(:,i+1) = sin(i*pi*time_reduced);
end
f = 1 + regress*coeff;      % the correction 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% use bootstrapping to estimate confidence intervals
%%%%%%%%%%%%%%% this version uses standard bootstrapping, NOT block
%%%%%%%%%%%%%%% bootstrapping, which is better

if verbose>1
    fprintf('\n')
    disp('Estimating uncertainties by standard bootstrapping')
end

fboot = zeros(Nsamples,Nboot);
for i=1:Nboot
    ind = randi(Npairs,Npairs,1);
    coeff_boot = M_weight(ind,:)\m_weight(ind);
    fboot(:,i) = 1 + regress*coeff_boot;
end
sigma_f = std(fboot,0,2);       % standard deviation of f
ave_f = mean(fboot,2);          % should be almost equal to f


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Display the main results

if verbose
    clf
    subplot(211)
    ycorr = y./f;
    plot(time_days,y,time_days,ycorr,'--')
    xlabel('time [days]')
    ylabel('amplitude')
    legend('observed y(t)','corrected y(t)','location','northeast')
    grid
    title(['Order: ',int2str(order),'    # levels: ',int2str(Nlevels), ...
        '    Condition number: ',int2str(condition_min)])
    
    subplot(212)
    h = plot(time_days,[f-sigma_f f+sigma_f],time_days,f,'--');
    set(h(3),'Color',get(h(2),'Color'));
    set(h(2),'Color',get(h(1),'Color'));
    xlabel('time [days]')
    ylabel('f(t)')
    title('Correction f \pm \sigma  :   y_{obs} = y_{corr} * f')
    grid
    
    drawnow
end

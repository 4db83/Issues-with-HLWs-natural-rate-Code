function varargout = plotacf(x,no_acf,ylims,plfg,FNT,plot_dims,Legnd)
% % % call as: 
% % %				[ac,pac] = plotacf(x);  with automatic axis limits and 50 acfs
% % % or as full call:
% % %				[ac,pac] = plotacf(x,No_acf,ylim_lower,plfg,FontSize)
% % % NOTE: ylim_lower is [-.2 -.2] or [-.5 -.5] or so. ylim_upper is always 1. This is just for lower
% % %				plot limit control.
% % %---------------------------------------------------------------------------------------------------
% % % x				= data vector (column)
% % % nac			=	no. acf values to return (<= length(x))
% % % plfg		= 1 do a plot
% % % acalpha	= alpha for CI in ACF/PACF plots
% % 


% check if series is an fints object
if isa(x, 'fints')
	x = fts2mat(x);
end

% remove nans if exist
x = rmvnan(x);
plot_dims0 = [.37 .2];

% set some default values
SetDefaultValue(2,'no_acf',30)
SetDefaultValue(4,'plfg',1)
SetDefaultValue(5,'FNT',[16 16])
SetDefaultValue(6,'plot_dims',[])
SetDefaultValue(7,'Legnd',1)

tickshrinkValue = 0.7;

FNT1 = FNT(1);
if length(FNT)==1
	FNT2 = FNT(1);
else
	FNT2 = FNT(2);
end

Ipd = isempty(plot_dims);

if Ipd
 	plot_dims = plot_dims0;
end

%%
% % function main( )
% % 
% % clear all; clc;
% % 
% % seed(12300111);
% % b1 = -3.5 ; b2 = -2;
% % bL = [1 +b1 +b2];
% % aL = [1 -.9 +0.2];
% % u = randn(2000,1);
% % x = filter(bL,aL,u);
% % 
% % 
% % no_acf	= 50;
% % plfg		= 1;
% % 

%%
acalpha	=.05;

% compue the ACFs and PACFs
[nx nr] = size(x);
ac0			= acf_f(x,no_acf+1);
ac			= ac0(2:end);
pac			= pacf_f(x,no_acf,ac0)';		

% some old commands were
% pac=pacf(x,npac,0); % 0 for no plot by pacf, compute npac values
% ac  = autocorr(x,npac,0);
% pac		= parcorr(x,no_acf,0);
% pac		= pac(2:end);

% lower bounds on the ACF/PACF plot
thr				=	2/sqrt(nx); % 2*sqrt(T)
FF				= 10;
mmac_pac	= min(min(ac),min(pac));
%plt_min		= min( [-thr -.2 min(0,floor(FF*mmac_pac)/FF)] );

xm_grid = ([-.2:-.2:-1])';
plt_min = xm_grid(sum(mmac_pac<xm_grid)+1);

%plt_min		= min( [-thr -.2 min(0,floor(FF*mmac_pac)/FF)] );

SetDefaultValue(3,'ylims',plt_min*ones(2,1));

if isempty(ylims)
	ylims	= plt_min*ones(2,1);
end

if length(ylims)==1
	ylims = ylims*ones(2,1);
end

% FF = 10;
% min_pac = floor(FF*min(pac))/FF;
% min_ac	= floor(FF*min(ac))/FF;

%% plotting.
%CLR = [.73 .83 .96];
CLR = [.7 .8 1];
%CLR = [.4 .7 1];
S = zeros(2,1);

if plfg == 1
% the acf
S(1) = subplot(121);
LL(1)= bar(ac,'FaceColor',CLR,'EdgeColor',CLR);	        % a stem plot with dot marker
BS_H= get(LL(1),'BaseLine');set(BS_H,'LineStyle','-','Linewidth',0.5)
hold on;
    xlim([.25 no_acf+0.5]);
		ylim([ylims(1) 1.001]);
LL(2)= hline_f(thr,'-.r');
		hline_f(-thr,'-.r')
%    hline_f(0,'k')
%		ylabel('ACF','FontSize',FNT);
if Legnd == 1
		L1 = legend(LL,{'ACF','95% CI'},'FontSize',FNT2,'Location','Northeast','Orientation','Horizontal');
		%set(L2,'units','normalized')
% 		L1_tmp = get(L1,'Position');
% 		set(L1,'Position',[.295 .4775 L1_tmp(3:4)]);
% 		get(L1,'Position');
% 		legend('boxoff');
end
%		xlabel('horizon','FontSize',FNT);
		setplot([.07 .33 plot_dims],FNT1);
% 		setplot([plot_dims(1:2)],FNT1);
%		moveylabel(7)
		setxticklabels([1 get(S(1),'XTick')],S(1));
		xlim([.0 no_acf+1]);
		setyticklabels(ylims(1):.2:1);
		setytick(1,FNT1);
		set(gca, 'Layer','top')
		xclr = get(gca,'XColor');
		h0 = hline(0);
		set(h0,'Color',xclr)
		tickshrink(tickshrinkValue)
		hold off;
    
% the pacf
S(2) = subplot(122);
LL(1)= bar(pac,'FaceColor',CLR,'EdgeColor',CLR);	        % a stem plot with dot marker
BS_H= get(LL(1),'BaseLine');set(BS_H,'LineStyle','-','Linewidth',0.5)
    xlim([.25 no_acf+0.5]);
		ylim([ylims(2) 1.001]);
    hold on;
LL(2)= hline_f(thr,'-.r');
		hline_f(-thr,'-.r')
		hline(0,'k')
%		ylabel('PACF','FontSize',FNT);
if Legnd == 1
		L2 = legend(LL,{'PACF','95% CI'},'FontSize',FNT2,'Location','Northeast','Orientation','Horizontal');
% 		set(L2,'units','normalized')
%  		L2_tmp = get(L2,'Position');
%  		set(L2,'Position',[.725 .4775 L2_tmp(3:4)]);
% % 		get(L2,'Position');
% 		legend('boxon');
end
%		xlabel('horizon','FontSize',FNT);
		setplot([.57 .33 plot_dims],FNT2);
% 		setplot([plot_dims(1:2)],FNT1);
%		moveylabel(7)
%		setyticklabels([-.4:.2:1])
		setxticklabels([1 get(S(2),'XTick')],S(2));
		xlim([.0 no_acf+1]);
		setyticklabels(ylims(2):.2:1);
		setytick(1,FNT1);
		set(gca, 'Layer','top')
		set(gca, 'Layer','top')
		xclr = get(gca,'XColor');
		h0 = hline(0);
		set(h0,'Color',xclr)
		tickshrink(tickshrinkValue)
    hold off;
end
%%
% sequeeze together a bit
if Ipd
	subplotsqueeze(gcf, [1 1 .94]);
else
	subplotsqueeze(gcf, [1 1 plot_dims(3)]);
end

if plfg == 2
% the acf
S(1) = subplot(211);
LL(1)= bar(ac,'FaceColor',CLR,'EdgeColor',CLR);	        % a stem plot with dot marker
		baseline_handle = get(LL(1),'BaseLine');
		set(baseline_handle,'LineStyle','None')
    xlim([.25 no_acf+0.5]);
		ylim([ylims(1) 1.001]);
		hold on;		
%		plot(thr,'-.r');
LL(2) = hline_f(thr,'-.r');
		hline_f(-thr,'-.r')
    hline(0.0,'k')
%		ylabel('ACF','FontSize',FNT);
if Legnd == 1
		L1 = legend(LL,{'ACF','95% CI'},'FontSize',FNT2,'Location','Best','Orientation','Horizontal');
		%set(L2,'units','normalized')
		L1_tmp = get(L1,'Position');
		set(L1,'Position',[.725 .827 L1_tmp(3:4)]);
		get(L1,'Position');
		legend('boxoff');
end
%		xlabel('horizon','FontSize',FNT);
		setplot([ 0.10 .55 .8 .3],FNT1);
%		moveylabel(2)
		setxticklabels([1 get(S(1),'XTick')],S(1));
		xlim([0 no_acf+1]);
		setyticklabels(ylims(1):.2:1);
		setytick(1,FNT1);
    hold off;
    
% the pacf
S(2) = subplot(212);
LL(1)= bar(pac,'FaceColor',CLR,'EdgeColor',CLR);	        % a stem plot with dot marker
    baseline_handle = get(LL(1),'BaseLine');
		set(baseline_handle,'LineStyle','None')
		xlim([.25 no_acf+0.5]);
		ylim([ylims(2) 1.001]);
    hold on;
%		plot(thr,'-.r');
LL(2)= hline_f(thr,'-.r');
		hline_f(-thr,'-.r')
		hline(0.0,'k')
		%hl1 = hline(0,'g');set(hl1,'linewidth',.15);
%    ylabel('PACF','FontSize',FNT);
if Legnd == 1
		L1 = legend(LL,{'PACF','95% CI'},'FontSize',FNT2,'Location','Best','Orientation','Horizontal');
		%set(L2,'units','normalized')
		L1_tmp = get(L1,'Position');
		set(L1,'Position',[.725 .457 L1_tmp(3:4)]);
		get(L1,'Position');
		legend('boxoff');
end
%		xlabel('horizon','FontSize',FNT);
		setplot([ 0.10 .18 .8 .3],FNT1);
%		moveylabel(2)
%		setyticklabels([-.4:.2:1])
		setxticklabels([1 get(S(2),'XTick')],S(2));
		xlim([0 no_acf+1]);
		setyticklabels(ylims(2):.2:1);
		setytick(1,FNT1);
    hold off;
end


% to modifiy specific handles
% % icrm = .2;
% % mtcks = [min_pac:icrm:1]'
% % setyticklabels(mtcks,S(2))
% % setytick(4,S(2))
% % setylabel(5,S(2))
% % setxlabel([10.05 .05],S(2))

% print2pdf('here')


if nargout > 0
	varargout{1} = ac;
	varargout{2} = pac;
	varargout{3} = S;
end	


end


function r=acf_f(x,nac)
		% ==========================================================================
		% SUPER FAST ACF COMPUTATION
		% r=acf(x,normflg)
		% ==========================================================================
		%  This function computes the acf r(k) k=0:length(x)
		%  CAUTION the first value r(1) is hat{gamma}(0) IE the first lag is 0
		% x = time series vector (column)
		% normflg 0 to divide by nr=length(x)
		%         1 to divide by nr*sample variance
		%         2 to divide by nr-k
		%         3 to divide by (nr-k)*sample variance
		% ==========================================================================
		[nr,nc]=size(x);
		if nc > nr
		    error('x must be a column vector');
		end
		
		x = x - repmat(mean(x),nr,1);
		
		r = conv(flipud(x),x);
		r = r(nr:end);
		r = r/nr;
		r = r/r(1);
		r = r(1:nac);
end

function pr=pacf_f(x,m,rho)
% ==============================================================================
% pr=pacf(x,m,plfg)
% ==============================================================================
% computes the pacf of the series in x
% ==============================================================================
% x = input series
% m = maximum lag to compute
% pa = pacf(j) j=1:m
% ==============================================================================

		% first compute all the sample correlations needed
		nx=length(x);
		
		
		pr(1)=rho(2); % lag 1 rho
		%pr=zeros(1,m);
		for k=2:m
		   pmat=toeplitz(rho(1:k));% lag 0 to k-1
		   rhovec=rho(2:k+1);		% lag 1 to k
		   phi=inv(pmat)*rhovec;
		   pr(k)=phi(k);
		end
end

function hhh=hline_f(y,in1,in2)
% ==============================================================================
% function h=hline(y, linetype, label)
% ==============================================================================
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001
% ==============================================================================

		if length(y)>1  % vector input
		    for I=1:length(y)
		        switch nargin
		        case 1
		            linetype='r:';
		            label='';
		        case 2
		            if ~iscell(in1)
		                in1={in1};
		            end
		            if I>length(in1)
		                linetype=in1{end};
		            else
		                linetype=in1{I};
		            end
		            label='';
		        case 3
		            if ~iscell(in1)
		                in1={in1};
		            end
		            if ~iscell(in2)
		                in2={in2};
		            end
		            if I>length(in1)
		                linetype=in1{end};
		            else
		                linetype=in1{I};
		            end
		            if I>length(in2)
		                label=in2{end};
		            else
		                label=in2{I};
		            end
		        end
		        h(I)=hline(y(I),linetype,label);
		    end
		else
		    switch nargin
		    case 1
		        linetype='r:';
		        label='';
		    case 2
		        linetype=in1;
		        label='';
		    case 3
		        linetype=in1;
		        label=in2;
		    end
		
		    
		    
		    
		    g=ishold(gca);
		    hold on
		
		    x=get(gca,'xlim');
		    h=plot(x,[y y],linetype,'Linewidth',0.5);
		    if ~isempty(label)
		        yy=get(gca,'ylim');
		        yrange=yy(2)-yy(1);
		        yunit=(y-yy(1))/yrange;
		        if yunit<0.2
		            text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
		        else
		            text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
		        end
		    end
		
		    if g==0
		    hold off
		    end
		    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
		end % else
		
		if nargout
		    hhh=h;
		end
		end

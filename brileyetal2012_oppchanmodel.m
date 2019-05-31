function [chan maa ac] = brileyetal2012_oppchanmodel(figs,params)
% Matlab implementation of the computational opponent-channel model
% described in:
%
% Evidence for opponent-process analysis of sound-source location in humans
% Journal of the Association for Research in Otolaryngology
%
% Paul M. Briley (1*), Pádraig T. Kitterick (1), A. Quentin Summerfield (1,2)
% 1 Department of Psychology, University of York
% 2 Hull York Medical School, University of York
% * email: paul.briley@york.ac.uk
%
% [chan maa ac] = brileyetal2012_oppchanmodel(figs,params)
% both input arguments are optional
%
% outputs:
% * chan
% chan contains two structures, one for the left, and one for the right,
% spatial channel (the tuning function for each channel is modeled with a
% cumulative Gaussian). Each structure contains the channel parameters (m -
% the mean of the underlying Gaussian, and thus the point of maximum slope
% of the cumulative Gaussian, and sd - the standard deviation of the
% underlying Gaussian, and thus a measure of the channel's sharpness of 
% tuning with smaller values meaning sharper tuning and larger values meaning
% broader tuning), the angles over which the channel has been evaluated (x) 
% and the channel tuning function (w).
%
% * maa
% maa contains the angles over which the channels have been evaluated (x),
% the model output (f), the rate-of-change of f with respect to x (df), the
% reciprocal of the rate-of-change (rdf), the minimum audible angle for a
% reference azimuth of 0 degrees (maa0), the constant scale factor used to
% convert rdf to MAA in degrees (c) and the predicted minimum audible
% angles (vals).
%
% * ac
% ac contains two structures, one for each auditory cortex.
% Each structure contains the channel weights used for that cortex
% (weights), the compression factor (comp) and the scale and shift factors.
% resp stores the location-shift response curves produced by the model.
% The shift-direction preference values are stored in out_minus_in
% (positive values indicate a larger response to outward-going, than
% inward-going, location shifts, consistent with the opponent-channel model).
%

if nargin<1; figs = 1; end % show figures?
if nargin<2; params = ''; end

%%%% default model parameters %%%%
m = [0 0]; % point of max slope, in degrees (i.e. mean of underlying Gaussian); [left channel, right channel]
sd = [46, 46]; % tuning sharpness (standard deviation of underlying Gaussian)
left_ac_weights = [0.4570, 0.5430]; % must sum to 1
right_ac_weights = [0.5049, 0.4951];
comp = 7; % compression
scale = 80.4103;
shift = 13.4519;
maa0 = 1.51; % minimum audible angle for a reference azimuth of 0 degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(params,'m'); m = params.m; end
if isfield(params,'sd'); sd = params.sd; end
if isfield(params,'left_ac_weights'); left_ac_weights = params.left_ac_weights; end
if isfield(params,'right_ac_weights'); right_ac_weights = params.right_ac_weights; end
if isfield(params,'comp'); comp = params.comp; end
if isfield(params,'scale'); scale = params.scale; end
if isfield(params,'shift'); shift = params.shift; end
if isfield(params,'maa0'); maa0 = params.maa0; end 
    
x = -90:90; % angles to evaluate the channel tuning curves over, in degrees
loc_shifts = -120:30:120; % location shifts
locs = -60:30:60; % absolute locations
labs = {'-60°','-30°','0°','+30°','+60°'}; % legend labels

cols = [231 120 23; 0 146 63; 0 147 221; 219 33 76; 40 22 111]; % plot colors for location-shift response curves and shift-direction preference bars
cols = cols./255;
linstyle = {'-s','-^','-o','-+','-d'}; % plot line styles for the location-shift response curves
chancols = [0 0 0; 166 166 166]; % plot colors for channel tuning curves
chancols = chancols./255;
chanstyle = {'-','--'}; % line styles for channel tuning curves
moutcol = [89; 89; 89];
moutcol = moutcol./255; % plot color for model output

chan{1} = prep_chan('left channel',x,m(1),sd(1),1); % channels
chan{2} = prep_chan('right channel',x,m(2),sd(2),2);
maa = get_maa(chan,maa0); % use channel tuning curves to predict MAAs
ac{1} = prep_ac('left AC',left_ac_weights,comp,scale,shift); % auditory cortices
ac{2} = prep_ac('right AC',right_ac_weights,comp,scale,shift);

for i = 1:length(ac) % for each auditory cortex...
    ac{i}.resp = nan(length(loc_shifts),length(locs)); % resp stores the location-shift response curves in the format: location shift x post-shift location
    for post_loc = 1:length(locs) % for each post-shift location...
        for pre_loc = 1:length(locs) % for each pre-shift location...
            out = 0; % model output
            for ii = 1:length(chan) % for each channel...
                post_resp = chan{ii}.w(x==locs(post_loc)); % channel output for post-shift location
                pre_resp = chan{ii}.w(x==locs(pre_loc)); % channel output for pre-shift location
                d = post_resp - pre_resp;
                if d>0 % a measurable EEG response is assumed to be produced only when the post-shift location produces a larger output from the channel than the pre-shift location
                    out = out + d.*ac{i}.weights(ii); 
                end
            end
            shift = locs(post_loc)-locs(pre_loc); % location shift
            if ac{i}.comp % compression
                out = 1-exp(-comp.*out);
            end
            out = out .* ac{i}.scale;
            out = out + ac{i}.shift;
            ac{i}.resp(loc_shifts==shift,post_loc) = out;
        end
    end
    ac{i} = shift_dir_pref(ac{i},loc_shifts,locs);
end

if figs    
    figure; % channel tuning curves
    curvs = zeros(1,length(chan));
    curvlabs = cell(1,length(chan));
    for i = 1:length(chan) % for each spatial channel...
        curvs(i) = plot(chan{i}.x,chan{i}.w,chanstyle{i},'Color',chancols(i,:),'linewidth',4);
        curvlabs{i} = chan{i}.name;
        hold all
    end
    title('Channel tuning curves');
    legend(curvs,curvlabs,'location','SouthEast');     
    axis([x(1) x(end) 0 1]);
    
    figure; % model output
    plot(maa.x,maa.f,'Color',moutcol,'linewidth',4);
    axis([maa.x(1) maa.x(end) min(maa.f) max(maa.f)]);
    hold all;
    plot([0 0],[min(maa.f) max(maa.f)],'k--');
    xlabel('x');
    ylabel('f(x)');
    title('Model output');    
    
    for i = 1:length(ac) % for each auditory cortex...
        figure; % location-shift response curves
        lins = zeros(1,length(locs));
        for ii = 1:length(locs) % for each post-shift location...
            lins(ii) = plot(loc_shifts,ac{i}.resp(:,ii),linstyle{ii},'Color',cols(ii,:));
            hold all;
        end
        title(sprintf('Location-shift response curves - %s',ac{i}.name));
        h = legend(lins,labs,'location','SouthEast');
        set(get(h,'title'),'String','Post-shift location:');
        axis([loc_shifts(1) loc_shifts(end) 0 100]);
        xlabel('Location shift (°)');
        ylabel('Response magnitude (nAm)');
    end
    
    for i = 1:length(ac) % for each auditory cortex...
        figure; % shift-direction preference bar chart
        bars = zeros(1,length(ac));
        bars(1) = bar(1,ac{i}.out_minus_in(1),'FaceColor',cols(2,:)); hold all;
        bars(2) = bar(2,ac{i}.out_minus_in(2),'FaceColor',cols(4,:));
        title(sprintf('Shift-direction preference - %s',ac{i}.name));
        axis([0.5 2.5 0 20]);
        ylabel('Out minus in (nAm)');        
        xlabel('Post-shift location (°)');
        set(gca,'XTick',[1 2]);
        set(gca,'XTickLabel',{'-30°','+30°'});
    end
end
if figs; fprintf(1,'\n\n'); end

function chan = prep_chan(name,x,m,sd,dir)
switch(dir)
    case 1 % peak to the left
        w = normcdf(x.*-1,m.*-1,sd);
    case 2 % peak to the right
        w = normcdf(x,m,sd);
    otherwise
        error('dir should be 1 (peak to the left) or 2 (peak to the right)');
end
    
chan.name = name;
chan.x = x;
chan.m = m;
chan.sd = sd;
chan.dir = dir;
chan.w = w;

function maa = get_maa(chan,maa0)
x = chan{1}.x;
f = (chan{2}.w-chan{1}.w)./-(chan{2}.w(1)-chan{1}.w(1)); % model output, f(x)
df = diff(f); % rate-of-change of f(x) with respect to x
rdf = 1./df; % reciprocal of df, larger values (i.e. smaller values of df) suggest larger MAAs
c = maa0./rdf(x==0); % constant used to convert the reciprocal of df to MAA in degrees
vals = rdf.*c; % predicted MAA values

maa.x = x;
maa.f = f;
maa.df = df;
maa.rdf = rdf;
maa.maa0 = maa0;
maa.c = c;
maa.vals = vals;

function ac = prep_ac(name,weights,comp,scale,shift)
if sum(weights)~=1; error('weights should total 1'); end
ac.name = name;
ac.weights = weights;
ac.comp = comp;
ac.scale = scale;
ac.shift = shift;

function ac = shift_dir_pref(ac,loc_shifts,locs)
% ac.resp stores the location-shift response curves in the format: location shift x post-shift location 
in = ac.resp(loc_shifts==(+30),locs==(-30));
out = ac.resp(loc_shifts==(-30),locs==(-30));
ac.out_minus_in(1) = out - in; % post-shift location = -30°

in = ac.resp(loc_shifts==(-30),locs==(+30));
out = ac.resp(loc_shifts==(+30),locs==(+30));
ac.out_minus_in(2) = out - in; % post-shift location = +30°

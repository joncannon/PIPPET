function out = PATIPPET_params(varargin)

p = inputParser;
addParameter(p,'n_streams',1,@isnumeric);
addParameter(p,'events_unit',[.25],@isnumeric);
addParameter(p,'variance_unit',[.0001],@isnumeric);
addParameter(p,'lambda_unit',[0.02],@isnumeric);
addParameter(p,'lambda_0',0.01,@isnumeric);
addParameter(p,'expected_cycles',4,@isnumeric);
addParameter(p,'expected_period',.25,@isnumeric);

addParameter(p,'event_times',[1],@isnumeric);
addParameter(p,'tmax',nan,@isnumeric);
addParameter(p,'dt',0.001,@isnumeric);
addParameter(p,'mu_0',[0, 1],@isnumeric);
addParameter(p,'C_0',[.0001,0; 0,.04],@isnumeric);
addParameter(p,'sigma',0.05,@isnumeric);
addParameter(p,'display', true, @isboolean);
addParameter(p,'x_list', [0], @isboolean);

addParameter(p,'sigma_2',0.02,@isnumeric);
addParameter(p,'true_speed',1,@isnumeric);

addParameter(p,'title', '');


parse(p,varargin{:})
out = p.Results;
out.streams = cell(1,p.Results.n_streams);

last_events = [];
last_expected = [];

if out.n_streams > 1 || iscell(out.events_unit)
    for j = 1:out.n_streams
        out.streams{j} = PATIPPET_stream_params(out.events_unit{j}, out.variance_unit{j}, out.lambda_unit{j}, out.lambda_0{j}, out.expected_cycles{j}, out.expected_period{j}, out.event_times{j});
        if ~isempty(out.streams{j}.event_times)
            last_event(end+1) = max(out.streams{j}.event_times);
        end
        last_expected(end+1) = max(out.streams{j}.e_means);
    end
else
    out.streams{1} = PATIPPET_stream_params(out.events_unit, out.variance_unit, out.lambda_unit, out.lambda_0, out.expected_cycles, out.expected_period, out.event_times);
    if ~isempty(out.streams{1}.event_times)
        last_events = out.streams{1}.event_times(end);
    end
    last_expected(end+1) = max(out.streams{1}.e_means);
end

if isnan(out.tmax)
    if ~isempty(last_events)
        out.tmax = max(last_events) + 0.2;
    else
        out.tmax = max(last_expected) + .2;
    end
end

for j = 1:out.n_streams
    out.streams{j}.expect_func = @(x_list) expectation_func(x_list, out.streams{j}.e_means, out.streams{j}.e_vars, out.streams{j}.e_lambdas, out.streams{j}.lambda_0);
end
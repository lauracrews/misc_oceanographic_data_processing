% Downloads the 1-minute bin-averaged Sikuliaq underway data from Coriolix.
%
% Get URLs for desired data and time ranges from https://coriolix.sikuliaq.alaska.edu/data/download/binned/
% Note that data are divided into categories in the API:
% -Flowthrough ("flow")
% -Atmospheric ("met")
% -Navigation ("nav")
% -Other ("other")
%
% Necessary to process each of these data categories separately (they each
% have their own metadata file mapping the parameter numbers to the
% instrument and data type). Copy and past the URLs for each type below
%
% Per https://coriolix.sikuliaq.alaska.edu/data/format/, the 1-minute data
% product includes a few quality control options prior to averaging, and
% returns the statistics from bin averaging.
% 
% Here we use the data quality b: "Best Bins (b) are derived only from data
% that have flag values of 0, 1, or 2 (i.e. no suspect or failed data)".
% 
% Available bin averaging statistics are: std, min, max, num, spotval, and
% median. Here we retain the mean and std
%
% Oct 2025 L. Crews

% Flowthrough variables: Here temperature at the seawater intake and salinity
flow_url = 'https://coriolix.sikuliaq.alaska.edu/erddap/tabledap/binned_default_flow_rolling.csv?time,latitude,longitude,parameter_19,parameter_03,parameter_23&time>=2025-10-17T00:00&time<=2025-11-12T23:59&orderBy("time")';
flow_metadata_url = 'https://coriolix.sikuliaq.alaska.edu/erddap/info/binned_default_flow_rolling/index.csv';
[S_flow, S_flow_info] = get_SKQ_underway_data(flow_url, flow_metadata_url);
S_flow_renamed = rename_underway_variables_to_SWIFT_compliance(S_flow);

% Met variables: Here wind speed and direction, air temperature, pressure,
% relative humidity, SW and LW radiation, PAR
met_url = 'https://coriolix.sikuliaq.alaska.edu/erddap/tabledap/binned_default_met_rolling.csv?time,latitude,longitude,parameter_30,parameter_28,parameter_23,parameter_22,parameter_17,parameter_11,parameter_09,parameter_08,parameter_07,parameter_06,parameter_05,parameter_04,parameter_03&time>=2025-10-17T00:00&time<=2025-11-12T23:59&orderBy("time")';
met_metadata_url = 'https://coriolix.sikuliaq.alaska.edu/erddap/info/binned_default_met_rolling/index.csv';
[S_met, S_met_info] = get_SKQ_underway_data(met_url, met_metadata_url);
S_met_renamed = rename_underway_variables_to_SWIFT_compliance(S_met);

% Navigation variables: Here COG, SOG, heading, water depth
nav_URL = 'https://coriolix.sikuliaq.alaska.edu/erddap/tabledap/binned_default_nav_rolling.csv?time,latitude,longitude,parameter_25,parameter_08,parameter_07,parameter_06&time>=2025-10-17T00:00&time<=2025-11-12T23:59&orderBy("time")';
nav_metadata_URL = 'https://coriolix.sikuliaq.alaska.edu/erddap/info/binned_default_nav_rolling/index.csv';
[S_nav, S_nav_info] = get_SKQ_underway_data(nav_URL, nav_metadata_URL);
S_nav_renamed = rename_underway_variables_to_SWIFT_compliance(S_nav);

% 'Other' variables: Here cholorphyll and nitrate
other_URL = 'https://coriolix.sikuliaq.alaska.edu/erddap/tabledap/binned_default_other_rolling.csv?time,latitude,longitude,parameter_30,parameter_25&time>=2025-10-17T00:00&time<=2025-11-12T23:59&orderBy("time")';
other_metadata_URL = 'https://coriolix.sikuliaq.alaska.edu/erddap/info/binned_default_other_rolling/index.csv';
[S_other, S_other_info] = get_SKQ_underway_data(other_URL, other_metadata_URL);
S_other_renamed = rename_underway_variables_to_SWIFT_compliance(S_other);

% combine renamed structs
SWIFT = combineStructs(S_flow_renamed, S_met_renamed);
SWIFT = combineStructs(SWIFT, S_nav_renamed);
SWIFT = combineStructs(SWIFT, S_other_renamed);

%%
function [S, S_info] = get_SKQ_underway_data(data_url, metadata_url)
    % Read metadata from url
    metadata = readtable(metadata_url, 'TextType','string');
    
    % Make a map of info for each parameter - note that metadata file
    % includes info for all parameters, including those we don't want
    metadata_map = build_metadata_map(metadata);  
    
    % Read the data as text, split apart into lines for each timestamp 
    rawText = webread(data_url, weboptions('ContentType','text'));
    lines = splitlines(string(rawText));
    lines(lines=="") = [];
    
    % Header and data lines
    data_types = lines(1); %Contains what parameters we've requested and downloaded
    
    % Split apart what data types are present
    data_types = split(data_types, ',');
    
    % First 3 are known: time, latitude, longitude; extract remaining 'parameter_[number]' columns
    param_nums = cellstr(data_types(4:end));

    % Make metadata info output struct for metadata for each parameter
    % we've actually requested
    S_info = make_metadata_out(param_nums, metadata_map);  

    %Extract units
    data_units = lines(2);
    data_units = cellstr(split(data_units, ',')); %Units for all data_types (including UTC, lon, lat)
    
    % Iterate through lines of data for each timestep
    % Extract time bin mean and std for each parameter
    data_lines = lines(3:end);
        
    %Make output struct for data
    S = setup_data_out(S_info, data_lines);
    
    %Extract data at each timestep
    for j = 1:numel(data_lines)
        cur_line = data_lines(j); %Text that corresponds to the current timestep
    
        % Matches: time,latitude,longitude,"{...}",...
        pat = '^([^,]+),([^,]+),([^,]+),(.*)$';
        toks = regexp(cur_line, pat, 'tokens', 'once');
    
        time_str  = (toks{1}); % i.e., 2025-10-19T00:00:00Z
        time = datetime(time_str, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z','TimeZone','UTC');
        
        %Write lat/lon/time for this line of data
        S(j).time = datenum(time);
        S(j).latitude = str2double(toks{2});
        S(j).longitude = str2double(toks{3});
    
        %Need to parse this JSON to get data for timestep
        data_this_line = string(toks{4});
    
        %From https://coriolix.sikuliaq.alaska.edu/data/format/, there are a
        %variety of average types avavilable. We'll keep b: best
    
        % Extract the entire JSON block
        patJSON = '"(\{""a"".*?\})"';                         
        json_blocks = regexp(data_this_line, patJSON, 'tokens');
    
        assert(numel(json_blocks) == numel(param_nums), 'Number of data columns does not match number of selected instruments. Did you check that all instruments selected on Coriolix are active?') %There should be one block per parameter
       
        %Iterate through each parameter (i.e., parameter_01, parameter_02)
        for k = 1:numel(param_nums)
            json_str = json_blocks{k}{1}; %JSON string for this parameter
            json_str = replace(json_str, '""','"'); % get rid of doubled quotes
        
            %Uses matlab function to unpack JSON into a struct
            J = jsondecode(json_str); % has fields a, b, c, fa, fb, fc, sa, sb, sc for different data qualities
        
            % Figure out what entries in J to keep
            param_num  = param_nums{k};
            map_to_means = metadata_map.(param_num).array_contents;
            idx_mean = find(strcmp(map_to_means, 'mean')); %Can optionally also retain the other available statistics: std, min, max, num, spotval, median
            idx_std = find(strcmp(map_to_means, 'std'));
    
            % Use best bin averaged value ("b") and keep mean & std only
            mean_val = J.b(idx_mean); %If a different data quality is preferred, update here (i.e., to J.a for all bins)
            std_val = J.b(idx_std);
    
            % Save data for this parameter
            S(j).(S_info(k).field_mean) = mean_val;
            S(j).(S_info(k).field_std) = std_val;
        end
    
    end
end


% From metadata table, make a map for each 'parameter_##' -> parse long_name, units, array_contents, sensor_description
function metadata_map = build_metadata_map(metadata)

    % All parameter variables in the metadata - may not have requested all
    % of these in the data url
    isParamVar = metadata.RowType=="variable" & startsWith(metadata.VariableName,"parameter_");
    param_nums = unique(metadata.VariableName(isParamVar)); %i.e. parameter_01, parameter_02

    metadata_map = struct(); %initialize
    for j = 1:numel(param_nums) %iterate through each parameter 
        
        % All rows of info for this parameter
        p = param_nums(j);
        rows = metadata.RowType=="attribute" & metadata.VariableName==p;

        % Get info for this parameter number
        long_name = metadata.Value(rows & metadata.AttributeName=="long_name");
        long_name = strrep(long_name, ' ', '_');
        long_name = strrep(long_name, '-', '');

        units    = metadata.Value(rows & metadata.AttributeName=="units");
        ac_raw   = metadata.Value(rows & metadata.AttributeName=="array_contents"); %mean, std, median, etc
        ac_list = strtrim(split(ac_raw, ','));
        sensor    = metadata.Value(rows & metadata.AttributeName=="sensor_description");
  
        % Add to map for this parameter p
        metadata_map.(p).long_name = long_name;
        metadata_map.(p).units = units;
        metadata_map.(p).sensor_description = sensor;        
        metadata_map.(p).array_contents = cellstr(ac_list).'; 
    end
end

% Make metadata info struct using the metadata_map (all parameters) and the
% parameters we've actually requested (param_nums)
function S_info = make_metadata_out(param_nums, metadata_map);
    
    struct('param_num',[],'long_name',[],'units',[],'sensor_description',[], ...
                       'field_mean',[],'field_std',[]);
    
    for k = 1:numel(param_nums)
        pnum = param_nums{k};
        meta = metadata_map.(pnum);
        base = string(meta.long_name);
    
        S_info(k).param_num = pnum; %Specific to met/flowthrough etc. I.e. there can be met param_01 and also flowthrough param_01
        S_info(k).long_name = base;
        S_info(k).units = string(meta.units);
        S_info(k).sensor_description = string(meta.sensor_description);
        S_info(k).field_mean = base + "_b_mean"; %Keep "b" best data quality mean 
        S_info(k).field_std = base + "_b_std"; %Keed std from time bin averaging
    end
end

%Initialize data struct
function S = setup_data_out(S_info, data_lines)

    base_template = struct('time',NaN,'latitude',NaN,'longitude',NaN);
    for k = 1:numel(S_info)
        base_template.(S_info(k).field_mean) = NaN;
        base_template.(S_info(k).field_std)  = NaN;
    end
    S = repmat(base_template, numel(data_lines), 1);
end

% Map the Coriolix long_names to SWIFT defaults to make structs that are compatible with
% existing SWIFT routines
function S = rename_underway_variables_to_SWIFT_compliance(S)
    flds = fields(S);
    for f = 1:length(flds)
        
        %Start making new name - retain _b_mean etc. 
        tok = regexp(flds{f}, '(_.*)$', 'tokens', 'once');
        if isempty(tok) %i.e. for lat, lon, time
            fld_ending = '';
        else
            fld_ending = tok{1}; %e.g., '_b_mean' or  '_b_std'
        end

        if strcmp(fld_ending, '_b_std');
            S = rmfield(S, flds{f});
            continue
        end

        fld_name = lower(flds{f});
        fld_name = strrep(fld_name, '_', ' ');  

        %Note that renameStructField does not work with S because it's a
        %struct array
        if contains(fld_name, 'longitude')
            new_var_name = 'lon';
                                
        elseif contains(fld_name, 'latitude')
            new_var_name = 'lat';
            
        elseif contains(fld_name, 'time')
            new_var_name = 'time';

        elseif contains(fld_name, 'surface temperature')
            new_var_name = 'watertemp';

        elseif contains(fld_name, 'water skin')
            new_var_name = 'waterskintemp';

        elseif contains(fld_name, 'salinity') %Actual long name is Thermosalinograph_Salinity so still need to change           
            new_var_name = 'salinity';

        elseif contains(fld_name, 'barometric pressure')            
            new_var_name = 'airpres';

        elseif contains(fld_name, 'air temperature')            
            new_var_name = 'airtemp';

        elseif contains(fld_name, 'relative humidity')            
            new_var_name = 'relhumidity';

        elseif contains(fld_name, 'precipitation rate')            
            new_var_name = 'rainrate';

        elseif contains(fld_name, ['radiation', strrep(fld_ending, '_', '')])  %Have to distinguish shortwave (called "radiation") from longwave radiation       
            new_var_name = 'solarrad';

        elseif contains(fld_name, 'longwave')            
            new_var_name = 'LWdown';

        elseif contains(fld_name, 'true wind direction') %Check if we actually want the relative +  COG for SWIFT compatibility            
            new_var_name = 'winddirT';

        elseif contains(fld_name, 'true wind speed')            
            new_var_name = 'windspd';

        elseif contains(fld_name, 'PAR')            
            new_var_name = 'PAR';

        elseif contains(fld_name, 'speed over ground')            
            new_var_name = 'driftspd';

        elseif contains(fld_name, 'course over ground')            
            new_var_name = 'driftdirT';
            
        else
            new_var_name = ''; %Don't rename;
        end

        if ~isempty(new_var_name)
            %Add the new field name
            [S.([new_var_name])] = S.(flds{f});
        end

        %Remove the old field name
        if ~strcmp(new_var_name, flds{f}) %i.e., time doesn;t change its name but we don't want to remove it
            S = rmfield(S, flds{f});
        end
    end
end

% Combine structs for flowthrough/met/etc. horizontally after  matching lat/lon/time
function S_comb = combineStructs(S1, S2)


    %Find common timesteps
    t1 = [S1.time];
    t2 = [S2.time];
    [common_times, ia, ib] = intersect(t1, t2, 'stable');

    % Keep only the overlapping timesteps
    S1 = S1(ia);
    S2 = S2(ib);
    
    % check same length
    if numel(S1) ~= numel(S2)
        error('Structs must have same length after aligning timesteps');
    end

    % check lat/lon/time actually match at common timesteps
    if any([S1.lat] ~= [S2.lat]) || any([S1.lon] ~= [S2.lon]) || any([S1.time] ~= [S2.time])
        error('lat/lon/time fields do not match exactly - check time ranges in API urls');
    end
    
    % combine fields (skip duplicates)
    S_comb = S1;
    flds = setdiff(fieldnames(S2), {'lat','lon','time'});
    for f = 1:numel(flds)
        [S_comb.(flds{f})] = deal(S2.(flds{f}));
    end
end
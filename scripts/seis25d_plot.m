%% seis25d_plot.m
% Consolidated SEIS25D_GQG_v2 plotting script.
%
% Edit the USER SETTINGS section, then press Run in MATLAB.


clear; close all; clc;

%% ========================================================================
% USER SETTINGS
% ========================================================================

% Main SEIS25D_GQG_v2 run folder. Change this path before plotting another
% result.
input_folder = ['C:\Users\sedar\Desktop\SEIS25D_GQG_github\' ...
    'SEIS25D_GQG_v2\runs\Small_TTI\output_20260804_145749'];

% Figures are saved here, inside input_folder when this is a relative name.
output_folder_name = 'Matlab_plots';

% Case-specific controls.
well_x = 2205;                 % x-location for vertical profiles [m]
depth_limits = [-2510, 100];   % depth range for plots [m]
x_limits = [];                 % example: [0, 4500]; use [] for automatic
top_boundary_zone_thickness = 0;

% Figure export quality.
dpi = 150;

% Plot switches.
make_model_maps        = false;
make_vertical_profiles = false;
make_gradients         = false;
make_diagnostics       = true;
make_spectral_data     = true;
make_data_heatmaps     = true;

% Spectral data plot controls. These read files such as GT0_1_0.txt,
% GTO_1_0.txt, and GO_1_1.txt, then plot real/imag versus receiver ID.
spectral_source_id = 1;
spectral_frequency = 1;

% Limit to selected parameters if wanted. Use {} for all discovered.
parameters_to_plot = {};       % example: {'C33', 'Q33'}

% Optional measured/true well profiles. Each file should have columns:
% depth, value. Leave empty if not needed.
true_profile_files = struct();
% true_profile_files.C33 = 'C:\path\to\S13_C33_QC_well.csv';
% true_profile_files.Q33 = 'C:\path\to\Qp_QC.txt';

%% ========================================================================
% PARAMETER NAMES AND UNIT LABELS
% ========================================================================

param_labels = parameter_labels();

%% ========================================================================
% MAIN
% ========================================================================

set_default_plot_style();

if ~isfolder(input_folder)
    error('input_folder does not exist:\n%s', input_folder);
end

if isabsolute_path(output_folder_name)
    output_folder = output_folder_name;
else
    output_folder = fullfile(input_folder, output_folder_name);
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

fprintf('Input folder: %s\n', input_folder);

pcbb_catalog = find_final_models(input_folder);
grad_catalog = find_gradient_files(input_folder);

fprintf('Found %d PCBB files.\n', numel(pcbb_catalog));
fprintf('Found %d gradient files.\n', numel(grad_catalog));

if make_model_maps && ~isempty(pcbb_catalog)
    plot_model_maps(input_folder, output_folder, pcbb_catalog, ...
        parameters_to_plot, param_labels, x_limits, depth_limits, ...
        top_boundary_zone_thickness, dpi);
end

if make_vertical_profiles && ~isempty(pcbb_catalog)
    plot_vertical_profiles(input_folder, output_folder, pcbb_catalog, ...
        parameters_to_plot, param_labels, true_profile_files, well_x, ...
        depth_limits, dpi);
end

if make_gradients && ~isempty(grad_catalog)
    plot_gradient_series(input_folder, output_folder, grad_catalog, ...
        parameters_to_plot, param_labels, x_limits, depth_limits, ...
        top_boundary_zone_thickness, dpi);
end

if make_diagnostics
    plot_diagnostics(input_folder, output_folder, dpi);
end

if make_spectral_data
    plot_spectral_data(input_folder, output_folder, spectral_source_id, spectral_frequency, dpi);
end

if make_data_heatmaps
    plot_data_heatmaps(input_folder, output_folder, dpi);
end

fprintf('Done.\n');

%% ========================================================================
% FILE DISCOVERY FUNCTIONS
% ========================================================================

function catalog = find_final_models(folder_name)
    files = dir(fullfile(folder_name, 'PCBB_*.dat'));
    catalog = empty_catalog();

    for k = 1:numel(files)
        item = parse_iteration_file(files(k), 'PCBB');
        if ~isempty(item)
            catalog(end+1) = item; %#ok<AGROW>
        end
    end

    catalog = sort_catalog(catalog);
end

function catalog = find_gradient_files(folder_name)
    files = dir(fullfile(folder_name, 'GRAD*.dat'));
    catalog = empty_catalog();

    for k = 1:numel(files)
        item = parse_iteration_file(files(k), 'GRAD');
        if ~isempty(item)
            catalog(end+1) = item; %#ok<AGROW>
        end
    end

    catalog = sort_catalog(catalog);
end

function item = parse_iteration_file(file_info, file_type)
    item = [];

    [~, name_only, ext] = fileparts(file_info.name);
    if ~strcmpi(ext, '.dat')
        return;
    end

    if strcmp(file_type, 'PCBB')
        family = 'PCBB';
        raw = erase(name_only, 'PCBB_');
    else
        parts0 = split(name_only, '_');
        if numel(parts0) < 4
            return;
        end
        family = char(parts0{1});
        raw = extractAfter(name_only, [family, '_']);
    end

    % Example raw strings:
    %   {{C33}_03.05_IT01}
    %   {C33}_03.05_IT01
    raw = erase(raw, {'{', '}'});
    parts = split(raw, '_');
    if numel(parts) < 3
        return;
    end

    item = struct();
    item.path = fullfile(file_info.folder, file_info.name);
    item.name = file_info.name;
    item.family = family;
    item.param = char(parts{1});
    item.freq = char(parts{2});
    item.itLabel = char(parts{3});
    item.itNum = parse_iteration_number(item.itLabel);
end

function catalog = empty_catalog()
    catalog = struct( ...
        'path', {}, ...
        'name', {}, ...
        'family', {}, ...
        'param', {}, ...
        'freq', {}, ...
        'itLabel', {}, ...
        'itNum', {});
end

function catalog = sort_catalog(catalog)
    if isempty(catalog)
        return;
    end

    param = {catalog.param}.';
    family = {catalog.family}.';
    freq_num = zeros(numel(catalog), 1);
    it_num = zeros(numel(catalog), 1);
    for k = 1:numel(catalog)
        freq_num(k) = str2double(catalog(k).freq);
        if isnan(freq_num(k))
            freq_num(k) = inf;
        end
        it_num(k) = catalog(k).itNum;
    end

    T = table(family, param, freq_num, it_num, (1:numel(catalog)).', ...
        'VariableNames', {'family', 'param', 'freq', 'iter', 'idx'});
    T = sortrows(T, {'family', 'param', 'freq', 'iter'});
    catalog = catalog(T.idx);
end

function n = parse_iteration_number(label)
    token = regexp(label, '\d+$', 'match', 'once');
    if isempty(token)
        n = 0;
    else
        n = str2double(token);
    end
end

function files = find_parameter_files(folder_name, prefix)
    files = struct();
    d = dir(fullfile(folder_name, [prefix, '_*.dat']));

    for k = 1:numel(d)
        token = regexp(d(k).name, '\{([^}]+)\}', 'tokens', 'once');
        if isempty(token)
            continue;
        end
        param = erase(token{1}, {'{', '}'});
        files.(param) = fullfile(d(k).folder, d(k).name);
    end
end

function data_files = find_data_files(folder_name)
    data_files = struct();
    data_files.go = dir(fullfile(folder_name, 'GO_*.txt'));
    data_files.gt0 = dir(fullfile(folder_name, 'GT0_*.txt'));
    data_files.gto = dir(fullfile(folder_name, 'GTO_*.txt'));
    data_files.residual = dir(fullfile(folder_name, 'out_resid*.txt'));
end

function files = find_spectral_files(folder_name, frequency)
    data_files = find_data_files(folder_name);
    all_files = [data_files.gt0; data_files.gto; data_files.go];
    files = all_files([]);

    for k = 1:numel(all_files)
        if data_file_matches_frequency(all_files(k).name, frequency)
            files(end+1) = all_files(k); %#ok<AGROW>
        end
    end
end

function tf = data_file_matches_frequency(file_name, frequency)
    token = regexp(file_name, '^(GO|GT0|GTO)_([-+]?\d*\.?\d+)_', 'tokens', 'once');
    if isempty(token)
        tf = false;
        return;
    end

    file_frequency = str2double(token{2});
    wanted_frequency = str2double(string(frequency));
    tf = abs(file_frequency - wanted_frequency) <= max(1e-8, abs(wanted_frequency) * 1e-8);
end

%% ========================================================================
% FILE READING FUNCTIONS
% ========================================================================

function grid = read_model_grid(file_name)
    A = readmatrix(file_name, 'FileType', 'text');
    grid.x = A(2:end, 1);
    grid.z = A(1, 2:end);
    grid.values = A(2:end, 2:end);
    grid.values(grid.values >= 1.0e20) = NaN;
end

function [z_values, profile] = read_profile_at_x(file_name, well_x)
    grid = read_model_grid(file_name);
    [~, ix] = min(abs(grid.x - well_x));
    z_values = grid.z;
    profile = grid.values(ix, :);
end

function [z_values, values] = read_true_profile(file_name, depth_limits)
    A = readmatrix(file_name);
    z_values = A(:, 1);
    values = A(:, 2);

    if ~isempty(depth_limits)
        keep = z_values >= depth_limits(1) & z_values <= depth_limits(2);
        z_values = z_values(keep);
        values = values(keep);
    end
end

function topo = read_topography(folder_name)
    topo = [];
    file_name = fullfile(folder_name, 'm_Top.dat');
    if ~isfile(file_name)
        return;
    end
    A = readmatrix(file_name, 'FileType', 'text');
    topo.x = A(:, end-1);
    topo.z = A(:, end);
end

function sr = read_sr_geometry(folder_name)
    sr = [];
    file_name = fullfile(folder_name, 'm_SR.dat');
    if ~isfile(file_name)
        file_name = fullfile(folder_name, 'm_S-R.dat');
    end
    if ~isfile(file_name)
        return;
    end

    A = readmatrix(file_name, 'FileType', 'text');
    is_source = A(:, 1) == 1;
    x = A(:, end-1);
    z = A(:, end);

    sr.x_src = x(is_source);
    sr.z_src = z(is_source);
    sr.x_rec = x(~is_source);
    sr.z_rec = z(~is_source);
end

function data = read_go_or_residual(file_name)
    A = readmatrix(file_name, 'FileType', 'text');

    if size(A, 2) >= 14
        data.sid = A(:, 4);
        data.scomp = A(:, 5);
        data.rid = A(:, 6);
        data.rcomp = A(:, 7);
        data.real = A(:, 12);
        data.imag = A(:, 13);
        data.amp = A(:, 14);
    else
        data.sid = A(:, 3);
        data.scomp = A(:, 4);
        data.rid = A(:, 5);
        data.rcomp = A(:, 6);
        data.real = A(:, 11);
        data.imag = A(:, 12);
        data.amp = hypot(data.real, data.imag);
    end
end

%% ========================================================================
% PLOTTING FUNCTIONS
% ========================================================================

function plot_model_maps(folder_name, output_folder, catalog, selected_params, labels, x_limits, z_limits, top_zone, dpi)
    true_models = find_parameter_files(folder_name, 'mT');
    start_models = find_parameter_files(folder_name, 'mI');
    groups = unique_groups(catalog);

    for g = 1:size(groups, 1)
        family = groups{g, 1}; %#ok<NASGU>
        param = groups{g, 2};
        freq = groups{g, 3};

        if ~should_plot_param(param, selected_params)
            continue;
        end

        items = select_group(catalog, groups(g, :));
        final_item = items(end);

        ref_file = '';
        ref_name = '';
        if isfield(true_models, param)
            ref_file = true_models.(param);
            ref_name = 'True';
        elseif isfield(start_models, param)
            ref_file = start_models.(param);
            ref_name = 'Start';
        end

        if isempty(ref_file)
            warning('Skipping %s because no mT/mI reference file was found.', param);
            continue;
        end

        ref_grid = read_model_grid(ref_file);
        final_grid = read_model_grid(final_item.path);

        diff_abs = final_grid.values - ref_grid.values;
        diff_pct = 100 * diff_abs ./ ref_grid.values;

        fig = figure('Color', 'w', 'Name', ['Model ', param, ' ', freq], ...
            'Position', [100, 80, 1100, 780]);
        tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        clim_ref = finite_minmax(ref_grid.values);

        ax = nexttile(tl, 1);
        plot_grid_panel(ax, ref_grid, label_for_param(labels, param), clim_ref, 'jet');
        title(ax, ['(a) ', ref_name], 'FontWeight', 'normal');
        overlay_topography_and_acquisition(ax, folder_name, top_zone);
        apply_map_limits(ax, x_limits, z_limits);

        ax = nexttile(tl, 2);
        plot_grid_panel(ax, final_grid, label_for_param(labels, param), clim_ref, 'jet');
        title(ax, ['(b) Final ', freq, ' Hz ', final_item.itLabel], 'FontWeight', 'normal');
        overlay_topography_and_acquisition(ax, folder_name, top_zone);
        apply_map_limits(ax, x_limits, z_limits);

        ax = nexttile(tl, 3);
        plot_grid_panel(ax, make_grid(final_grid.x, final_grid.z, diff_pct), '%', [-10, 10], 'redblue');
        title(ax, '(c) Percent difference', 'FontWeight', 'normal');
        overlay_topography_and_acquisition(ax, folder_name, top_zone);
        apply_map_limits(ax, x_limits, z_limits);

        ax = nexttile(tl, 4);
        vm = robust_abs_limit(diff_abs, 95, 1);
        plot_grid_panel(ax, make_grid(final_grid.x, final_grid.z, diff_abs), label_for_param(labels, param), [-vm, vm], 'redblue');
        title(ax, '(d) Absolute difference', 'FontWeight', 'normal');
        overlay_topography_and_acquisition(ax, folder_name, top_zone);
        apply_map_limits(ax, x_limits, z_limits);

        out_png = fullfile(output_folder, ['model_compare_', safe_name(param), '_', safe_name(freq), '.png']);
        save_png(fig, out_png, dpi);
        close(fig);
    end
end

function plot_vertical_profiles(folder_name, output_folder, catalog, selected_params, labels, true_profiles, well_x, z_limits, dpi)
    true_models = find_parameter_files(folder_name, 'mT');
    start_models = find_parameter_files(folder_name, 'mI');
    groups = unique_groups(catalog);

    for g = 1:size(groups, 1)
        param = groups{g, 2};
        freq = groups{g, 3};

        if ~should_plot_param(param, selected_params)
            continue;
        end

        items = select_group(catalog, groups(g, :));
        profiles_for_xlim = {};

        fig = figure('Color', 'w', 'Name', ['Profile ', param, ' ', freq], ...
            'Position', [150, 80, 680, 820]);
        ax = axes(fig); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

        if isfield(true_models, param)
            [z, p] = read_profile_at_x(true_models.(param), well_x);
            plot(ax, p, z, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True');
            profiles_for_xlim{end+1} = p; %#ok<AGROW>
        end

        if isfield(true_profiles, param) && isfile(true_profiles.(param))
            [z, p] = read_true_profile(true_profiles.(param), z_limits);
            plot(ax, p, z, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Well/QC');
            profiles_for_xlim{end+1} = p; %#ok<AGROW>
        end

        if isfield(start_models, param)
            [z, p] = read_profile_at_x(start_models.(param), well_x);
            plot(ax, p, z, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Start');
            profiles_for_xlim{end+1} = p; %#ok<AGROW>
        end

        cmap = hsv(numel(items));
        for k = 1:numel(items)
            [z, p] = read_profile_at_x(items(k).path, well_x);
            plot(ax, p, z, '-', 'LineWidth', 1.3, 'Color', cmap(k, :), ...
                'DisplayName', ['Iter ', num2str(items(k).itNum)]);
            profiles_for_xlim{end+1} = p; %#ok<AGROW>
        end

        xlabel(ax, label_for_param(labels, param));
        ylabel(ax, 'Depth [m]');
        if ~isempty(z_limits)
            ylim(ax, z_limits);
        end
        set_profile_xlim(ax, profiles_for_xlim);

        title(ax, sprintf('%s at x = %g m, %s Hz', param, well_x, freq), ...
            'Interpreter', 'none', 'FontWeight', 'normal');
        legend(ax, 'show', 'Location', 'eastoutside');

        out_png = fullfile(output_folder, ...
            sprintf('profile_%s_%s_x%g.png', safe_name(param), safe_name(freq), well_x));
        save_png(fig, out_png, dpi);
        close(fig);
    end
end

function plot_gradient_series(folder_name, output_folder, catalog, selected_params, labels, x_limits, z_limits, top_zone, dpi)
    groups = unique_groups(catalog);
    grad_folder = fullfile(output_folder, 'Grad');
    if ~isfolder(grad_folder)
        mkdir(grad_folder);
    end

    for g = 1:size(groups, 1)
        family = groups{g, 1};
        param = groups{g, 2};
        freq = groups{g, 3};

        if ~should_plot_param(param, selected_params)
            continue;
        end

        items = select_group(catalog, groups(g, :));
        n = numel(items);
        ncols = ceil(sqrt(n));
        nrows = ceil(n / ncols);

        all_values = [];
        for k = 1:n
            grid = read_model_grid(items(k).path);
            all_values = [all_values; grid.values(:)]; %#ok<AGROW>
        end
        vm = robust_abs_limit(all_values, 99, 1);

        fig = figure('Color', 'w', 'Name', [family, ' ', param, ' ', freq], ...
            'Position', [80, 60, 1200, 820]);
        tl = tiledlayout(fig, nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

        for k = 1:n
            ax = nexttile(tl, k);
            grid = read_model_grid(items(k).path);
            plot_grid_panel(ax, grid, label_for_param(labels, param), [-vm, vm], 'redblue');
            overlay_topography_and_acquisition(ax, folder_name, top_zone);
            apply_map_limits(ax, x_limits, z_limits);
            title(ax, sprintf('(%c) %s', 'a' + k - 1, items(k).itLabel), ...
                'FontWeight', 'normal');
        end

        title(tl, sprintf('%s Gradient at %s Hz for %s', family, freq, param), ...
            'Interpreter', 'none', 'FontWeight', 'normal');

        out_png = fullfile(grad_folder, [family, '_', safe_name(param), '_', safe_name(freq), '.png']);
        save_png(fig, out_png, dpi);
        close(fig);
    end
end

function plot_diagnostics(folder_name, output_folder, dpi)
    diag_files = [dir(fullfile(folder_name, 'out_diag*.txt')); dir(fullfile(folder_name, 'out_diagt*.txt'))];
    if isempty(diag_files)
        fprintf('No out_diag*.txt files found.\n');
        return;
    end

    for f = 1:numel(diag_files)
        file_name = fullfile(diag_files(f).folder, diag_files(f).name);
        A = readmatrix(file_name, 'FileType', 'text');
        A = A(all(isfinite(A), 2), :);
        if size(A, 2) < 3
            continue;
        end

        freq = A(:, 1);
        iter = A(:, 2);
        fcost = A(:, 3);
        has_grad = size(A, 2) >= 4;
        has_resid = size(A, 2) >= 5;
        has_time = size(A, 2) >= 6;

        fig = figure('Color', 'w', 'Name', ['Diagnostics ', diag_files(f).name], ...
            'Position', [100, 80, 1050, 780]);
        tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

        ax1 = nexttile(tl, 1); hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
        ax2 = nexttile(tl, 2); hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
        ax3 = nexttile(tl, 3); hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on');
        ax4 = nexttile(tl, 4); hold(ax4, 'on'); box(ax4, 'on'); grid(ax4, 'on');
        ax5 = nexttile(tl, 5); hold(ax5, 'on'); box(ax5, 'on'); grid(ax5, 'on');

        freqs = unique(freq, 'stable');
        cmap = frequency_colormap(max(numel(freqs), 2));

        for k = 1:numel(freqs)
            mask = freq == freqs(k);
            [itk, order] = sort(iter(mask));
            ck = fcost(mask); ck = ck(order);
            base = ck(1);
            if base == 0
                base = eps;
            end

            semilogy(ax1, itk, 100 * ck / base, '-o', ...
                'LineWidth', 1.3, 'MarkerSize', 5, 'Color', cmap(k, :), ...
                'DisplayName', sprintf('f = %.2f Hz', freqs(k)));
            plot(ax2, itk, 100 * (base - ck) / base, '-d', ...
                'LineWidth', 1.3, 'MarkerSize', 4, 'Color', cmap(k, :));

            if has_grad
                gk = A(mask, 4); gk = gk(order);
                semilogy(ax3, itk, gk, '-s', 'LineWidth', 1.3, ...
                    'MarkerSize', 4, 'Color', cmap(k, :));
            end
            if has_resid
                rk = A(mask, 5); rk = rk(order);
                plot(ax4, itk, rk, '-^', 'LineWidth', 1.3, ...
                    'MarkerSize', 4, 'Color', cmap(k, :));
            end
            if has_time
                tk = A(mask, 6); tk = tk(order);
                plot(ax5, itk, tk, '-o', 'LineWidth', 1.3, ...
                    'MarkerSize', 4, 'Color', cmap(k, :));
            end
        end

        xlabel(ax1, 'Iterations'); ylabel(ax1, 'Normalized \phi [%]');
        xlabel(ax2, 'Iterations'); ylabel(ax2, 'Cumulative \Delta\phi [%]');
        xlabel(ax3, 'Iterations'); ylabel(ax3, '||gradient||');
        xlabel(ax4, 'Iterations'); ylabel(ax4, 'Residual RMS');
        xlabel(ax5, 'Iterations'); ylabel(ax5, 'Time per iteration [min]');
        title(ax1, '(a)', 'FontWeight', 'normal');
        title(ax2, '(b)', 'FontWeight', 'normal');
        title(ax3, '(c)', 'FontWeight', 'normal');
        title(ax4, '(d)', 'FontWeight', 'normal');
        title(ax5, '(e)', 'FontWeight', 'normal');
        legend(ax1, 'show', 'Location', 'best');

        out_png = fullfile(output_folder, ['diagnostic_', erase(diag_files(f).name, '.txt'), '.png']);
        save_png(fig, out_png, dpi);
        close(fig);
    end
end

function plot_data_heatmaps(folder_name, output_folder, dpi)
    files = find_data_files(folder_name);
    chosen = [];
    if ~isempty(files.residual)
        chosen = files.residual(1);
    elseif ~isempty(files.gt0)
        chosen = files.gt0(1);
    elseif ~isempty(files.gto)
        chosen = files.gto(1);
    elseif ~isempty(files.go)
        chosen = files.go(1);
    end

    if isempty(chosen)
        fprintf('No residual, GT0, GTO, or GO files found for heatmaps.\n');
        return;
    end

    file_name = fullfile(chosen.folder, chosen.name);
    data = read_go_or_residual(file_name);

    sids = unique(data.sid, 'stable');
    rids = unique(data.rid, 'stable');
    scomps = unique(data.scomp, 'stable');
    rcomps = unique(data.rcomp, 'stable');

    ncols = numel(scomps) * numel(rcomps);
    fig = figure('Color', 'w', 'Name', ['Heatmap ', chosen.name], ...
        'Position', [80, 80, max(900, 320*ncols), 780]);
    tl = tiledlayout(fig, 3, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

    col = 0;
    for is = 1:numel(scomps)
        for ir = 1:numel(rcomps)
            col = col + 1;
            mask = data.scomp == scomps(is) & data.rcomp == rcomps(ir);

            Mreal = sid_rid_matrix(data.sid(mask), data.rid(mask), data.real(mask), sids, rids);
            Mimag = sid_rid_matrix(data.sid(mask), data.rid(mask), data.imag(mask), sids, rids);
            Mamp = sid_rid_matrix(data.sid(mask), data.rid(mask), data.amp(mask), sids, rids);

            ax = nexttile(tl, col);
            imagesc(ax, sids, rids, Mreal.'); axis(ax, 'xy'); set(ax, 'YDir', 'reverse');
            colormap(ax, blue_white_red(256)); colorbar(ax); title(ax, sprintf('s%d r%d real', scomps(is), rcomps(ir)));
            xlabel(ax, 'Source ID'); ylabel(ax, 'Receiver ID');

            ax = nexttile(tl, ncols + col);
            imagesc(ax, sids, rids, Mimag.'); axis(ax, 'xy'); set(ax, 'YDir', 'reverse');
            colormap(ax, blue_white_red(256)); colorbar(ax); title(ax, 'imag');
            xlabel(ax, 'Source ID'); ylabel(ax, 'Receiver ID');

            ax = nexttile(tl, 2*ncols + col);
            imagesc(ax, sids, rids, Mamp.'); axis(ax, 'xy'); set(ax, 'YDir', 'reverse');
            colormap(ax, parula); colorbar(ax); title(ax, 'amplitude');
            xlabel(ax, 'Source ID'); ylabel(ax, 'Receiver ID');
        end
    end

    out_png = fullfile(output_folder, ['heatmap_', erase(chosen.name, '.txt'), '.png']);
    save_png(fig, out_png, dpi);
    close(fig);
end

function plot_spectral_data(folder_name, output_folder, source_id, frequency, dpi)
    spectral_files = find_spectral_files(folder_name, frequency);
    if isempty(spectral_files)
        fprintf('No GT0/GTO/GO files found for frequency %g.\n', frequency);
        return;
    end

    loaded = struct('name', {}, 'data', {});
    receiver_components = [];

    for k = 1:numel(spectral_files)
        file_name = fullfile(spectral_files(k).folder, spectral_files(k).name);
        data = read_go_or_residual(file_name);
        source_mask = data.sid == source_id;
        if ~any(source_mask)
            continue;
        end
        loaded(end+1).name = erase(spectral_files(k).name, '.txt'); %#ok<AGROW>
        loaded(end).data = data;
        receiver_components = [receiver_components; unique(data.rcomp(source_mask), 'stable')]; %#ok<AGROW>
    end

    if isempty(loaded)
        fprintf('No spectral rows found for source ID %g at frequency %g.\n', source_id, frequency);
        return;
    end

    receiver_components = unique(receiver_components, 'stable');
    ncols = numel(receiver_components);

    fig = figure('Color', 'w', 'Name', sprintf('Spectral source %g freq %g', source_id, frequency), ...
        'Position', [80, 80, max(900, 380*ncols), 720]);
    tl = tiledlayout(fig, 2, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for c = 1:ncols
        rcomp = receiver_components(c);
        ax_real = nexttile(tl, c); hold(ax_real, 'on'); box(ax_real, 'on'); grid(ax_real, 'on');
        ax_imag = nexttile(tl, ncols + c); hold(ax_imag, 'on'); box(ax_imag, 'on'); grid(ax_imag, 'on');

        for k = 1:numel(loaded)
            data = loaded(k).data;
            mask = data.sid == source_id & data.rcomp == rcomp;
            if ~any(mask)
                continue;
            end

            receiver_id = data.rid(mask);
            real_part = data.real(mask);
            imag_part = data.imag(mask);
            [receiver_id, order] = sort(receiver_id);
            real_part = real_part(order);
            imag_part = imag_part(order);

            plot(ax_real, receiver_id, real_part, '-', 'LineWidth', 1.2, ...
                'DisplayName', loaded(k).name);
            plot(ax_imag, receiver_id, imag_part, '-', 'LineWidth', 1.2, ...
                'DisplayName', loaded(k).name);
        end

        title(ax_real, sprintf('Receiver component %g', rcomp), 'FontWeight', 'normal');
        ylabel(ax_real, 'Real');
        ylabel(ax_imag, 'Imaginary');
        xlabel(ax_imag, 'Receiver ID');
        legend(ax_real, 'show', 'Location', 'best');
    end

    title(tl, sprintf('Spectral data, source %g, frequency %g', source_id, frequency), ...
        'FontWeight', 'normal');
    out_png = fullfile(output_folder, sprintf('spectral_source%g_freq%g.png', source_id, frequency));
    save_png(fig, out_png, dpi);
    close(fig);
end

%% ========================================================================
% SMALL PLOTTING HELPERS
% ========================================================================

function plot_grid_panel(ax, grid, colorbar_label, clim, cmap_name)
    imagesc(ax, grid.x, grid.z, grid.values.');
    set(ax, 'YDir', 'normal');

    if strcmpi(cmap_name, 'redblue')
        colormap(ax, blue_white_red(256));
    else
        colormap(ax, cmap_name);
    end

    if ~isempty(clim) && numel(clim) == 2 && all(isfinite(clim)) && clim(2) > clim(1)
        caxis(ax, clim);
    end

    cb = colorbar(ax);
    cb.Label.String = colorbar_label;
    cb.FontSize = 14;

    xlabel(ax, 'X-distance [m]');
    ylabel(ax, 'Depth [m]');
    box(ax, 'on');
end

function overlay_topography_and_acquisition(ax, folder_name, top_zone)
    hold(ax, 'on');

    topo = read_topography(folder_name);
    if ~isempty(topo)
        plot(ax, topo.x, topo.z, '-c', 'LineWidth', 2);
        if top_zone > 0
            fill(ax, [topo.x; flipud(topo.x)], ...
                [topo.z; flipud(topo.z - top_zone)], ...
                'c', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        end
    end

    sr = read_sr_geometry(folder_name);
    if ~isempty(sr)
        scatter(ax, sr.x_rec, sr.z_rec, 18, 'g', 'filled', 'o', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.25, ...
            'MarkerFaceAlpha', 0.75, 'DisplayName', 'Receivers');
        scatter(ax, sr.x_src, sr.z_src, 70, 'r', 'filled', 'p', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.4, ...
            'DisplayName', 'Sources');
    end
end

function apply_map_limits(ax, x_limits, z_limits)
    if ~isempty(x_limits)
        xlim(ax, x_limits);
    end
    if ~isempty(z_limits)
        ylim(ax, z_limits);
    end
end

function set_profile_xlim(ax, profiles)
    vals = [];
    for k = 1:numel(profiles)
        p = profiles{k};
        vals = [vals; p(:)]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
    if isempty(vals)
        return;
    end
    xmin = min(vals);
    xmax = max(vals);
    if xmin == xmax
        pad = 0.15 * max(1, abs(xmin));
    else
        pad = 0.05 * (xmax - xmin);
    end
    xlim(ax, [xmin - pad, xmax + pad]);
end

function save_png(fig, file_name, dpi)
    try
        exportgraphics(fig, file_name, 'Resolution', dpi, 'BackgroundColor', 'white');
    catch
        print(fig, file_name, '-dpng', ['-r', num2str(dpi)]);
    end
    fprintf('Saved: %s\n', file_name);
end

%% ========================================================================
% GENERAL HELPERS
% ========================================================================

function set_default_plot_style()
    set(0, 'DefaultAxesFontName', 'Times New Roman', ...
        'DefaultTextFontName', 'Times New Roman', ...
        'DefaultLegendFontName', 'Times New Roman', ...
        'DefaultColorbarFontName', 'Times New Roman', ...
        'DefaultAxesFontSize', 14, ...
        'DefaultTextFontSize', 16, ...
        'DefaultLegendFontSize', 14, ...
        'DefaultColorbarFontSize', 14, ...
        'DefaultAxesLineWidth', 1.0);
end

function labels = parameter_labels()
    labels = struct();
    labels.rho = 'rho [kg/m^3]';
    labels.C11 = 'C_{11} [Pa]';
    labels.C13 = 'C_{13} [Pa]';
    labels.C33 = 'C_{33} [Pa]';
    labels.C44 = 'C_{44} [Pa]';
    labels.C66 = 'C_{66} [Pa]';
    labels.Q11 = 'Q_{11}';
    labels.Q13 = 'Q_{13}';
    labels.Q33 = 'Q_{33}';
    labels.Q44 = 'Q_{44}';
    labels.Q66 = 'Q_{66}';
    labels.Vp = 'V_p [m/s]';
    labels.Vs = 'V_s [m/s]';
    labels.Qp = 'Q_p';
    labels.Qs = 'Q_s';
end

function label = label_for_param(labels, param)
    if isfield(labels, param)
        label = labels.(param);
    else
        label = param;
    end
end

function tf = should_plot_param(param, selected_params)
    tf = isempty(selected_params) || any(strcmp(param, selected_params));
end

function groups = unique_groups(catalog)
    if isempty(catalog)
        groups = {};
        return;
    end
    all_groups = cell(numel(catalog), 3);
    for k = 1:numel(catalog)
        all_groups{k, 1} = catalog(k).family;
        all_groups{k, 2} = catalog(k).param;
        all_groups{k, 3} = catalog(k).freq;
    end
    [~, ia] = unique(join(string(all_groups), "|", 2), 'stable');
    groups = all_groups(sort(ia), :);
end

function items = select_group(catalog, group_row)
    keep = false(numel(catalog), 1);
    for k = 1:numel(catalog)
        keep(k) = strcmp(catalog(k).family, group_row{1}) && ...
                  strcmp(catalog(k).param, group_row{2}) && ...
                  strcmp(catalog(k).freq, group_row{3});
    end
    items = catalog(keep);
    [~, order] = sort([items.itNum]);
    items = items(order);
end

function grid = make_grid(x, z, values)
    grid.x = x;
    grid.z = z;
    grid.values = values;
end

function lim = finite_minmax(values)
    vals = values(isfinite(values));
    if isempty(vals)
        lim = [];
        return;
    end
    vmin = min(vals);
    vmax = max(vals);
    if vmax <= vmin
        pad = 0.15 * max(1, abs(vmin));
        vmin = vmin - pad;
        vmax = vmax + pad;
    end
    lim = [vmin, vmax];
end

function vm = robust_abs_limit(values, percentile_value, default_value)
    vals = abs(values(isfinite(values)));
    if isempty(vals)
        vm = default_value;
        return;
    end
    vm = prctile(vals, percentile_value);
    if ~isfinite(vm) || vm == 0
        vm = default_value;
    end
end

function M = sid_rid_matrix(sid, rid, values, sids, rids)
    M = nan(numel(sids), numel(rids));
    for i = 1:numel(sids)
        for j = 1:numel(rids)
            keep = sid == sids(i) & rid == rids(j);
            if any(keep)
                M(i, j) = mean(values(keep), 'omitnan');
            end
        end
    end
end

function cmap = blue_white_red(n)
    if nargin < 1
        n = 256;
    end
    n1 = floor(n/2);
    n2 = n - n1;
    blue_to_white = [linspace(0, 1, n1).', linspace(0, 1, n1).', ones(n1, 1)];
    white_to_red = [ones(n2, 1), linspace(1, 0, n2).', linspace(1, 0, n2).'];
    cmap = [blue_to_white; white_to_red];
end

function cmap = frequency_colormap(n)
    try
        cmap = turbo(n);
    catch
        cmap = parula(n);
    end
end

function name = safe_name(text_value)
    name = regexprep(char(text_value), '[^\w.-]', '_');
end

function tf = isabsolute_path(path_value)
    path_value = char(path_value);
    tf = startsWith(path_value, filesep) || ...
         ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once'));
end

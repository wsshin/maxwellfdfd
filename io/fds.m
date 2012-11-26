%% fds
% Advanced interface to the lightlabsFDS solver.

%% Syntax
%  [E, H, err] = fds(omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E, J, max_iters, err_thresh, view_progress)

%% Description
% |[E, H, err] = fds(omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, J, max_iters, err_thresh, view_progress)| 
% solves the time-harmonic Maxwell's equation by communicating with the lightlabsFDS 
% servers. The specific equation solved is
%
% $$ \nabla \times \mu^{-1} \nabla \times E - \omega^2 \epsilon E = -i
% \omega J. $$
%
% Here, the curl operator assumes a default grid spacing of 1. Note also that
% the values of |mu| and |epsilon| are assumed to be absolute (they will not be
% multiplied by any constants).

%% Input parameters
% * |omega| -- angular frequency, complex scalar.
% * |d_prim|, |d_dual| -- length factors for the grid at the primary and dual 
% grid points, respectively. These parameters must be 3-element cell arrays.
% * |s_prim|, |s_dual| -- scaling factors for the grid at the primary and dual 
% grid points, respectively. These parameters must be 3-element cell arrays.
% where every element is a vector with corresponding lengths |xx|, |yy|, and 
% |zz|.
% * |mu, epsilon, E, J| -- permeability, permittivity, initial electric field, and current sources.
% These parameters must by 3-element cell arrays where every element is a 
% 3-dimensional array of size |xx| by |yy| by |zz|. Each array corresponds to
% the x-, y-, or z-component of the vector field.
% * |max_iters| -- maximum number of iterations to run the solver. 
% Must be a positive integer.
% * |err_thresh| -- Threshold error at which the solver terminates. 
% Must be a real positive number.
% * |display_progress| -- option of how to display the progress of the solver.
% Must be either 'plot', 'text', or 'none'.

%% Output parameters
% * |E|, |H| -- electric and magnetic field vectors, respectively. Both of these
% are 3-element cell arrays, and correspond to the x-, y-, and z-components of
% the vector fields respectively.
% * |err| -- the final error in the solution, real positive scalar.

%% Example: Point source with no pml.
%   shape = [100 200 40]; % The shape of the simulation.
%	Jx = zeros(shape);
%	Jx(50, 100, 20) = 1; % Single point source.
%	[E, H, success] = fds(  0.06, ...
%	                        {ones(100,1), ones(200,1), ones(40,1)}, ...
%	                        {ones(100,1), ones(200,1), ones(40,1)}, ...
%	                        {ones(100,1), ones(200,1), ones(40,1)}, ...
%	                        {ones(100,1), ones(200,1), ones(40,1)}, ...
%	                        {ones(shape), ones(shape), ones(shape)}, ...
%	                        {ones(shape), ones(shape), ones(shape)}, ...
%	                        {zeros(shape), zeros(shape), zeros(shape)}, ...
%	                        {Jx, zeros(shape), zeros(shape)}, ...
%	                        5000, 1e-6, 'text');
%
% See additional examples in the |test_fds.m| file.

function [E, H, err] = fds(omega, d_prim, d_dual, s_prim, s_dual, ...
                            mu, epsilon, E, J, ...
                            max_iters, err_thresh, view_progress)

    %% Make sure hdf5 compression is available.

    if ~H5Z.filter_avail('H5Z_FILTER_DEFLATE') || ...
        ~H5ML.get_constant_value('H5Z_FILTER_CONFIG_ENCODE_ENABLED') || ...
        ~H5ML.get_constant_value('H5Z_FILTER_CONFIG_DECODE_ENABLED') || ...
        ~H5Z.get_filter_info('H5Z_FILTER_DEFLATE')
        error('HDF5 gzip filter not available!') 
    end

    %% TODO: Make sure we can communicate with the lightlabsFDS server.

    %% Verify inputs.

    % Check omega.
    if numel(omega) ~= 1 
        error('OMEGA must be a scalar.');
    end

    % Check shapes of mu, epsilon, J.
    % Specifically, make sure all component fields have shape xx-yy-zz.
    shape = size(epsilon{1}); % Make sure all 3D fields have this shape.
    if any([numel(mu), numel(epsilon), numel(E), numel(J)] ~= 3)
        error('All 3D vector fields (MU, EPSILON, E, J) must have three cell elements.');
    end
    fields = [mu, epsilon, E, J];
    for k = 1 : length(fields)
        if any(size(fields{k}) ~= shape)
            error('All 3D vector fields (MU, EPSILON, E, J) must have the same shape.');
        end
    end

    % Check shapes of d_prim, d_dual, s_prim, and s_dual.
    % Specifically, each array of each must have length xx, yy, and zz respectively.
    if any([numel(d_prim), numel(d_dual), numel(s_prim), numel(s_dual)] ~= 3)
        error('D_PRIM, D_DUAL, S_PRIM, and S_DUAL must each have three cell elements.');
    end
    for k = 1 : 3
        d_prim{k} = d_prim{k}(:);
        d_dual{k} = d_dual{k}(:);
        s_prim{k} = s_prim{k}(:);
        s_dual{k} = s_dual{k}(:);
        if (length(d_prim{k}) ~= shape(k)) || (length(d_dual{k}) ~= shape(k) || ...
            length(s_prim{k}) ~= shape(k)) || (length(s_dual{k}) ~= shape(k))
            error('The lengths of D_PRIM, D_DUAL, S_PRIM, and S_DUAL vectors must be xx, yy, and zz, in that order.')
        end
    end

    % Make sure the value for max_iters is valid.
    if mod(max_iters,1) ~= 0 || max_iters <= 0 || ~isreal(max_iters)
        error('MAX_ITERS must be a positive integer.');
    end

    % Make sure the value for err_thresh is valid.
    if ~isfloat(err_thresh) || err_thresh <= 0 || err_thresh >= 1 || ~isreal(err_thresh)
        error('ERR_THRESH must be a positive number between 0 and 1.');
    end

    % Make sure the value for view_progress is valid.
    if all(~strcmp(view_progress, {'plot', 'text', 'none'}))
        error('VIEW_PROGRESS must be either ''plot'', ''text'', or ''none''.');
    end


    %% Construct the exportable hdf5 file.
    % Make all cells use the same write-out function.
        
    % Choose a filename. TODO: Randomize filename to allow for parallelization.
    filename = ['.fds.', strrep(num2str(clock), ' ', ''), '.h5'];
    status_file = [filename, '.status'];
    out_file = [filename, '.out'];

    % Open the hdf5 file, use read-write mode.
    file = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

    % Write omega to the input file. 
    h5write_complex(file, 'omega', omega)

    % Write d_prim, d_dual, s_prim, and s_dual to the input file.
    h5write_field(file, 'd_prim', d_prim);
    h5write_field(file, 'd_dual', d_dual);
    h5write_field(file, 's_prim', s_prim);
    h5write_field(file, 's_dual', s_dual);

    % Write mu, epsilon, E, and J to the input file.
    h5write_field(file, 'mu', mu);
    h5write_field(file, 'epsilon', epsilon);
    h5write_field(file, 'E', E);
    h5write_field(file, 'J', J);

    % Write out max_iters and err_thresh.
    hdf5write(filename, 'max_iters', int64(max_iters), 'WriteMode', 'append');
    hdf5write(filename, 'err_thresh', double(err_thresh), 'WriteMode', 'append');

    H5F.close(file) % Close file, flushing to storage.


    %% Send the file to server.
%     servername = 'raven1:21';
    servername = 'raven1:5988';
    try
        f = ftp(servername, 'fds-user', 'maxwell'); % Connect to ftp server.
        mput(f, filename); % Upload simulation file.
    catch
        error('Could not connect and upload the simulation file.'); 
    end

    delete(filename); % remove file on client-side.


    %% Show simulation progress.

    f_stat = ftp(servername, 'fds-user', 'maxwell'); % Used to get status.
    err_len = 0; % Monitors if more error values have been added or not.

    while true % Keep on updating the status until the simulation completes.

        % Check if the simulation has completed.
        if ~isempty(dir(f_stat, out_file)) % See if the *.out file is present.
            sim_done = true;
        else 
            sim_done = false;
        end

        status = dir(f_stat, status_file);
        if ~isempty(status) % A status file exists.
            mget(f_stat, status_file); % Download it.

            % Attempt to read the status file.
            % This may fail because the status file that we downloaded may
            % have been only partially written to.
            try 
                err = hdf5read(status_file, '/err');
            end

            % If more err values are obtained, refresh progress plot/text.
            if (length(err) > err_len)
                err_len = length(err);
                plot_progress(err, err_thresh, view_progress); 
            end

            delete(status_file); % remove file on client-side.
        end

        if sim_done
            break % Break out of loop, since simulation has completed.
        else
            pause(0.2); % Wait before requerying status.
        end
    end
    close(f_stat); % Close ftp session.


    %% Retrieve final solution.
    
    % Note that we use the initial ftp session here.
    % This actually guarantees that we will only download the fully-written 
    % output file.

    close(f);
%     f = ftp('raven1:21', 'fds-user', 'maxwell'); % Connect to ftp server.
    f = ftp('raven1:5988', 'fds-user', 'maxwell'); % Connect to ftp server.
    mget(f, out_file);
    close(f); % Close ftp session.
    
    % Read out data.
    E = read_from_h5(out_file, 'E');
    H = read_from_h5(out_file, 'H');
    success = hdf5read([filename '.out'], '/success');
    success = success.Data;
    err = hdf5read([filename '.out'], '/err');
    delete(out_file); % remove file on client-side.

    % Plot the final progress.
    plot_progress(err, err_thresh, view_progress, '(final)'); 

    % Warn user if solver did not succeed.
    if ~success
        warning('Did not converge!');
    end


%% Write a complex 3D vector-field to an hdf5 file.
function h5write_field(file, name, data)
    xyz = 'xyz'; 
    for k = 1 : numel(data)
        h5write_complex(file, [name, '_', xyz(k)], data{k});
    end


%% Write a complex-valued 1D or 3D array to an hdf5 file.
function h5write_complex(file, name, data)
 
    % Format data.
    if (ndims(data) == 2) && (size(data,2) == 1) % 1D data detected.
        N = 1;
        data = data;
        dims = numel(data);
        chunk_dims = dims;

    elseif ndims(data) == 3 % 3D data detected.
        N = 3;
        data = permute((data), [ndims(data):-1:1]); % Make data row-major.
        dims = fliplr(size(data)); % Size of the array.
        chunk_dims = [1 dims(2:3)]; % Should heavily affect compression.

    else
        error('Only 1-D or 3-D arrays accepted.')
    end


    % Create the dataspace.
    space = H5S.create_simple(N, dims, []);

    % Set dataspace properties.
    dcpl = H5P.create('H5P_DATASET_CREATE');
    H5P.set_deflate(dcpl, 1); % Deflation level: 0 (none) to 9 (most).
    H5P.set_chunk(dcpl, chunk_dims);

    % Create dataset and write to file.
    dset = H5D.create(file, [name, '_real'], 'H5T_IEEE_F64BE', space, dcpl); % Real part.
    H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', real(data));
    dset = H5D.create(file, [name, '_imag'], 'H5T_IEEE_F64BE', space, dcpl); % Imaginary part.
    H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', imag(data));

    % Close resources.
    H5P.close(dcpl);
    H5D.close(dset);
    H5S.close(space);


%% Extracts 3D field information from an hdf5 file.
function [field] = read_from_h5(file, name)
    xyz = 'xyz';
    for k = 1 : 3
        field{k} = hdf5read(file, ['/', name, '_', xyz(k), '_real']) + ...
                i * hdf5read(file, ['/', name, '_', xyz(k), '_imag']);
        field{k} = permute(field{k}, [ndims(field{k}):-1:1]);
    end


%% Plot the progress of the simulation.
function plot_progress(x, err_thresh, option, varargin)
    if isempty(varargin)
        extra_text = '';
    else
        extra_text = varargin{1};
    end

    if strcmp(option, 'plot')
        % Progress plot.
        semilogy(x, '.-');
        title(extra_text);
        ylabel('relative error');
        xlabel('iteration');
        a = axis;
        hold on
        semilogy([a(1), a(2)], err_thresh * [1 1], 'k--'); % Error threshold line.
        axis([a(1), a(2), err_thresh/10, a(4)])
        hold off

    elseif strcmp(option, 'text')
        fprintf('iter: %d, err: %e %s\n', length(x), x(end), extra_text);
    
    elseif strcmp(option, 'none')
        % Do nothing.

    else
        error('Invalid option.');
    end

    drawnow;


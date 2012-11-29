function write_input(filenamebase, osc, grid3d, s_factor_cell, eps_node_array, eps_cell, mu_cell, J_cell, tol, maxit)

% Check arguments.
chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');
chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GK.count], [1 0]), ...
	'"s_factor_cell" should be %d-by-%d cell array whose each element is row vector with complex elements', Axis.count, GK.count);
for w = Axis.elems
	for g = GK.elems
		chkarg(istypesizeof(s_factor_cell{w,g}, 'complex', [1 grid3d.N(w)]), '"s_factor_cell{%d,%d} should be length-%d.', int(w), int(g), grid3d.N(w));
	end
end
chkarg(istypesizeof(eps_node_array, 'complex', grid3d.N), ...
	'"eps_node_array" should be %d-by-%d-by-%d array with complex elements.', grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(eps_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"eps_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(mu_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"mu_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(J_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"J_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(tol, 'real') && tol > 0, '"tol" should be positive.');
chkarg(istypesizeof(maxit, 'int') && maxit > 0, '"maxit" should be positive integer.');

% Make sure hdf5 compression is available.
if ~H5Z.filter_avail('H5Z_FILTER_DEFLATE') || ...
	~H5ML.get_constant_value('H5Z_FILTER_CONFIG_ENCODE_ENABLED') || ...
	~H5ML.get_constant_value('H5Z_FILTER_CONFIG_DECODE_ENABLED') || ...
	~H5Z.get_filter_info('H5Z_FILTER_DEFLATE')
	error('HDF5 gzip filter is not available.') 
end

%% Construct the exportable hdf5 file.
% Make all cells use the same write-out function.

% Choose a filename. TODO: Randomize filename to allow for parallelization.
filename = ['./', filenamebase, '.h5'];
if exist(filename, 'file')
	delete(filename);
end

% Write simple fields to the input file.

% MATLAB uses column-major arrays, while C uses row-major arrays.  Both the
% types are stored in linearized arrays.  Therefore, to import an array stored
% in MATLAB correctly in C, the fastest-varying index in the C array should be
% the same as the fastest-varying index in the MATLAB array. Since
% BC[Naxis][Nsign] is used in C, the fastest-varying index should be the index
% for signs, so MATLAB arrays should be prepared so that signs vary along the
% column direction.
h5create(filename, '/omega', 1);
h5write(filename, '/omega', double(osc.in_omega0()));

h5create(filename, '/lambda', 1);
h5write(filename, '/lambda', double(osc.in_L0()));

h5create(filename, '/L0', 1);
h5write(filename, '/L0', double(osc.unit.value(PhysQ.L)));

h5create(filename, '/N', Axis.count, 'Datatype', 'int32');
% h5create(filename, '/N', Axis.count);
h5write(filename, '/N', int32(grid3d.N.'));

h5create(filename, '/bc', [Sign.count Axis.count], 'Datatype', 'int32');
% h5create(filename, '/bc', [Sign.count Axis.count]);
h5write(filename, '/bc', int32(subsindex(grid3d.bc.')));  % not int(grid3d.bc.')

h5create(filename, '/Npml', [Sign.count Axis.count], 'Datatype', 'int32');
% h5create(filename, '/Npml', [Sign.count Axis.count]);
h5write(filename, '/Npml', int32(grid3d.Npml.'));

h5create(filename, '/tol', 1);
h5write(filename, '/tol', double(tol));

h5create(filename, '/maxit', 1, 'Datatype', 'int32');
% h5create(filename, '/maxit', 1);
h5write(filename, '/maxit', int32(maxit));

% Write complex values.
e_ikL = exp(-1i * (grid3d.kBloch .* grid3d.L));
h5create(filename, '/e_ikL', [2 Axis.count]);
h5write(filename, '/e_ikL', expand_complex(e_ikL.'));

for w = Axis.elems
	datasetname = ['/', char(w), '_prim'];
	real_array = grid3d.lall{w,GK.prim}.';
	h5create(filename, datasetname, length(real_array));
	h5write(filename, datasetname, real_array);
	
	for g = GK.elems
		g_str = char(g);
		datasetname = ['/d', char(w)', '_', g_str(1:4)];
		complex_array = grid3d.dl{w,g}.';
		h5create(filename, datasetname, [2 length(complex_array)]);
		h5write(filename, datasetname, expand_complex(complex_array));

		datasetname = ['/s', char(w)', '_', g_str(1:4)];
		complex_array = s_factor_cell{w,g}.';
		h5create(filename, datasetname, [2 length(complex_array)]);
		h5write(filename, datasetname, expand_complex(complex_array));
	end
end

use_petsc = true;
if use_petsc
	eps_array = cell2array(eps_cell, Axis.count);
	if isreal(eps_array)
		eps_array = complex(eps_array);
	end
% 	PetscBinaryWrite([filenamebase, '.eps'], eps_array(:), 'indices', 'int32', 'precision', 'float64');
	PetscBinaryWrite([filenamebase, '.eps'], eps_array(:));

	J_array = cell2array(J_cell, Axis.count);
	if isreal(J_array)
		J_array = complex(J_array);
	end
% 	PetscBinaryWrite([filenamebase, '.J'], 1i*J_array(:), 'indices', 'int32', 'precision', 'float64');
	PetscBinaryWrite([filenamebase, '.J'], J_array(:));

% 	J_array2 = PetscBinaryRead([filenamebase, '.J'], 'complex', true);
% 	norm(1i*J_array(:)-J_array2)
% 	eps_array2 = PetscBinaryRead([filenamebase, '.eps'], 'complex', true);
% 	eps_array2 = PetscBinaryRead([filenamebase, '.eps'], true);
% 	norm(eps_array-eps_array2);
else
	cl = 5;  % compression level (0-9, where 0 means no compression)
	dims = [2 Axis.count grid3d.N];

	% h5create(filename, '/eps_node', dims, 'Deflate', cl, 'ChunkSize', dims);
	% h5write(filename, '/eps_node', ...
	% 	expand_complex(cell2array({eps_node_array, eps_node_array, eps_node_array}, Axis.count)));

	h5create(filename, '/eps', dims, 'Deflate', cl, 'ChunkSize', dims);
% 	h5create(filename, '/eps', dims);
	h5write(filename, '/eps', expand_complex(cell2array(eps_cell, Axis.count)));

	% h5create(filename, '/mu', dims, 'Deflate', cl, 'ChunkSize', dims);
	% h5write(filename, '/mu', expand_complex(cell2array(mu_cell, Axis.count)));

	h5create(filename, '/J', dims, 'Deflate', cl, 'ChunkSize', dims);
% 	h5create(filename, '/J', dims);
	h5write(filename, '/J', expand_complex(cell2array(J_cell, Axis.count)));
end

% % Write complex arrays to the input file.
% file = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
% for w = Axis.elems
% 	for g = GK.elems
% 		g_str = char(g);
% 		h5write_complex_array1d(file, ['d', char(w)', '_', g_str(1:4)], grid3d.dl{w,g});
% 		h5write_complex_array1d(file, ['s', char(w)', '_', g_str(1:4)], s_factor_cell{w,g});
% 	end
% end
% 
% N = grid3d.N;
% h5write_complex_array3d(file, 'mu', mu_cell, grid3d.N);
% h5write_complex_array3d(file, 'eps', eps_cell, N);
% h5write_complex_array3d(file, 'J', J_cell, N);
% 
% H5F.close(file) % Close file, flushing to storage.



% %% Write a complex-valued 1D or 3D array to an hdf5 file.
% function h5write_complex_array1d(file, datasetname, data_row)
% chkarg(istypesizeof(datasetname, 'char', [1 0]), '"datasetname" should be string.');
% chkarg(istypesizeof(data_row, 'complex', [1 0]), '"data_row" should be row vector with complex elements.');
% 
% data = [real(data_row); imag(data_row)];
% dims = size(data);
% nD = ndims(data);
% data = permute(data, nD:-1:1);
% dims = wrev(dims);
% 
% % Create the dataspace.
% space = H5S.create_simple(nD, dims, []);
% 
% % Set dataspace properties.
% dcpl = H5P.create('H5P_DATASET_CREATE');
% H5P.set_deflate(dcpl, 5); % Deflation level: 0 (none) to 9 (most).
% dims_chunk = dims;
% % dims_chunk(1) = 1;
% H5P.set_chunk(dcpl, dims_chunk);
% 
% % Create dataset and write to file.
% dset = H5D.create(file, datasetname, 'H5T_IEEE_F64LE', space, dcpl);
% H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
% 
% % Close resources.
% H5P.close(dcpl);
% H5D.close(dset);
% H5S.close(space);
% 
% 
% %% Write a complex-valued 1D or 3D array to an hdf5 file.
% function h5write_complex_array3d(file, datasetname, data_cell, N3d)
% chkarg(istypesizeof(datasetname, 'char', [1 0]), '"datasetname" should be string.');
% chkarg(istypesizeof(N3d, 'int', [1 Axis.count]), '"N3d" should be length-%d row vector.', Axis.count);
% chkarg(istypesizeof(data_cell, 'complexcell', [1 Axis.count], N3d), ...
% 	'"data_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
% 	Axis.count, N3d(Axis.x), N3d(Axis.y), N3d(Axis.z)); 
% 
% data = [data_cell{Axis.x}(:), data_cell{Axis.y}(:), data_cell{Axis.z}(:)];
% data = data.';
% data = data(:);
% data = [real(data), imag(data)];
% dims = [2, Axis.count, N3d];
% data = reshape(data, dims);
% n_dims = ndims(data);
% data = permute(data, n_dims:-1:1);
% dims = wrev(dims);
% 
% 
% % Create the dataspace.
% space = H5S.create_simple(n_dims, dims, []);
% 
% % Set dataspace properties.
% dcpl = H5P.create('H5P_DATASET_CREATE');
% H5P.set_deflate(dcpl, 5); % Deflation level: 0 (none) to 9 (most).
% dims_chunk = dims;
% % dims_chunk(1) = 1;
% H5P.set_chunk(dcpl, dims_chunk);
% 
% % Create dataset and write to file.
% dset = H5D.create(file, datasetname, 'H5T_IEEE_F64LE', space, dcpl);
% H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
% 
% % Close resources.
% H5P.close(dcpl);
% H5D.close(dset);
% H5S.close(space);

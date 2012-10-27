function eps_face_cell = harmonic_mean_eps_node(eps_node)

chkarg(istypesizeof(eps_node, 'complex', zeros(1, Axis.count)), ...
	'"eps_node" should be %dD array with complex elements.', Axis.count);

eps_face_cell = cell(1, Axis.count);
eps_face_cell{Axis.x} = 2./(1./eps_node(1:end-1, 1:end-1, 1:end-1) + 1./eps_node(2:end, 1:end-1, 1:end-1));
eps_face_cell{Axis.y} = 2./(1./eps_node(1:end-1, 1:end-1, 1:end-1) + 1./eps_node(1:end-1, 2:end, 1:end-1));
eps_face_cell{Axis.z} = 2./(1./eps_node(1:end-1, 1:end-1, 1:end-1) + 1./eps_node(1:end-1, 1:end-1, 2:end));

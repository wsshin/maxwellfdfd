function node_array = expand_node_array(grid3d, node_array)

% Extend eps and mu to the ghost points considering the boundary conditions.
for w = Axis.elems
	ind_gn = {':',':',':'};
	ind_gp = {':',':',':'};	
	if grid3d.bc(w) == BC.p
		ind_gn{w} = grid3d.N(w);
		ind_gp{w} = 1;
	else
		ind_gn{w} = 1;
		ind_gp{w} = grid3d.N(w);
	end
	node_array = cat(int(w), node_array(ind_gn{:}), node_array, node_array(ind_gp{:}));
end

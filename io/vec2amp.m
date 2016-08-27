function Famp = vec2amp(F_cell)

chkarg(istypesizeof(F_cell, 'Scalar3dcell', [1 Axis.count]), ...
	'"E_cell should be length-%d cell array with Scalar3d as elements"', Axis.count);

Fx = F_cell{Axis.x};
grid3d = Fx.grid3d;
osc = Fx.osc;
physQcell = Fx.physQcell;
name = Fx.name(1:end-2);

li = grid3d.l(:,GT.dual);
[Xi, Yi, Zi] = ndgrid(li{:});

fi = NaN([grid3d.N, Axis.count]);
for w = Axis.elems
	Fw = F_cell{w};
	[fw, l] = Fw.data_ghost_expanded();
	[X, Y, Z] = ndgrid(l{:});
	fi(:,:,:,w) = interpn(X, Y, Z, fw, Xi, Yi, Zi);	
end

famp = sqrt(sum(abs(fi).^2, Axis.count+1));
assert(all(size(famp)==grid3d.N));

% Attach extra points.
for w = Axis.elems
	famp = attach_extra_Famp(famp, w, grid3d);
end

gt_array = [GT.dual GT.dual GT.dual];

Famp = Scalar3d(famp, grid3d, gt_array, osc, physQcell, ['|', name, '|']);


function array = attach_extra_Famp(array, w, grid3d)
% This makes Painter2d.prep_data() estimates the field values at the boundary of
% the computational domain inaccurate.  For example, suppose that we have Ez on
% the xy-plane with a PEC boundary at the x-boundaries.  Then, Ez = 0 at x = 0,
% but for Eamp we have nonzero |Ez| before and after x = 0, leading to a nonzero
% estimate of |Ez| at x = 0.  We cannot handle this case correctly in the
% current visualization framework...  It might be better to estimate Eamp and
% the corners of Yee's cells rather than centers?

ind_n = {':', ':', ':'};
ind_p = {':', ':', ':'};
bc_w = grid3d.bc(w);
if bc_w == BC.p
	ind_n{w} = grid3d.N(w);
	ind_p{w} = 1;
else  % bc_d == BC.e or BC.m
	ind_n{w} = 1;
	ind_p{w} = grid3d.N(w);	
end
array = cat(int(w), array(ind_n{:}), array, array(ind_p{:}));  % Bloch phases in S are ignored in amplitude

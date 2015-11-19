%% MaxwellFDFD Documentation
%
% <html>
% <script src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript">
% </script>
% </html>

%% Introduction
% MaxwellFDFD is a MATLAB-based solver package of Maxwell's equations.  It
% solves the equations by the finite-difference frequency-domain (FDFD) method,
% and hence the name MaxwellFDFD.  The resulting solution completely describes
% the the interaction of electromagnetic (EM) waves with objects of interest,
% such as waveguides, antennas, photonic crystals, and photovoltaic cells.
%
% The FDFD method solves the frequency-domain form of Maxwell's equations, which
% are
%
% \(
% \nabla \times \mathbf{E}(\mathbf{r}) = -i \, \omega \, \mu(\mathbf{r},\mathbf{\omega}) \, \mathbf{H}(\mathbf{r}) - \mathbf{M}(\mathbf{r}),\\
% \nabla \times \mathbf{H}(\mathbf{r}) = i \, \omega \, \varepsilon(\mathbf{r},\mathbf{\omega}) \, \mathbf{E}(\mathbf{r}) + \mathbf{J}(\mathbf{r}),
% \)
%
% where 
%
% * \(\mathbf{E}(\mathbf{r})\) and \(\mathbf{H}(\mathbf{r})\) are the solution
% electric and magnetic fields of the EM waves that we want to obtain.
% * \(\mathbf{J}(\mathbf{r})\) and \(\mathbf{M}(\mathbf{r})\) are the electric
% and magnetic current source densities at a position \(\mathbf{r}\) that
% emanate EM waves.
% * \(\omega\) is the oscillation frequency of the current sources and EM waves.
% * \(\varepsilon(\mathbf{r},\omega)\) and \(\mu(\mathbf{r},\omega)\) are the
% electric permittivity and magnetic permeability of the object at a position
% \(\mathbf{r}\).
%
% To construct the frequency-domain Maxwell's equations, you first need to
% specify the frequency \(\omega\), and then you need to place objects and
% current sources in your simulation domain.

%% Examples
% See various problems that MaxwellFDFD can solve in <topic/gallery2d.html 2D
% Example Gallery> and <topic/gallery3d.html 3D Example Gallery>.

%% Main Features
% * Built-in frequency-dependent dielectric constants
% (\(\varepsilon(\omega)/\varepsilon_0\)) for commonly used nanophotonic
% materials (e.g., \(\mathrm{Si}\), \(\mathrm{SiO_2}\), \(\mathrm{Ag}\),
% \(\mathrm{Au}\)) from trusted references such as Palik, Johnson and Christy,
% and CRC Handbook.
% * Various object shapes: box, sphere, cylinder, etc.
% * Various types of sources: point source, line source, plane wave source
% (with an oblique emission support), etc.
% * Total-field/scattered-field (TF/SF) method.
% * PEC, PMC, and periodic (or Bloch) boundary condition.
% * Perfectly matched layer (PML) absorbing boundary.
% * Power flux calculation.
% * Direct solver (for small problems) and iterative solver (for large
% problems).
% * Visualization of objects and sources.
% * Visualization of solution fields.

%% User Guides
% * <topic/basic.html Basic Usage>
% * <topic/intermediate.html Intermediate Usage (for building complicated systems)>
% * <topic/fd3d.html Using FD3D for Large Problems>

%% User Guides (to be written)
% * <topic/material.html Assigning Material Parameters>
% * <topic/shape.html Placing Objects and Sources>
% * <topic/vis.html Visualizing Solutions>
% * <topic/flux.html Measuring Power Flux>
% * <topic/tfsf.html Total-Field/Scattered-Field Method>
% * <topic/bc.html Using Boundary Conditions as Symmetry Planes>
% * <topic/diel.html Adding Dielectric Constants>
% * <topic/pml.html Changing PML Parameters>
% * <topic/mathematica.html Exporting Solutions for Mathematica>
% * <topic/trouble.html Troubleshooting>

%% List of Functions and Classes
% * <comp/components.html List of Functions and Classes>


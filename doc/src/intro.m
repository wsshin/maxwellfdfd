%% Introduction to MaxwellFDFD
%
% <html>
% <script src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript">
% </script>
% </html>
%
% MaxwellFDFD is a MATLAB-based solver package of the frequency-domain Maxwell's
% equations.  By solving the equations, we obtain the solution that completely
% describes the the interaction of electromagnetic (EM) waves with EM objects of
% interest, such as waveguides, antennas, photonic crystals, and photovoltaic
% cells.
%
% The frequency-domain Maxwell's equations are differential equations described
% as
%
% \(
% \nabla \times \mathbf{E}(\mathbf{r}) = -i \, \omega \, \mu(\mathbf{r},\mathbf{\omega}) \, \mathbf{H}(\mathbf{r}) - \mathbf{M}(\mathbf{r}),\\
% \nabla \times \mathbf{H}(\mathbf{r}) = i \, \omega \, \varepsilon(\mathbf{r},\mathbf{\omega}) \, \mathbf{E}(\mathbf{r}) + \mathbf{J}(\mathbf{r}),
% \)
%
% where 
%
% * \(\mathbf{E}(\mathbf{r})\) and \(\mathbf{H}(\mathbf{r})\) are the solution
% electric and magnetic fields of the EM waves that we want to obtain,
% * \(\mathbf{J}(\mathbf{r})\) and \(\mathbf{M}(\mathbf{r})\) are the electric
% and magnetic current source densities at a position \(\mathbf{r}\) that
% emanate EM waves,
% * \(\omega\) is the oscillation frequency of the current sources,
% * \(\varepsilon(\mathbf{r},\omega)\) and \(\mu(\mathbf{r},\omega)\) are the
% electric permittivity and magnetic permeability of the EM object at a position
% \(\mathbf{r}\).
%
% To construct the frequency-domain Maxwell's equations, you first need to
% decide the frequency \(\omega\) you are interested in, and then you need to
% place EM objects, which are described in terms of shapes and material
% parameters \(\varepsilon(\mathbf{r},\omega)\) and \(\mu(\mathbf{r},\omega)\),
% in your simulation domain.  Here
%


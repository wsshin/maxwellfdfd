MaxwellFDFD installation instruction
====================================
System requirments
------------------
- Operating System

	Windows, LINUX, OS X, or any operating system where MATLAB can be installed.

- MATLAB

	Version 2010a or later.  (Older versions are not supported due to the lack of some object-oriented programming features that are heavily used in MaxwellFDFD.)  No Toolbox is needed.


How to install MaxwellFDFD
--------------------------
1. Extract the package tarball (or zip file); it will create a directory `maxwellfdfd/`.

2. Move the `maxwellfdfd/` to wherever you want to store the package.

3. In the MATLAB startup directory, open `startup.m` file.  (Create the file if there is not.)  

	On Mac and Windows, the default path to the startup directory is `~/Documents/MATLAB/`, but the exact location depends on your setting.  To locate the startup directory, execute `userpath` in MATLAB Command Window; it returs a list of directories separated by semicolons (:), and the first directory in the list is the startup directory.

4. In `startup.m` file opened in Step 3, copy and paste the following lines:

		maxwellfdfd_root = 'INSTALLATION DIRECTORY'
		addpath(maxwellfdfd_root);
		addpath([maxwellfdfd_root, filesep, 'base']);
		addpath([maxwellfdfd_root, filesep, 'diff']);
		addpath([maxwellfdfd_root, filesep, 'dielconst']);
		addpath([maxwellfdfd_root, filesep, 'grid']);
		addpath([maxwellfdfd_root, filesep, 'integ']);
		addpath([maxwellfdfd_root, filesep, 'io']);
		addpath([maxwellfdfd_root, filesep, 'material']);
		addpath([maxwellfdfd_root, filesep, 'modesolver']);
		addpath([maxwellfdfd_root, filesep, 'petsc']);
		addpath([maxwellfdfd_root, filesep, 'shape']);
		addpath([maxwellfdfd_root, filesep, 'source']);
		addpath([maxwellfdfd_root, filesep, 'vis']);

	and change `INSTALLATION DIRECTORY` in the first line to the directory path you chose in Step 2.  (Make sure that the directory path is enclosed by single quotes.)

5. Restart MATLAB.  

	Alternatively, you can execute the above created startup.m script in MATLAB Command Window.

6. Test MaxwellFDFD.

	In MATLAB, open the installed `maxwellfdfd/` directory, open the subdirectory `example/2d/`.  There are 2D examples.  Try to run the simplest one, e.g., `pointsrc_2d.m`.  If it runs without errors and shows a field plot, the program is installed correctly.

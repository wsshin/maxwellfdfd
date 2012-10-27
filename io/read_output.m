function [E_cell, H_cell] = read_output(filenamebase)
chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');

E_array = h5read([filenamebase, '.E.h5'], '/E');


function 

function S = HDMMF_load(name)
%HDMMF_load Load the reference PS/PF data shipped with the HDMMF suite.
%   S = HDMMF_load('NMMF1') loads PF / PS / draw_pf from
%   NMMF1_Reference_PSPF_data.mat located NEXT TO THIS FILE.
%
%   PORT NOTE (2026-09-27): the original problems used a relative-path
%   load('NMMF1_Reference_PSPF_data.mat'), which only worked when the
%   current directory happened to be the problem folder - it broke the
%   platform GUI reference-set display and any batch run from another
%   directory. This helper resolves the absolute path via mfilename, so
%   the data is found regardless of the working directory.

    f = fullfile(fileparts(mfilename('fullpath')), ...
                 [name '_Reference_PSPF_data.mat']);
    S = load(f);
end

function blkStruct = slblocks
% Loads core VPC library into SIMULINK.
% This function specifies that the library should appear in the Library
% Browser and be cached in the browser repository.

% Specify library to load
Browser.Library = 'simlib_vpc'; % filename
Browser.Name = 'VPC Library';   % name that appears in SIMULINK lib-browser

blkStruct.Browser = Browser;
end
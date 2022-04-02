function sl_customization(cm)
% Change order of loaded libraries in the Simulink Library Browser.

% Set VPC library of higher order than Simulink
cm.LibraryBrowserCustomizer.applyOrder({'VPC Library',-2});

end
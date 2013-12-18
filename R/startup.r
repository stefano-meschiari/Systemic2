# This file is executed automatically by the GUI. Add any commands that you want to run at the beginning of a session here.

cat(sprintf('\nHello master %s!\nThis is Console %.4f speaking\nStefano Meschiari, >>>http://www.stefanom.org\nPlease cite >>>http://adsabs.harvard.edu/abs/2009PASP..121.1016M when you use this software.\n?* Type quickstart() to get started.\n?* Type acknowledgments() to view a list of acknowledgments.\n?* Type tutorial() to view a brief tutorial.\n\n',
	Sys.getenv("USER"), SYSTEMIC.VERSION))

if (file.exists("../private/startup.r"))
	source("../private/startup.r")

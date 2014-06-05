# This file is executed automatically by the GUI. Add any commands that you want to run at the beginning of a session here.

# Default options

# Autosave kernels in the autosave/ directory every 240 seconds
options(systemic.autosave = TRUE)
options(systemic.autosave.interval = 240)

# Change the palette here (see R/colors.r for themes)
systemic.plot.theme('tomorrow')

# Set the default graphic device


# This is the startup message. You can make it wittier if you'd like.
cat(sprintf('Hello master %s!\nThis is Console %.4f speaking, %s.\n(c) 2013-2014 Stefano Meschiari, >>>http://www.stefanom.org\nPlease cite >>>http://adsabs.harvard.edu/abs/2009PASP..121.1016M when you use this software.\n?* Type quickstart() to get started.\n?* Type tutorial() to view a brief tutorial.\n?* Type acknowledgments() to view the license and a list of acknowledgments.\n',
	Sys.getenv("USER"), SYSTEMIC.VERSION, R.version.string))

# Window title
.gui.title(sprintf("Systemic %.4f", SYSTEMIC.VERSION))

# Set the default CRAN repository
options(repos=c(CRAN="http://cran.us.r-project.org"))

if (file.exists("../private/startup.r"))
	source("../private/startup.r")

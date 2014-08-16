# Tomorrow palette
# https://github.com/chriskempson/tomorrow-theme
systemic.theme.tomorrow <- c(
    '#000000', # black
    '#4271ae', # blue
    '#c82829', # red
    '#718c00', # green
    '#8959a8', # purple
    '#3e999f', # aqua
    '#f5871f', # orange
    '#eab700', # yellow
    '#8f5536', # brown
    '#d6d6d6', # dark Gray
    '#d33682'  # bubblegum
)

systemic.theme.tomorrow.face <- adjustcolor(systemic.theme.tomorrow, 0.75)

# Default R palette
systemic.theme.r <- palette("default")
systemic.theme.r.face <- palette("default")

systemic.palette <- systemic.theme.tomorrow
systemic.palette.face <- systemic.theme.tomorrow.face
